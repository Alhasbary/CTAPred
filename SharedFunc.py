"""
Compound-Target Activity Prediction (CTAPred) Tool

This script serves as the core framework for the CTAPred tool, which is designed to facilitate the prediction of protein 
targets for natural products. 

Key Features:
    1. **System Requirements Check**: Ensures the necessary versions of dependencies (e.g., RDKit, MolVS, DuckDB) are 
       available, preventing compatibility issues.
    2. **Helper Functions**: Provides utility functions that assist with logging, validation, and general data handling tasks.
    3. **Update Dataset Versions**: Handles the processing and standardization of new dataset versions (e.g., ChEMBL, COCONUT), 
       transforming them into compressed Parquet files for efficient storage and future use.
    4. **Create CTA Reference Dataset**: Generates and manages the CTA reference dataset, a crucial component that links 
       chemical compounds to their bioactivity data, enabling prediction of their potential protein targets.
    5. **Target Prediction**: Leverages the CTA reference dataset along with similarity-based models to predict protein targets for novel compounds.
    6. **Parallel Processing**: Optimizes tasks such as SMILES string preprocessing through parallel execution, significantly 
       speeding up data transformation and predictions.
    7. **Efficient Data Storage**: Saves processed datasets in Parquet format using efficient compression algorithms to reduce 
       storage size while maintaining accessibility for future analyses.

Prepared By: Alhasbary
Date: 10/10/2024
"""

import argparse
import glob
import time
import os
import gc

import numpy as np
import pandas as pd

import sqlite3
import duckdb

from joblib import Parallel, delayed

import molvs
import rdkit
from rdkit import Chem
from rdkit.Avalon.pyAvalonTools import GetAvalonFP
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


import fastparquet
import joblib
import tqdm

###############################################  Check System Requirements  ###############################################
# This section verifies the system's installed versions of essential libraries before execution. It ensures that the 
# required versions of RDKit, MolVS, DuckDB, Fastparquet, NumPy, pandas, Joblib, tqdm, and SQLite3 are present.
#
# The required versions below reflect the minimum needed for specific functionality to run correctly. If the installed
# versions of these libraries are older than the specified versions, some functions may behave unexpectedly, or the
# script may fail entirely due to API changes or missing features in older versions.
#
# If newer versions of these dependencies introduce breaking changes, update this script to reflect those changes and
# adjust the version requirements accordingly.
#
# Ensure the following versions or newer are installed:
#   - RDKit: >= 2020.09.1 (important for molecule processing functions)
#   - MolVS: >= 0.1.1 (used for molecular standardization)
#   - DuckDB: >= 1.1.1 (required for efficient querying with parquet files)
#   - Fastparquet: >= 0.5.0 (for writing and reading parquet files)
#   - NumPy: >= 1.21.5 (needed for numerical operations)
#   - Pandas: >= 1.3.5 (used for handling tabular data)
#   - Joblib: >= 1.1.0 (essential for parallel processing)
#   - tqdm: >= 4.66.4 (for progress bars and tracking execution progress)
#   - SQLite3: >= 2.6.0 (used for interacting with the SQLite database)
###########################################################################################################################

rdkversion = rdkit.__version__.split(".")
if rdkversion < ["2020", "09", "1"]:
    raise ValueError("need an RDKit version >= 2020.09.1")

molvsversion = molvs.__version__.split(".")
if molvsversion < ["0", "1", "1"]:
    raise ValueError("need a MolVS version >= 0.1.1")

duckdbversion = duckdb.__version__.split(".")
if duckdbversion < ["1", "1", "1"]:
    raise ValueError("need a DuckDB version >= 1.1.1")

fastparquetversion = fastparquet.__version__.split(".")
if fastparquetversion < ["0", "5", "0"]:
    raise ValueError("need a Fastparquet version >= 0.5.0")

numpyversion = np.__version__.split(".")
if numpyversion < ["1", "21", "5"]:
    raise ValueError("need a NumPy version >= 1.21.5")

pandasversion = pd.__version__.split(".")
if pandasversion < ["1", "3", "5"]:
    raise ValueError("need a pandas version >= 1.3.5")

joblibversion = joblib.__version__.split(".")
if joblibversion < ["1", "1", "0"]:
    raise ValueError("need a Joblib version >= 1.1.0")

tqdmversion = tqdm.__version__.split(".")
if tqdmversion < ["4", "66", "4"]:
    raise ValueError("need a tqdm version >= 4.66.4")

sqlite3version = sqlite3.version.split(".")
if sqlite3version < ["2", "6", "0"]:
    raise ValueError("need a SQLite3 version >= 2.6.0")


##########################################       Helper Functions        ##################################################
# This section contains utility functions that assist with various tasks required throughout the script.
###########################################################################################################################
                
def assert_query_resource_files(log, args):
    """
    Make sure the existence of the query list files
    """

    args.allQueryLists = glob.glob(args.inputPath+"/*_smiles.csv")
    if len(args.allQueryLists)==0:
        print("Cannot find any query list file!", file=log)
        print("The query list file must be a CSV file. The tool can processes more than one query list.", file=log)
        print("File names should follow the pattern: QueryListN_smiles.csv, where N is an integer.", file=log)
        print("Each CSV file must contain a 'smiles' column for SMILES strings and a 'np_id' column for compound IDs.", file=log)            
        return False
       
    args.allQueryLists = sorted(args.allQueryLists)
    args.datNames = []
    for p in args.allQueryLists:
        args.datNames.append(os.path.basename(p).split("_")[0])
    return True
                        
def tstamp():
    t = time.localtime()
    return "%d/%d/%d %d:%.2d.%d" % ((t[1], t[2], t[0]) + t[3:6])

class Tee(object):
    def __init__(self, *args):
        self.files = [f for f in args]

    def write(self, s):
        for f in self.files:
            f.write(s)

    def flush(self):
        for f in self.files:
            f.flush()

def int_range(min_value, max_value):
    def int_range_type(value):
        ivalue = int(value)
        if ivalue < min_value or ivalue > max_value:
            raise argparse.ArgumentTypeError("Value must be between %s and %s. For the Avalon fingerprint, the value "
                                             "is recommended to be divisible by 8." % (min_value, max_value))
        return ivalue

    return int_range_type

def float_range(min_value, max_value):
    def float_range_type(value):
        fvalue = float(value)
        if fvalue < min_value or fvalue > max_value:
            raise argparse.ArgumentTypeError("Value must be between %s and %s" % (min_value, max_value))
        return fvalue
    return float_range_type

##########################################      Update Dataset Versions        #############################################
# This section handles the processing and transformation of updated versions of datasets, 
# such as ChEMBL and COCONUT. It standardizes the data into efficient, compressed Parquet files for storage.
#
# Key functionalities include:
#   - Preprocessing the data to ensure consistency and remove any invalid or problematic entries.
#   - Converting the cleaned data into a Parquet format, optimized for both storage and retrieval.
#
# These steps ensure that the latest versions of the datasets are properly formatted and ready for use in analysis.
############################################################################################################################

def refine_chembl_database(log, dataPath, destinationPath, compression='brotli'):
    """
    Refine the ChEMBL database and save as parquet files.
    
    This function processes the ChEMBL SQLite database and extracts key tables, 
    transforming them into compressed Parquet files. 
    
    :param log: Log file to record progress.
    :param dataPath: Full filepath of the SQLite file containing the ChEMBL database.
    :param destinationPath: Full filepath for saving the Parquet files.
    :param compression: Compression algorithm for Parquet files. Default is 'brotli'.
    :return: Parquet files stored in the 'refined_chembl' folder.
    """    
    try:
        time_s = time.time()
        # Connect to the ChEMBL SQLite database
        with sqlite3.connect(dataPath) as conn:
            
            # List of tables and corresponding SQL queries
            tables_sql = {
                'activities': '''SELECT activity_id, assay_id, molregno, standard_relation, standard_type, 
                                 standard_value, standard_units, activity_comment, data_validity_comment, 
                                 potential_duplicate, pchembl_value 
                                 FROM activities 
                                 WHERE standard_value IS NOT NULL 
                                 AND pchembl_value IS NOT NULL 
                                 AND standard_units = 'nM' 
                                 AND standard_relation = '=' 
                                 AND potential_duplicate = 0 
                                 AND standard_type IN ('IC50', 'Ki', 'Kd', 'EC50') 
                                 AND (data_validity_comment IS NULL OR data_validity_comment = 'manually validated') 
                                 AND (activity_comment IS NULL OR LOWER(activity_comment) NOT IN 
                                     ('not active', 'inactive', 'inconclusive', 'undetermined')) 
                                 AND standard_value <= 10000''',
                'target_dictionary': '''SELECT tid, target_type, pref_name, organism, chembl_id FROM target_dictionary
                                        WHERE LOWER(target_type) IN ('single protein', 'protein complex')''',
                'assays': '''SELECT assay_id, assay_type, description, tid, confidence_score FROM assays 
                             WHERE assay_type = 'B' ''',
                'molecule_dictionary': 'SELECT molregno, chembl_id FROM molecule_dictionary',
                'compound_structures': 'SELECT molregno, canonical_smiles FROM compound_structures WHERE canonical_smiles IS NOT NULL',
                'compound_properties': 'SELECT molregno, np_likeness_score FROM compound_properties',
                'component_sequences': 'SELECT component_id, accession, sequence FROM component_sequences WHERE accession IS NOT NULL',
                'protein_classification': 'SELECT protein_class_id, protein_class_desc FROM protein_classification',
                'component_class': 'SELECT component_id, protein_class_id, comp_class_id FROM component_class',
                'target_components': 'SELECT targcomp_id, tid, component_id FROM target_components',
            }

            # Loop through each table, execute the query, and write the result to Parquet
            print(f"Reading ChEMBL data from {dataPath}", file=log)
            for table, sql in tables_sql.items():
                parquet_file_path = os.path.join(destinationPath, f"{table}.parquet")

                print(f"Processing table: {table}", file=log)
                df = pd.read_sql_query(sql, conn)
                
                # To focus on relevant compounds instead of processing all compounds
                if table == 'activities':
                    relevant_compounds = df.molregno.unique()
                    

                # Additional processing for compound_structures table
                if table == 'compound_structures':
                    print(f"Processing ChEMBL data for SMILES cleaning...", file=log)
                    print(f"Cleaning {len(relevant_compounds)} SMILES out of {len(df.molregno.unique())} SMILES found in ChEMBL.", file=log)
                    df = df[df['molregno'].isin(relevant_compounds)].reset_index(drop=True)
                    cleaned_smi = Parallel(n_jobs=-1, backend="multiprocessing")(
                        delayed(processMS)(row['molregno'], row['canonical_smiles'], verbose=args.verbose)  
                        for _, row in tqdm.tqdm(df.iterrows(), total=len(df), desc="ChEMBL processMS")
                    )

                    clean_smiles = pd.DataFrame(cleaned_smi, columns=['molregno', 'clean_smiles'])
                    df['clean_smiles'] = clean_smiles['clean_smiles']

                    # Remove problematic compounds
                    df = df[df['clean_smiles'] != 'issue'].reset_index(drop=True)
                    df.dropna(inplace=True)
                    print(f"{len(df.molregno.unique())} SMILES have bean successfully cleaned.", file=log)
                    
                # Save to Parquet
                print(f"Transforming {table} into compressed Parquet file and Saving it to {parquet_file_path} ...", file=log)
                df.to_parquet(parquet_file_path, compression=compression)
            
            print("Refined-ChEMBL Parquet database created successfully!", file=log)

    except sqlite3.DatabaseError as db_err:
        print(f"Database error occurred: {db_err}", file=log)
    except IOError as io_err:
        print(f"I/O error occurred: {io_err}", file=log)
    except Exception as e:
        print(f"An error occurred: {e}", file=log)
    finally:
        elapsed_time = time.time() - time_s
        print(f"ChEMBL file processing completed in {elapsed_time:.2f} seconds.", file=log)
   
def refine_NPs_data(log, dataPath, destinationPath, verbose, compression='brotli'):
    """
    Refine the NPs data and save as parquet files.
    
    This function processes the NPs data (e.g., COCONUT SMILES) and preprocesses them,
    transforming the data into compressed parquet files. 

    :param log: Log file to record progress.
    :param dataPath: Full filepath of the NPs file containing the COCONUT SMILES.
    :param destinationPath: Full filepath for saving the Parquet files.
    :param compression: Compression algorithm for Parquet files. Default is 'brotli'.
    :return: Parquet files stored in the 'destinationPath' folder.
    """
    try:
        time_s = time.time()
        # Read the CSV file containing NPs data
        print(f"Reading NPs data from {dataPath}", file=log)
        NPs = pd.read_csv(dataPath, usecols=['identifier', 'canonical_smiles'])   
        NPs.dropna(inplace=True)     
        print(f"NPs data contains {len(NPs.identifier.unique())} SMILES.", file=log)
        #NPs= NPs[:10]

        # Process each row to clean SMILES using multiprocessing
        print(f"Processing NPs data for SMILES cleaning...", file=log)
        cleaned_smi = Parallel(n_jobs=-1, backend="multiprocessing")(
            delayed(processMS)(row['identifier'], row['canonical_smiles'], verbose=verbose)  
            for _, row in tqdm.tqdm(NPs.iterrows(), total=len(NPs), desc="NPs processMS")
        )

        # Create DataFrame of cleaned SMILES
        clean_smiles = pd.DataFrame(cleaned_smi, columns=['molregno', 'clean_smiles'])
        NPs['clean_smiles'] = clean_smiles['clean_smiles']
        NPs['molregno'] = clean_smiles['molregno']
        
        # Remove problematic compounds (where clean_smiles is 'issue') and drop NaN values
        NPs = NPs[NPs['clean_smiles'] != 'issue'].reset_index(drop=True)
        NPs.dropna(inplace=True)
        print(f"{len(NPs.molregno.unique())} SMILES have bean successfully cleaned.", file=log)

        # Save the cleaned NPs data to Parquet
        parquet_file_path = os.path.join(destinationPath, "NPs_clean_smiles.parquet")        
        print(f"Transforming cleaned NPs SMILES into compressed Parquet file and Saving it to {parquet_file_path} ...", file=log)
        NPs.to_parquet(parquet_file_path, compression=compression)
        
    except FileNotFoundError as fnf_error:
        print(f"Error: File not found - {fnf_error}", file=log)
    except Exception as e:
        print(f"An error occurred: {e}", file=log)
    finally:
        elapsed_time = time.time() - time_s
        print(f"NPs file processing completed in {elapsed_time:.2f} seconds.", file=log)
      
def initialize_CTA_files(log, destinationPath, compression='brotli'):
    """
    Initialize the CTA metadata and data files, setting up the necessary structure 
    and data types, and saving them as parquet files. This function ensures the required CTA files are created 
    only once if they do not already exist.

    :param log : Log file to record progress.  
    :param destinationPath: Path where the CTA files will be saved. Default is 'refined_chembl'.
    :param compression: Compression method for Parquet files. Default is 'brotli'.
    """    
    # Define file paths for CTA_meta and CTA_data Parquet files
    parquet_meta_path = os.path.join(destinationPath, "CTA_meta.parquet")
    parquet_data_path = os.path.join(destinationPath, "CTA_data.parquet")

    try:
        # Define appropriate data types for CTA_meta
        CTA_meta_data = {
            'dataset_id': pd.Series([], dtype='int64'),
            'dataset': pd.Series([], dtype='object'),
            'fingerprint': pd.Series([], dtype='object'),
            'nBits': pd.Series([], dtype='int64'),
            'radius': pd.Series([], dtype='int64'),
            'Tc': pd.Series([], dtype='float64'),
            'standard_value': pd.Series([], dtype='float64'),
            'targets': pd.Series([], dtype='int64'),
            'compounds': pd.Series([], dtype='int64'),
            'records': pd.Series([], dtype='int64')
        }
        # Create DataFrame for CTA_meta with the defined data types
        CTA_meta = pd.DataFrame(CTA_meta_data)
        CTA_meta.to_parquet(parquet_meta_path, compression=compression)
        print(f"Saved initialized CTA_meta to {parquet_meta_path}", file=log)

        # Define appropriate data types for CTA_data
        CTA_data_data = {
            'dataset_id': pd.Series([], dtype='int64'),
            'uniprot': pd.Series([], dtype='object'),
            'target_chembl_id': pd.Series([], dtype='object'),
            'molregno': pd.Series([], dtype='object')
        }
        # Create DataFrame for CTA_data with the defined data types
        CTA_data = pd.DataFrame(CTA_data_data)
        CTA_data.to_parquet(parquet_data_path, compression=compression)
        print(f"Saved initialized CTA_data to {parquet_data_path}", file=log)

    except Exception as e:
        print(f"An error occurred: {e}", file=log)

def set_ChEMBL_as_first_option_in_CTA_table(log, destinationPath, compression='brotli'):
    """
    Set the new version of ChEMBL as the first option in the CTA reference dataset table.
    
    :param log : Log file to record progress.  
    :param destinationPath: Path where the CTA files will be saved. Default is 'refined_chembl'.
    :param compression: Compression method for Parquet files. Default is 'brotli'.
    """    
    # Retrieve the ChEMBL dataset
    ChEMBL = extractDB(log)    
    if not ChEMBL.empty:
        results_df = ChEMBL[['molregno', 'uniprot', 'target_chembl_id']].drop_duplicates(ignore_index=True)
        targets = len(results_df.uniprot.unique())
        compounds = len(results_df.molregno.unique())
        records = len(results_df)
        
        # Load the CTA tables from the Parquet files
        CTA_meta_file_path = f"{destinationPath}/CTA_meta.parquet"
        CTA_data_file_path = f"{destinationPath}/CTA_data.parquet"
        CTA_meta = pd.read_parquet(CTA_meta_file_path)
        CTA_data = pd.read_parquet(CTA_data_file_path) 
        
        # Remove existing entries with dataset_id == 1 (ChEMBL)
        if not CTA_meta.empty:
            CTA_meta = CTA_meta[CTA_meta.dataset_id != 1]
            CTA_data = CTA_data[CTA_data.dataset_id != 1]
        
        # Prepare a new row for the CTA_meta DataFrame
        new_dataset_id = 1  # ChEMBL will always be assigned ID 1
        new_row = pd.DataFrame({
            'dataset_id': [new_dataset_id],
            'dataset': ['ChEMBL'],
            'fingerprint': [None],
            'nBits': [None],
            'radius': [None],
            'Tc': [None],
            'standard_value': [10000],
            'targets': [targets],
            'compounds': [compounds],
            'records': [records]
        })
        
        # Append the new row to the CTA_meta DataFrame
        CTA_meta = pd.concat([CTA_meta, new_row], ignore_index=True)

        # Save the updated CTA_meta DataFrame
        CTA_meta.to_parquet(CTA_meta_file_path, compression=compression)
        
        # Append new ChEMBL data to the CTA_data DataFrame
        results_df['dataset_id'] = new_dataset_id
        CTA_data = pd.concat([CTA_data, results_df], ignore_index=True)

        # Save the updated CTA_data DataFrame
        CTA_data.to_parquet(CTA_data_file_path, compression=compression)
        
##########################################      Create CTA Reference Dataset        ############################################
# This section focuses on the creation of the Compound-Target Activity (CTA) reference dataset, 
# which serves as a key resource for target prediction models. 
#
# Key tasks include:
#   - Aggregating and preprocessing compound-target activity data from various sources.
#   - Standardizing and filtering the data to ensure high-quality associations.
#   - Storing the processed data in Parquet format for optimized access and usage in modeling tasks.
#
# The CTA reference dataset is a crucial component in the development and evaluation of target prediction algorithms..
###############################################################################################################################

def extractDB(log, standard_value=10000):
    """
    Extract High-Quality Data using DuckDB and Parquet files.
    :param log: Log file for recording extraction information.
    :param standard_value: Activity values in nM measurement.
    :return: DataFrame containing the selected data.
    """
    try:
        # Define the path to the Parquet folder
        parquet_folder = "Data/refined_chembl"

        # Create a connection to DuckDB
        con = duckdb.connect()

        # Construct the query using parquet_scan for each relevant table
        df = con.query(f"""
            SELECT DISTINCT
                cs.molregno, md.chembl_id AS molecule_chembl_id, 
                cs.clean_smiles, act.standard_value, 
                act.standard_type, act.pchembl_value, ass.confidence_score, 
                ass.description, td.organism, td.target_type, td.pref_name, 
                td.chembl_id AS target_chembl_id, cseq.accession AS uniprot, 
                pc.protein_class_desc, cp.np_likeness_score
            FROM 
                parquet_scan('{os.path.join(parquet_folder, 'compound_structures.parquet')}') cs
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'molecule_dictionary.parquet')}') md
            ON 
                cs.molregno = md.molregno
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'compound_properties.parquet')}') cp
            ON 
                cp.molregno = md.molregno
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'activities.parquet')}') act
            ON 
                md.molregno = act.molregno
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'assays.parquet')}') ass
            ON 
                act.assay_id = ass.assay_id
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'target_dictionary.parquet')}') td
            ON 
                ass.tid = td.tid
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'target_components.parquet')}') tc
            ON 
                tc.tid = td.tid
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'component_sequences.parquet')}') cseq
            ON 
                cseq.component_id = tc.component_id
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'component_class.parquet')}') cclass
            ON 
                cseq.component_id = cclass.component_id
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'protein_classification.parquet')}') pc
            ON 
                cclass.protein_class_id = pc.protein_class_id
            WHERE 
                ass.confidence_score IN (7, 9)
            AND 
                act.standard_value <= {standard_value}
        """).df()

        # Close the DuckDB connection
        con.close()

        # Remove duplicate rows except for 'protein_class_desc'
        selected_bioactivity = df.drop_duplicates(subset=df.columns.difference(['protein_class_desc']))

        return selected_bioactivity

    except Exception as e:
        print(f"An error occurred: {e}")
        exit(0)

def retrieveCTA(dataset_id):

    """
    Extract High-Quality Data using DuckDB and Parquet files.
    :param dataset_id: CTA ID to be retrieved.
    :return: DataFrame containing the selected data.
    """
    try:
        # Define the path to the Parquet folder
        parquet_folder = "Data/refined_chembl"

        # Create a connection to DuckDB
        con = duckdb.connect()

        # Construct the query using parquet_scan for each relevant table
        df = con.query(f"""
            SELECT DISTINCT
                cs.molregno, md.chembl_id AS molecule_chembl_id, 
                cs.clean_smiles, act.standard_value, 
                act.standard_type, act.pchembl_value, ass.confidence_score, 
                ass.description, td.organism, td.target_type, td.pref_name, 
                td.chembl_id AS target_chembl_id, cseq.accession AS uniprot, 
                pc.protein_class_desc, cp.np_likeness_score
            FROM 
                parquet_scan('{os.path.join(parquet_folder, 'compound_structures.parquet')}') cs
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'molecule_dictionary.parquet')}') md
            ON 
                cs.molregno = md.molregno
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'compound_properties.parquet')}') cp
            ON 
                cp.molregno = md.molregno
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'activities.parquet')}') act
            ON 
                md.molregno = act.molregno
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'assays.parquet')}') ass
            ON 
                act.assay_id = ass.assay_id
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'target_dictionary.parquet')}') td
            ON 
                ass.tid = td.tid
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'target_components.parquet')}') tc
            ON 
                tc.tid = td.tid
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'component_sequences.parquet')}') cseq
            ON 
                cseq.component_id = tc.component_id
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'component_class.parquet')}') cclass
            ON 
                cseq.component_id = cclass.component_id
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'protein_classification.parquet')}') pc
            ON 
                cclass.protein_class_id = pc.protein_class_id
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'CTA_data.parquet')}') CTA_data
            ON 
                td.chembl_id = CTA_data.target_chembl_id
            JOIN 
                parquet_scan('{os.path.join(parquet_folder, 'CTA_meta.parquet')}') CTA_meta
            ON 
                CTA_data.dataset_id = CTA_meta.dataset_id
            WHERE 
                ass.confidence_score IN (7, 9)
            AND 
                cs.molregno = CTA_data.molregno
            AND 
                CTA_meta.dataset_id = {dataset_id}
        """).df()
        
        # Close the DuckDB connection
        con.close()

        # Remove duplicate rows except for 'protein_class_desc'
        CTA = df.drop_duplicates(subset=['molregno', 'uniprot'])

        return CTA

    except Exception as e:
        print(f"An error occurred: {e}")
        exit(0)
         
def processMS(smi_id, smi, verbose=False):
    """
    Process Molecular Structures with the standardization tool.
    
    This function standardizes a given SMILES string by performing several 
    preprocessing steps. Then, the final SMILES string is generated with stereochemistry preserved.

    :param smi_id: Identifier for the SMILES string (e.g., molecule ID).
    :param smi: SMILES string to be processed.
    :param verbose: Whether to print error messages (default is True).
    
    :return: Tuple containing SMILES ID and the standardized SMILES string. 
             If any issues occur during processing, 'issue' is returned as the SMILES string.
    """
    
    # Initialize the result with default 'issue' values
    smiles = 'issue'
    standard = molvs.Standardizer(prefer_organic=True)
    
    try:
        mol = Chem.MolFromSmiles(smi)        
        if mol is not None:
            # Cleane and standardize using MolVS
            std_mol = standard.fragment_parent(mol, skip_standardize=False)            
            # Generate a SMILES string with stereochemistry
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        else:
            if verbose:            
                print(f"\nError: Invalid SMILES for molecule {smi_id}")
    except Exception as e:
        if verbose: 
            print(f"\nError processing molecule {smi_id}: {e}")
    
    # Return the SMILES ID and the standardized SMILES string (or 'issue' if any step failed)
    return smi_id, smiles     

                
def gen_fp(smi, fingerprint='ecfp', nBits=2024, radius=2, verbose=False):
    """
    Generate a molecular fingerprint from a SMILES string using RDKit.

    Parameters:
    smi (str): SMILES string representing the molecule.
    fingerprint (str): Type of fingerprint to generate ('ecfp', 'fcfp', 'maccs', or 'avalon').
    nBits (int): Length of the generated fingerprint (default: 2024).
    radius (int): Radius for ECFP/FCFP fingerprints (default: 2).
    verbose (bool): If True, print error messages (default: True).

    Returns:
    str: Generated fingerprint or 'issue' if an error occurs.
    """
    
    # Default result in case of an error
    fp = 'issue'
    
    try:
        # Convert SMILES to RDKit Mol object with sanitization
        mol = Chem.MolFromSmiles(smi, sanitize=True)
        if mol:         
            try: 
                # Generate the selected fingerprint type
                if fingerprint == 'avalon':
                    fp = GetAvalonFP(mol, nBits=nBits)
                elif fingerprint == 'ecfp':
                    fp = GetMorganFingerprintAsBitVect(mol, nBits=nBits, radius=radius)
                elif fingerprint == 'fcfp':
                    fp = GetMorganFingerprintAsBitVect(mol, nBits=nBits, radius=radius, useFeatures=True)
                elif fingerprint == 'maccs':
                    fp = MACCSkeys.GenMACCSKeys(mol)
            except Exception as e:
                if verbose:
                    print(f"Error generating {fingerprint} fingerprint for SMILES {smi}: {e}")
        else:
            if verbose:
                print(f"Error: Invalid SMILES {smi}")
    except Exception as e:
        if verbose:
            print(f"Error processing SMILES {smi}: {e}")
    
    # Return the fingerprint or 'issue'
    return fp


def gen_fps(log, args, df, db_name):
    """
    Generate molecular fingerprints for the given DataFrame of SMILES strings.

    Parameters:
    log (file object): Log file for recording progress.
    args (Namespace): Arguments containing fingerprint type, number of jobs, and fingerprint parameters.
    df (DataFrame): DataFrame containing the molecular data with columns 'clean_smiles' and 'molregno'.

    Returns:
    DataFrame: DataFrame containing the generated fingerprints and corresponding 'molregno'.
    """
              
    # Remove duplicate SMILES and reset index
    df = df.drop_duplicates(subset=['molregno', 'clean_smiles'], ignore_index=True)
    
    print(f'Generating fingerprints for {len(df)} unique {db_name} compounds ...', file=log)
    start_time = time.time()
    
    # Generate fingerprints using parallel processing
    fps = Parallel(n_jobs=int(args.n_jobs), backend="multiprocessing")(
        delayed(gen_fp)(row['clean_smiles'], 
                        fingerprint=args.fingerprint, nBits=args.nBits, radius=args.radius, verbose=args.verbose) 
        for _, row in tqdm.tqdm(df.iterrows(), total=len(df), desc=f"{db_name} fingerprints"))
    
    # Create a DataFrame for fingerprints
    df_fps = pd.DataFrame({'fp': fps, 'molregno': df['molregno']})
    
    # Filter out rows with 'issue' in the fingerprint
    df_fps = df_fps[df_fps['fp'] != 'issue'].reset_index(drop=True)
    
    elapsed_time = time.time() - start_time
    print(f'Fingerprint generation completed in {elapsed_time:.4f} seconds.', file=log)
    
    return df_fps
      
def createCTA(log, args, output_dir):
    """
    Generate the CTA reference dataset
    :param log: Log file for recording processing information.
    :param args: argparse object containing program arguments.
    """
    try:    
        # Define the path to the Parquet folder
        parquet_folder = "Data/refined_chembl"

        # Create a connection to DuckDB
        con = duckdb.connect()

        # Construct the query using parquet_scan for each relevant table
        CTA_meta = con.query(f"""
            SELECT DISTINCT
                dataset_id, dataset, fingerprint , nBits, radius, 
                Tc, standard_value, targets, compounds, records
            FROM 
                parquet_scan('{os.path.join(parquet_folder, 'CTA_meta.parquet')}') """).df()

        # Close the DuckDB connection
        con.close()                

        found = False
        for _, row in CTA_meta.iterrows():
            if (row['standard_value'] == args.standard_value) and \
            (row['fingerprint'] == args.fingerprint) and \
            (row['nBits'] == args.nBits) and \
            (row['Tc'] == args.Tc): 
                found = True
                dataset_id = int(row['dataset_id'])
                break

        if not found:    
            if len(CTA_meta) > 0:
                print('\nNote: There are existing reference datasets with the following parameters:\n')
                CTA_meta.sort_values(by=['dataset_id'],inplace=True)
                print(CTA_meta.to_string(index=False))                
                inp = input('\nWould you like to create a new dataset with the current parameters? ([y]/n): ').strip().lower() or 'y'
                
            else:
                inp = 'y'
          
        else:            
            inp = input('\nA reference dataset with the specified parameters already exists. \nWould you like to export it to a CSV file? ([y]/n): ').strip().lower() or 'y'
            if inp == 'y':
                CTA = retrieveCTA(dataset_id)                
                
                # Save to csv file 
                fname = os.path.join(output_dir, f'CTA_fp_{args.fingerprint}_nBits_{args.nBits}_activity_{args.standard_value}_Tc_{args.Tc}.csv')
                CTA.to_csv(fname, index=False, float_format='%.4f')   
                print(f'\nSaved in {fname}')            
                exit(0)
            
        # Handle the user's choice
        if inp == 'y':
            # User chose to create a new dataset with the current parameters
            print('\nProceeding with the new parameter options to create a CTA dataset. \nRetrieving ChEMBL (less than 60 sec)...', file=log)
            
            # Preprocess all ChEMBL SMILES generate fingerprints by RDKit.
            clean_ChEMBL_structures = extractDB(log, args.standard_value)  
            #clean_ChEMBL_structures = clean_ChEMBL_structures[:20000]
            ChEMBL = clean_ChEMBL_structures[['molregno', 'clean_smiles', 'uniprot', 'target_chembl_id']].drop_duplicates(ignore_index=True)
            ChEMBL.dropna(inplace=True)                   
            ChEMBL_fps = gen_fps(log, args, ChEMBL[['molregno', 'clean_smiles']], 'ChEMBL').drop_duplicates(ignore_index=True)
            ChEMBL = pd.merge(ChEMBL, ChEMBL_fps, on='molregno', how='inner')
            ChEMBL.dropna(inplace=True)                   
                    
            # Preprocess all NPs SMILES generate fingerprints by RDKit. 
            parquet_file_path = "Data/NPs_clean_smiles.parquet"
            clean_NPs_structures = pd.read_parquet(parquet_file_path)
            #clean_NPs_structures = clean_NPs_structures[:10000]
            clean_NPs_structures.dropna(inplace=True)  
            NPs_fps = gen_fps(log, args, clean_NPs_structures, 'NPs').drop_duplicates(ignore_index=True)
                                        
            # Extracting unique UniProt identifiers for ChEMBL targets interacting with similar compounds
            results_df, targets, compounds, records, elapsed_time = CTA_conductSS(log, args, NPs_fps, ChEMBL, args.Tc)                          
            
            if not results_df.empty:
                # Update CTA_meta table  
               
                # Load the CTA_meta DataFrame from the Parquet file
                parquet_file_path = "Data/refined_chembl/CTA_meta.parquet"
                CTA_meta = pd.read_parquet(parquet_file_path)

                # Calculate the new dataset_id
                if not CTA_meta.empty:
                    new_dataset_id = CTA_meta['dataset_id'].max() + 1
                else:
                    new_dataset_id = 2  # If the dataframe is empty, start dataset_id from 2. ChEMBL will always be assigned ID 1                

                # Create a new row as a DataFrame
                new_row = pd.DataFrame({
                    'dataset_id': [new_dataset_id],
                    'dataset': [f'CTA_{new_dataset_id}'],
                    'fingerprint': [args.fingerprint],
                    'nBits': [args.nBits],
                    'radius': [args.radius],
                    'Tc': [args.Tc],
                    'standard_value': [args.standard_value],
                    'targets': [targets],
                    'compounds': [compounds],
                    'records': [records]
                })

                # Append the new row to the CTA_meta DataFrame
                CTA_meta = pd.concat([CTA_meta, new_row], ignore_index=True)

                # Save the updated DataFrame back to the Parquet file
                CTA_meta.to_parquet(parquet_file_path, compression='brotli')

                # Update CTA_data table    
                parquet_file_path = "Data/refined_chembl/CTA_data.parquet"
                CTA_data = pd.read_parquet(parquet_file_path) 
                results_df['dataset_id'] =  new_dataset_id
                CTA_data = pd.concat([CTA_data, results_df], ignore_index=True)
                CTA_data.to_parquet(parquet_file_path,compression='brotli')
                
                CTA = retrieveCTA(int(new_dataset_id))
                print('\nThe CTA tables have been successfully updated with the new reference dataset.')     
                    
        else:
            # User decided not to proceed
            print('\nProcess aborted by the user.', file=log)
    
    except NameError as e:
        print(f"An error occurred: {e}")
        exit(0)        
                  
def select_CTA_reference_dataset(log, args):
    """
    Select the CTA reference dataset
    :param log: Log file for recording processing information.
    :param args: argparse object containing program arguments.
    """
    try:    
        # Define the path to the Parquet folder
        parquet_folder = "Data/refined_chembl"

        # Create a connection to DuckDB
        con = duckdb.connect()

        # Construct the query using parquet_scan for each relevant table
        CTA_meta = con.query(f"""
            SELECT DISTINCT
                dataset_id, dataset, fingerprint , nBits, radius,  
                Tc, standard_value, targets, compounds, records
            FROM 
                parquet_scan('{os.path.join(parquet_folder, 'CTA_meta.parquet')}') """).df()

        # Close the DuckDB connection
        con.close()
                
        CTA = None     # Initialize CTA variable

        if len(CTA_meta) > 0:
            print('Available CTA datasets were created with the following parameters:\n')
            CTA_meta.sort_values(by=['dataset_id'],inplace=True)
            print(CTA_meta.to_string(index=False))                
            inp = input('\nWould you like to select one of the existing datasets? Enter the dataset_id to proceed, or press [Enter] to abort: ').strip().lower()
            
            if inp.isdigit():
                dataset_id = int(inp)
                print(f'\nUser chose an existing dataset by entering dataset_id = {dataset_id}. \nRetrieving the reference dataset (less than 60 sec)...', file=log)
                CTA = retrieveCTA(dataset_id)
            else:
                print('\nProcess aborted by the user.', file=log)
                exit(0)
        else:
            print('\nThere are no CTA datasets yet. Please use "generate_CTA.py" to generate one.\n')
            exit(0)

    except Exception as e:
        print(f"An error occurred: {e}")
        exit(0)  
        
    return CTA, dataset_id    

    
def find_similar_targets(args, row, reference_fps, Tc):
    """
    Find targets with Tanimoto similarity above a given threshold (Tc) for a given fingerprint.

    Parameters:
    args (Namespace): Arguments containing relevant configurations.
    row (Series): Row containing the fingerprint to compare.
    reference_fps (DataFrame): DataFrame containing reference fingerprints and associated targets.
    Tc (float): Tanimoto coefficient threshold for similarity.

    Returns:
    list: List of tuples containing the 'uniprot' and 'target_chembl_id' of similar targets.
    """
    # Calculate Tanimoto similarity between the given fingerprint and reference fingerprints
    ss_results = pd.Series(Chem.DataStructs.BulkTanimotoSimilarity(row['fp'], list(reference_fps['fp'])), index=reference_fps.index)
    
    # Filter results based on the similarity threshold Tc
    similar_targets = reference_fps.loc[ss_results[ss_results >= Tc].index, ['molregno', 'uniprot', 'target_chembl_id']].apply(tuple, axis=1).tolist()
    
    return similar_targets
   
    
def CTA_conductSS(log, args, query_fps, reference_fps, Tc):
    """
    Perform chemical similarity searches and identify target interactions based on a Tanimoto coefficient (Tc).

    Parameters:
    log (file object): Log file to record progress and results.
    args (Namespace): Arguments containing configurations such as fingerprint type, nBits, radius, and more.
    query_fps (DataFrame): DataFrame containing query fingerprints.
    reference_fps (DataFrame): DataFrame containing reference fingerprints and target information.
    Tc (float): Tanimoto coefficient threshold for similarity.

    Returns:
    list: A list of unique UniProt IDs corresponding to similar targets.
    """

    # Description of the similarity search process
    desc = f"""
    Performing chemical similarity searches with {args.nBits}-bit {args.fingerprint} fingerprints 
    and a radius of {args.radius}. Targets that interact with compounds having a Tc of {args.Tc} or higher 
    will be considered in the CTA. The search involves {len(reference_fps.uniprot.unique())} unique ChEMBL targets.
    """
    print(desc, file=log)

    # Remove duplicate fingerprints
    query_fps = query_fps.drop_duplicates()
    
    # Start timing the process
    start_time = time.time()
    
    # Perform similarity search in parallel
    results = Parallel(n_jobs=args.n_jobs, backend="multiprocessing")(
        delayed(find_similar_targets)(args, row, reference_fps, Tc) for _, row in tqdm.tqdm(query_fps.iterrows(), total=len(query_fps), desc="CTA conductSS")
    )
    
    # Flatten the results list since each row might return multiple similar targets
    results_flat = [item for sublist in results for item in sublist]

    # Convert to DataFrame
    results_df = pd.DataFrame(results_flat, columns=['molregno', 'uniprot', 'target_chembl_id']).drop_duplicates(subset=['molregno', 'uniprot'], ignore_index=True)

    targets = len(results_df.uniprot.unique())
    compounds = len(results_df.molregno.unique())
    records = len(results_df)
    
    
    # Calculate elapsed time
    elapsed_time = time.time() - start_time
    
    # Log the results and completion status
    if not results_df.empty:
        print(f"Similarity search completed in {elapsed_time:.2f} seconds", file=log)
        print(f"The dataset comprised {records} bioactivity records, covering {compounds} unique compounds and {targets} unique targets", file=log)
    else:
        print(f"No targets identified with the given similarity threshold. Process completed in {elapsed_time:.2f} seconds.", file=log)

    return results_df, targets, compounds, records, elapsed_time

    
##########################################     Predict Potential Protein Targets     #################################################
# This section focuses on the identification of potential protein targets for compounds using the CTA reference dataset.
#
# Key tasks include:
#   - Utilizing similarity-based or machine learning models to predict likely protein targets.
#   - Leveraging the curated CTA dataset to find associations between compounds and protein targets.
#   - Outputting ranked predictions for further validation and analysis.
#
# The predictions derived from this section are instrumental in expanding knowledge of compound-target interactions and 
# supporting drug discovery initiatives.
#######################################################################################################################################
        
def preprocess_Query_SMILES(log, args, query_list):
    """
    Process Query list Molecular Structures.
    :param log: Log file for recording processing information.
    :param args: argparse object containing program arguments.
    :param query_list: DataFrame containing the query list data.
    :return: DataFrame containing the cleaned query list data.
    """
    try:
        time_s = time.time()    
        
        print('Number of compounds in the query list file:', query_list.shape[0], file=log)
          
                    
        cleaned_smi = Parallel(n_jobs=int(args.n_jobs), backend="multiprocessing")(
            delayed(processMS)(row['np_id'], row['smiles'], verbose=args.verbose) for _, row in tqdm.tqdm(query_list.iterrows(), total=len(query_list), desc="Query processMS"))               
       
        clean_smiles = pd.DataFrame(cleaned_smi, columns=['np_id', 'clean_smiles'])                    
        # Remove rows containing problematic compounds
        clean_smiles = clean_smiles[clean_smiles['clean_smiles'] != 'issue'].reset_index(drop=True)    
        clean_smiles.dropna(inplace=True)        
        
        print('SMILES in the query list after preprocessing:', clean_smiles.shape[0], file=log)
        
        elapsed_time = time.time() - time_s
        print(f'Preprocess query SMILES strings completed in {elapsed_time:.2f} seconds.', file=log)   
        
        return clean_smiles
        
    except NameError as e:
        print(f"An error occurred: {e}")
        exit(0)


def Query_conductSS(log, args, query_fps, reference_fps):
    desc = f'''\nPerform chemical similarity searches. 
    The cleaned SMILES strings for each compound are converted into {args.nBits}-bit {args.fingerprint} fingerprints 
    with a radius of {args.radius}.'''    
    
    print(desc, file=log)
    query_fps = query_fps.drop_duplicates()
    
    print(f'\nStart chemical similarity searches for {len(query_fps)} smiles\n', file=log)
    time_s = time.time()

    # Initialize lists to store relevant information
    np_ids = []
    molregnos = []
    target_chembl_ids = []
    uniprots = []
    scores = []
    
    # Group reference_fps by 'uniprot'
    grouped_fps = reference_fps.groupby('uniprot')
    
    if isinstance(args.top_k, list):
        k = max(args.top_k)
    else:        
        k = args.top_k
    
    # Iterate over each query
    for idx, row in tqdm.tqdm(query_fps.iterrows(), total=len(query_fps), desc="Query conductSS"):
        
        # Calculate similarity scores for each group (target)
        for target, fps in grouped_fps:
            # Calculate similarity scores
            ss_results = pd.Series(Chem.DataStructs.BulkTanimotoSimilarity(row['fp'], list(fps['fp'])), index=fps.index)
            
            # Get indices of top N similarity scores
            top_indices = ss_results.nlargest(k, keep='all').index
            
            # Store relevant information
            np_ids.extend([row['np_id']] * len(top_indices))
            molregnos.extend(fps.loc[top_indices, 'molregno'].tolist())
            target_chembl_ids.extend(fps.loc[top_indices, 'target_chembl_id'].tolist())
            uniprots.extend(fps.loc[top_indices, 'uniprot'].tolist())
            scores.extend(ss_results[top_indices].tolist())

    # Create DataFrame with relevant information
    df_result = pd.DataFrame({
        'np_id': np_ids,
        'molregno': molregnos,
        'target_chembl_id': target_chembl_ids,
        'uniprot': uniprots,
        'score': scores
    })

    if not df_result.empty:
        identified_targets = len(df_result.uniprot.unique())
        elapsed_time = time.time() - time_s
        print(f'Similarity search completed in {elapsed_time:.2f} seconds. {identified_targets} unique potential targets identified."', file=log)     
    else:
        print('No similarities were found!', file=log)

    return df_result


def rankTargets(log, args, reference_dataset, dest):
    """
    Retrieving ranked potential targets.
    :param log: Log file for recording search information.
    :param args: argparse object containing program arguments.
    :param reference_dataset: DataFrame containing the results of conductSS.
    :param dest: File name for saving the potential targets for each SMILES string.
    """ 
          
    desc = f'\nRanking the potential targets for each SMILES string.\n'
    print(desc, file=log)
    time_s = time.time()    
    
    if reference_dataset.empty:
        print(f'The provided CTA reference dataset is empty!', file=log)
        return 
    
        
    # Group reference_dataset by ['uniprot', 'np_id']
    grouped = reference_dataset.groupby(['uniprot', 'np_id'])
    
    if isinstance(args.top_k, list):
        k_values = args.top_k
    else:        
        k_values = [args.top_k]
    
    print(f"k_values = {k_values}")

    # Loop over each k-value
    for top_k in k_values :
        # Initialize lists to store relevant information for each k-value
        np_ids = []
        target_chembl_ids = []
        uniprots = []
        mean_scores = []

        # Iterate over each group (['uniprot', 'np_id'])
        for _, group_df in tqdm.tqdm(grouped, desc=f"k = {top_k}", total=len(grouped)):
            # Select the top K scores within the group       
           
            if top_k == 1:
                # Special handling for k=1, just take the first (top) score
                top_K_group = group_df[['score']].max()  
            else:
                # For k > 1, use nlargest to get the top K scores, keeping ties
                top_K_group = group_df.nlargest(top_k, 'score', keep='all')[['score']]
                
            #top_K_group = group_df.nlargest(top_k, 'score', keep='all')

            # Calculate the mean of the top K scores for the group
            mean_score = top_K_group['score'].mean()

            # Append relevant information to lists
            np_ids.append(group_df['np_id'].iloc[0])
            target_chembl_ids.append(group_df['target_chembl_id'].iloc[0])
            uniprots.append(group_df['uniprot'].iloc[0])
            mean_scores.append(mean_score)

        # Create DataFrame with relevant information for the current k-value
        df = pd.DataFrame({
            'np_id': np_ids,
            'target_chembl_id': target_chembl_ids,
            'uniprot': uniprots,
            'mean_score': mean_scores,
        })

        # Calculate ranks based on the mean score, handling each 'np_id' separately
        df['rank'] = df.groupby('np_id')['mean_score'].rank(ascending=False, method='min').astype(int)

        # Sort values by np_id and mean_score
        df.sort_values(by=['np_id', 'mean_score'], inplace=True, ascending=[True, False])

        # Save the ranking to a CSV file
        df.to_csv(f'{dest}_{top_k}.csv', index=False, float_format='%.4f')

    print(f'\nResults:', file=log)
    print(f'\tIdentified {len(reference_dataset.uniprot.unique())} potential targets for {len(reference_dataset.np_id.unique())} SMILES strings.', file=log)
    print(f"\tAdditional details are stored in {dest}_{k_values}.csv", file=log)

    print(f'Ranking the potential targets took {time.time() - time_s:.4f} seconds.', file=log)

