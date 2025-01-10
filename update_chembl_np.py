"""
Compound-Target Activity Prediction (CTAPred) Tool - Dataset Update Script

    This script is part of the CTAPred tool. It manages the processing and transformation of updated versions of key datasets 
    such as ChEMBL and COCONUT (natural products database). The data is standardized into compressed Parquet files for efficient 
    storage and use in target prediction models. The script allows for updating one or both datasets as required, ensuring they 
    are prepared for further analysis and stored in the 'refined_chembl' folder.

Inputs:
   - --ChEMBL_data : Full path to the ChEMBL SQLite database file [Required].
   - --NPs_data    : Full path to the Natural Products (NPs) data CSV file [Required].
   - --destination : Directory path for saving the refined Parquet files [Optional, default='Data'].
   - --ChEMBL      : Flag to update the ChEMBL database.
   - --NPs         : Flag to update the Natural Products data.

Outputs: 
   - 'refined_chembl' folder in the specified destination directory, containing the refined ChEMBL Parquet files.
   - 'NPs_clean_smiles.parquet' file in the specified destination directory, containing the refined NPs data.

Usage:
   - For help:
     python update_versions_ChEMBL_NPs.py -h
   - To run the program:
     python update_versions_ChEMBL_NPs.py --ChEMBL_data=FullPathToChEMBL --NPs_data=FullPathToNPs --destination=FullPathToSaveParquet

Usage Example:
    1. Download the ChEMBL SQLite dataset from: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/
       Extract using: tar -xvzf chembl_34_sqlite.tar.gz 

    2. Download the COCONUT dataset from: https://coconut.naturalproducts.net/download/
       Extract using: tar -xvzf coconut-10-2024.csv.zip

    3. To update both datasets:
       python update_versions_ChEMBL_NPs.py --ChEMBL_data=chembl_34/chembl_34_sqlite/chembl_34.db --NPs_data=coconut_complete-10-2024.csv --ChEMBL --NPs

    4. To update only ChEMBL:
       python update_versions_ChEMBL_NPs.py --ChEMBL_data=chembl_34/chembl_34_sqlite/chembl_34.db --ChEMBL

Prepared By: Alhasbary
Date: 10/10/2024
"""

import os
import pandas as pd
import time
import argparse
import sys
from SharedFunc import *

def parse_args():   
    parser = argparse.ArgumentParser(description="Preprocess and update refined versions of ChEMBL and NPs data.")
    
    parser.add_argument("--ChEMBL_data", action="store", dest='ChEMBL_dataPath', type=str, 
                        help='Full filepath of the SQLite file containing the ChEMBL database to be processed and refined.')
    
    parser.add_argument("--NPs_data", action="store", dest='NPs_dataPath', type=str,
                        help='Full filepath of the CSV file containing the latest version of natural products (NPs) data to be processed.')
                        
    parser.add_argument("--ChEMBL", action="store_true", default=False,
                        help="Flag to update the ChEMBL database. Default is False.")
                        
    parser.add_argument("--NPs", action="store_true", default=False,
                        help="Flag to update NPs data. Default is False.")
    
    parser.add_argument("--destination", action="store", dest='destinationPath', type=str, 
                        default='Data', help='Full filepath for saving the refined Parquet files. Default is "Data".')

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # Create the destination folder if it does not exist
    os.makedirs(args.destinationPath, exist_ok=True) 

    # Begin logging: create or append to the log file
    logPath = os.path.join(args.destinationPath, "update_chembl_np_log.txt")
    _logFile = open(logPath, 'a')  
    log = Tee(sys.stdout, _logFile) 

    # Log the script execution details
    print("\n", file=log)
    print(" ".join(sys.argv), file=log)  
    print(f"Start time: {tstamp()}", file=log)  

    # Ensure the output directory for refined ChEMBL data exists
    output_dir = os.path.join(args.destinationPath, "refined_chembl")
    os.makedirs(output_dir, exist_ok=True)  

    # Define file path for the CTA_meta Parquet file and check if it exists
    parquet_meta_path = os.path.join(output_dir, "CTA_meta.parquet")
    if not os.path.exists(parquet_meta_path):
        # If CTA_meta does not exist, create both CTA_meta and CTA_data files
        initialize_CTA_files(log, output_dir)
        
    # Process and refine the ChEMBL database, saving the output as Parquet files
    if args.ChEMBL:
        if args.ChEMBL_dataPath:
            #refine_chembl_database(log, args.ChEMBL_dataPath, output_dir)
            set_ChEMBL_as_first_option_in_CTA_table(log, output_dir)
        else:
            print("Please provide the full file path of the ChEMBL SQLite file using the --ChEMBL_data option.")

    # Process and refine the NPs data, saving the output as Parquet files
    if args.NPs:
        if args.NPs_dataPath:
            refine_NPs_data(log, args.NPs_dataPath, args.destinationPath, args.verbose)
        else:
            print("Please provide the CSV file with the latest NPs data using the --NPs_data option.")

    # Log the script end time
    print(f"End time: {tstamp()}", file=log)
    
