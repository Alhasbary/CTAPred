# Compound-Target Activity Prediction (CTAPred) Tool

CTAPred is an open-source command-line tool for predicting potential protein targets for natural products. It utilizes a two-stage process that integrates molecular fingerprinting and similarity-based searches to identify candidate drug targets.

## Key Features

### 1. **extractDB: Bioactivity Data Extraction**
Retrieves high-quality bioactivity data from ChEMBL with custom filters and activity thresholds, ensuring only relevant records are included.

### 2. **processMS: Molecular Structure Standardization**

Standardizes molecular structures using MolVS, ensuring consistent SMILES representation by removing hydrogens, normalizing metals, reionizing, and canonicalizing tautomers while preserving stereochemistry.

### 4. **createCTA: Compound-Target Activity Dataset Construction**
Generates a CTA reference dataset by performing chemical similarity searches between ChEMBL and natural product compounds. SMILES are converted into a user-selected fingerprint (Avalon, Morgan ECFP, Morgan FCFP, or MACCS). Compounds showing a Tanimoto coefficient above a specified threshold linked to their protein targets in the CTA dataset.

### 5. **rankTargets: Ranking Potential Targets**

Ranks protein targets in the CTA dataset based on the average similarity scores between the query compound and the top user-defined number of similar compounds per target.

## Parameters
CTAPred offers customizable parameters for detailed analysis:

| Parameter        | Description                                                                                                                                          
|------------------|---------------------------------------------------------------------------------------------------------------------|
| `sv`             | Activity threshold in nM. Default is `10000` nM                                                                     |
| `fingerprint`    | Fingerprint type (avalon, ecfp, fcfp, or maccs). Defaults: `maccs` for CTA creation, `ecfp` for target prediction.  |
| `nBits`          | Number of fingerprint bits (avalon, ecfp, or fcfp). Defaults: `166` for CTA creation, `2048` for target prediction. |
| `radius`         | Radius for ECFP/FCFP fingerprints (2 or 3). Default is `2`.                                                         |
| `Tc`             | Similarity threshold for Tanimoto coefficient (0.1–1.0). Default is `0.80`                                          |
| `k`              | Top-k reference compounds for ranking. Default is`1`                                                                |
| `n_jobs`         | CPU cores for processing. Default is `-1` (uses all available cores).                                               |
| `verbose`        | Flag for detailed error messages during SMILES preprocessing. Default is `False`.                                   |
| `inputPath`      | Path for input data files with NP query lists. Default is `input`                                                   |
| `outputPath`     | Path to save results and CTA data. Default is `output`                                                              |

## System Requirements

- Python: >= 3.7
- RDKit: >= 2020.09.1 
- MolVS: >= 0.1.1 
- DuckDB: >= 1.1.1 
- Fastparquet: >= 0.5.0 
- NumPy: >= 1.21.5 
- Pandas: >= 1.3.5 
- Joblib: >= 1.1.0 
- tqdm: >= 4.66.4 
- SQLite3: >= 2.6.0 

## Installation

1. Download the repository and extract files.
2. Open a command-line interface.
3. Run the Python scripts as outlined below.

## Usage

          
### Identifying Potential Targets


Inputs:
   - Optional parameters.
   - `--input`: Full path to the directory containing SMILES string files (e.g., QueryList1_smiles.csv) [Optional].
              The directory can include multiple lists, which the script will process sequentially. File names should follow the pattern: QueryListN_smiles.csv, where N is an integer.
              Each CSV file must contain a 'smiles' column for SMILES strings and a 'np_id' column for compound IDs.
   - `--output`: Full path to save the results [Optional].

Output:
   - CSV files with predicted targets.

Usage:
   - For help:
     ```
     python predict_targets.py -h
     ```
   - To run the program:
     ```
     python predict_targets.py [parameters] --input=FullPathToInputFiles --output=FullPathToSaveResults
     ```

Example:
     ```
   python predict_targets.py --fingerprint 'ecfp' --nBits 1024 --k 1 3 
     ```
   
### Creating a Compound-Target Activity (CTA) Dataset


Inputs:
   - Optional parameters.
   - `--output`: Full path to save results [Optional].

Outputs:
   - Updated `CTA_data.parquet` and `CTA_meta.parquet` in `Data/mini_chembl/`.
   - A CSV file in `output/CTA_datasets/`.

Usage:
   - For help:
     ```
     python generate_CTA.py -h
     ```
   - To execute the script:
     ```
     python generate_CTA.py [parameters] --output=FullPathToSaveResults
     ```

Example:
     ```
   python generate_CTA.py --fingerprint 'ecfp' --nBits 1024 --radius 2 --sv 1000 --Tc 0.30
     ```
 
### Updating ChEMBL and NP Datasets

Inputs:
   - --ChEMBL_data : Path to ChEMBL SQLite database [Required].
   - --NPs_data    : Path to NP data CSV [Required].
   - --destination : Directory path for saving refined Parquet files [Optional, default='Data'].
   - --ChEMBL      : Flag for updating ChEMBL.
   - --NPs         : Flag for updating Natural Products data.

Outputs: 
   - 'refined_chembl' folder containing Parquet files.
   - 'NPs_clean_smiles.parquet' with refined NP data.

Usage:
   - For help:
     python update_versions_ChEMBL_NPs.py -h
   - To run the program:
     python update_versions_ChEMBL_NPs.py --ChEMBL_data=FullPathToChEMBL --NPs_data=FullPathToNPs --destination=FullPathToSaveParquet

Usage Example:

1. Download ChEMBL SQLite dataset from [ChEMBL Downloads](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/).
   Extract using:
     ``` 
     tar -xvzf chembl_34_sqlite.tar.gz 
     ```

2. Download COCONUT dataset from [COCONUT Downloads](https://coconut.naturalproducts.net/download/).
  Extract using: 
     ```
     tar -xvzf coconut-10-2024.csv.zip
     ```

3. To update both datasets:
     ```
     python update_versions_ChEMBL_NPs.py --ChEMBL_data=chembl_34/chembl_34_sqlite/chembl_34.db --NPs_data=coconut_complete-10-2024.csv --ChEMBL --NPs
     ```

4. To update only ChEMBL:
     ```
     python update_versions_ChEMBL_NPs.py --ChEMBL_data=chembl_34/chembl_34_sqlite/chembl_34.db --ChEMBL
     ```
       
    
Feel free to explore CTAPred with your query compound lists. Please ensure that each query list is provided as a CSV file, with smiles strings in a column named 'smiles' and compound IDs in a column named 'np_id'.
