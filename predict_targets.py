"""
This script is part of the Compound-Target Activity Prediction (CTAPred) Program, designed to identify potential target(s) for chemical compounds.

Inputs:
   - Optional parameters list.
   - --input: Full path to the directory containing SMILES string files (e.g., QueryList1_smiles.csv) [Optional].
              The directory can include multiple lists, which the script will process sequentially. File names should follow the pattern: QueryListN_smiles.csv, where N is an integer.
              Each CSV file must contain a 'smiles' column for SMILES strings and a 'np_id' column for compound IDs.
   - --output: Full path to save the results [Optional].

Output:
   - A CSV file for each query list (QueryListN_potential_targets_based_on_fp_with_k_value_of_*.csv) identifying potential targets based on the mean similarity scores of the top k similar reference compounds.

Usage:
   - For help:
     python predict_targets.py -h
   - To run the program:
     python predict_targets.py [parameters] --input=FullPathToInputFiles --output=FullPathToSaveResults

Example:
   python predict_targets.py --fingerprint 'ecfp' --nBits 1024 --k 1 3 

Prepared By: Alhasbary
Date: 10/10/2024
"""

import argparse
import sys
from SharedFunc import *


def parse_args():           
    parser = argparse.ArgumentParser(description="Predict potential protein target(s).")
    parser.add_argument("--fingerprint", action="store", default='ecfp', dest="fingerprint",
                        choices=['avalon', 'ecfp', 'fcfp', 'maccs'],
                        help="Desired fingerprint type (avalon, ecfp, fcfp, or maccs). Default='ecfp'")
                        
    parser.add_argument("--nBits", action="store", type=int_range(8, 2048), default=2048, dest="nBits",
                        help="Number of bits parameter that specifies the length of the generated (avalon, ecfp, "
                             "or fcfp)fingerprint. Default=2048")
                             
    parser.add_argument("--radius", action="store", type=int, default=2, choices=[2, 3], dest="radius",
                        help="Desired radius value for Morgan ECFP/FCFP fingerprints (2, 3). Default=2")
                        
    parser.add_argument("--k", action="store", type=int, default=1, dest="top_k", nargs='*',
                        help="Desired value for 'top-k' reference compounds")
                        
    parser.add_argument("--n_jobs", action="store", type=int, default=-1, dest="n_jobs",
                        help="Number of CPU cores to use. Default=-1 to use all available CPU cores.")
                        
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Enable this flag to print error messages during SMILES preprocessing. Default is False.")
                    
    parser.add_argument("--input", action="store", dest='inputPath', type=str, default='input',
                        help='Full filepath of the data files, which contains the NP query lists. Default is input.')
                        
    parser.add_argument("--output", action="store", dest='outputPath', type=str,
                        default='output', help='Full filepath, in which the potential targets is saved. Default is output')
    args = parser.parse_args()
    return args
    
    
def main():
    # get inputs
    args = parse_args()

    # create a folder to save the results if it doesn't exist
    if not os.path.exists(args.outputPath):
        os.makedirs(args.outputPath)
    # begin log
    logPath = os.path.join(args.outputPath, "predict_targets_log.txt")
    _logFile = open(logPath, 'a')
    log = Tee(sys.stdout, _logFile)

    print("\n", file=log)
    print(" ".join(sys.argv), file=log)
    print("Start time: %s " % (tstamp()), file=log)
        
    if not assert_query_resource_files(log, args):
        print("\nFinish time: %s" % (tstamp()), file=log)
        _logFile.close()
        exit()    
    
    # Determine the value of the 'radius' and nBits field based on the fingerprint type
    if args.fingerprint in ['avalon', 'maccs']:
        args.radius = None
        if args.fingerprint == 'maccs':
            args.nBits = 116    
            
    # Retrieved the CTA reference dataset 
    CTA , CTA_id = select_CTA_reference_dataset(log, args)
    #CTA = CTA[:1000]
    CTA = CTA[['molregno', 'clean_smiles', 'uniprot', 'target_chembl_id']].drop_duplicates(ignore_index=True)
    CTA.dropna(inplace=True)                   
    CTA_fps = gen_fps(log, args, CTA[['molregno', 'clean_smiles']], 'Reference').drop_duplicates(ignore_index=True)
    CTA = pd.merge(CTA, CTA_fps, on='molregno', how='inner')
    CTA.dropna(inplace=True)     

    print(f"\nStart processing the {args.datNames} Query sets:", file=log)       
    
    for db_name in args.datNames:
        print(f"\n==================\nData set: {db_name}\n==================\n", file=log)
        smiles = os.path.join(args.inputPath, db_name+"_smiles.csv")
        if not os.path.exists(smiles):
            print("Error, file not found", smiles)
            print("It must be a csv file which contains smiles strings in a column named 'smiles' and compound id in "
                  "a column named 'np_id'")
            break
    
        # Read the smiles file into a pandas DataFrame
        query_list = pd.read_csv(smiles)
        clean_query_structures = preprocess_Query_SMILES(log, args, query_list)
        
          
        clean_query_structures.rename(columns={'np_id': 'molregno'}, inplace=True)
                      
        query_fps = gen_fps(log, args, clean_query_structures, 'Query').drop_duplicates(ignore_index=True)
        
        # Chemical similarity searches between two data sets.
        query_fps.rename(columns={'molregno': 'np_id'}, inplace=True)
        ChEMBL_similarity = Query_conductSS(log, args, query_fps[['np_id', 'fp']], CTA[['molregno', 'target_chembl_id', 'uniprot', 'fp']])    
        
                    
        if len(ChEMBL_similarity) > 0:            
           
            # Ensure directory exists
            output_dir = os.path.join(args.outputPath, "potential_targets/")
            os.makedirs(output_dir, exist_ok=True)
            # Retrieve potential targets for each smiles string according to the mean similarity scores of K nearest neighbors.
            potential_targets_path = os.path.join(output_dir, f"{db_name}_using_dataset_id_{CTA_id}_include_with_k_value_of")
            rankTargets(log, args, ChEMBL_similarity, potential_targets_path)   
            
            del ChEMBL_similarity
            gc.collect()
                     
    print("\nFinish time: %s " % (tstamp()), file=log)
    _logFile.close()


if __name__ == '__main__':
    main()

 