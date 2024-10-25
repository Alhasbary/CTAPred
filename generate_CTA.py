"""
This script is part of the Compound-Target Activity Prediction (CTAPred) Program, designed to generate a CTA dataset for use as a reference in similarity-based target identification.

Inputs:
   - List of optional parameters.
   - --output: Full path to save the results [Optional].

Outputs:
   - Updated "CTA_data.parquet" and "CTA_meta.parquet" files located in the "Data/mini_chembl/" directory.
   - A CSV file saved in the "output/CTA_datasets/" directory.

Usage:
   - For help:
     python generate_CTA.py -h
   - To execute the script:
     python generate_CTA.py [parameters] --output=FullPathToSaveResults

Example:
   python generate_CTA.py --fingerprint 'ecfp' --nBits 1024 --radius 2 --sv 1000 --Tc 0.30

Prepared By: Alhasbary
Date: 10/10/2024
"""

import sys
from SharedFunc import *


def parse_args():           
    parser = argparse.ArgumentParser()
    parser.add_argument("--sv", type=float_range(0.01, 10000.0), default=10000, action="store", dest="standard_value",
                        help="Activity values in nM measurement. default=10000 nM.")
                        
    parser.add_argument("--Tc", action="store", type=float_range(0.1, 1.0), default=0.80, dest="Tc",
                        help="Desired value for CTA 'Tc' similarity threshold (0.1-1.0). default=0.85.")
                        
    parser.add_argument("--fingerprint", action="store", default='maccs', dest="fingerprint",
                        choices=['avalon', 'ecfp', 'fcfp', 'maccs'],
                        help="Desired fingerprint type (avalon, ecfp, fcfp, or maccs). default='ecfp'")
                        
    parser.add_argument("--nBits", action="store", type=int_range(8, 2048), default=2048, dest="nBits",
                        help="Number of bits parameter that specifies the length of the generated (avalon, ecfp, "
                             "or fcfp)fingerprint. default=2048")
                             
    parser.add_argument("--radius", action="store", type=int, default=2, choices=[2, 3], dest="radius",
                        help="Desired radius value for Morgan ECFP/FCFP fingerprints (2, 3). default=2")
                        
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Enable this flag to print error messages during SMILES preprocessing. Default is False.")
                        
    parser.add_argument("--n_jobs", action="store", type=int, default=-1, dest="n_jobs",
                        help="Number of CPU cores to use. default=-1 to use all available CPU cores.")
                        
    parser.add_argument("--output", action="store", dest='outputPath', type=str,
                        default='output', help='Full filepath, in which the CTA data is saved. Default is output')
    args = parser.parse_args()
    return args
    
    
def main():
    # get inputs
    args = parse_args()

    # create a folder to save the results if it doesn't exist
    if not os.path.exists(args.outputPath):
        os.makedirs(args.outputPath)
    # begin log
    logPath = os.path.join(args.outputPath, "generate_CTA_log.txt")
    _logFile = open(logPath, 'a')
    log = Tee(sys.stdout, _logFile)

    print("\n", file=log)
    print(" ".join(sys.argv), file=log)
    print("Start time: %s " % (tstamp()), file=log)
             
           
    # Ensure directory exists
    output_dir = os.path.join(args.outputPath, "CTA_datasets/")
    os.makedirs(output_dir, exist_ok=True)

    # Determine the value of the 'radius' and nBits field based on the fingerprint type
    if args.fingerprint in ['avalon', 'maccs']:
        args.radius = None
        if args.fingerprint == 'maccs':
            args.nBits = 116    
                        
    createCTA(log, args, output_dir)    
    print("\nFinish time: %s " % (tstamp()), file=log)
    _logFile.close()

if __name__ == '__main__':
    main()

 