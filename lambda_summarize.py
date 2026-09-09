#!/usr/bin/env python3
#shift-enter to run interactively
import os
import gzip
import glob

# define cgmap directory and summary file path
summary_file="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/04_controls/scripts/lambda_conversion_CO_MARM1.txt"
cgmap_dir = "/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/04_controls"

# Use all lambda cgmap files
files = glob.glob(os.path.join(cgmap_dir, "*CO_MARM1*_lambda_10x.CGmap.gz"))

# read the CGmap file for each sample with wildcard
for cgmap_file in files:
    sample = os.path.basename(cgmap_file).split("_lambda_10x.CGmap.gz")[0]  # Extract sample name from CGmap prefix
    print(f"Processing file: {sample}")
  
  # initialize variables to calculate averages
    sum_lambda = 0
    sum_rows = 0

  # process CGmap files
    with gzip.open(cgmap_file, "rt") as f: # open zipped files as text
        for line in f:
            columns = line.strip().split()
            methylation = float(columns[5]) # methylation proportion CGmap
            
            sum_lambda += methylation
            sum_rows += 1
    
    avg_lambda = sum_lambda / sum_rows

    # calculate the conversion rate
    conversion_rate = 1 - avg_lambda
    
    # append results with write "a" to the summary file
    with open(summary_file, "a") as summary_out:
        summary_out.write(f"File: {sample},  MC Average: {avg_lambda:.2f}, Lambda Conversion Rate: {conversion_rate:.2f}\n")

print("Finished writing summary file.")
