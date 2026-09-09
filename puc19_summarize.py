#!/usr/bin/env python3
#shift-enter to run interactively
import os
import gzip
import glob

# define directories and summary file path
summary_file="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/04_controls/scripts/puc19_conversion_CO_MARM1.txt"
cgmap_dir = "/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/04_controls"

# Use all puc19 cgmap files
files = glob.glob(os.path.join(cgmap_dir, "CO_MARM1*_puc19_10x.CGmap.gz"))

# read the CGmap file for each sample with wildcard
for cgmap_file in files:
    sample = os.path.basename(cgmap_file).split("_puc19_10x.CGmap.gz")[0]  # Extract sample name from CGmap prefix
    print(f"Processing file: {sample}")
    
# initialize variables to calculate averages
    sum_cg = 0
    count_cg = 0
    sum_other = 0
    count_other = 0

    # Process CGmap files
    with gzip.open(cgmap_file, "rt") as f:
        for line in f:
            columns = line.strip().split()
            context = columns[3]            # dinucleotide context
            methylation = float(columns[5]) # methylation proportion in CGmap
            
            if context == "CG":    # Context is only CGs
                sum_cg += methylation
                count_cg += 1
            else:                  # Non-CG context
                sum_other += methylation
                count_other += 1
    
    # # Skip files with no valid rows for both contexts
    # if count_cg == 0 and count_other == 0:
    #     print(f"No valid rows in file: {cgmap_file}")
    #     continue

    # calculate both averages for valid samples
    avg_cg = sum_cg / count_cg
    avg_other = sum_other / count_other

    # calculate the puc19 CpG conversion rate
    conversion_rate_cg = 1 - avg_cg

    # append results with write "a" to the summary file
    with open(summary_file, "a") as summary_out:
        summary_out.write(f"File: {sample}, MC Avg CGs: {avg_cg:.2f}, MC Avg Other: {avg_other:.2f}, pUC19 Conversion Rate: {conversion_rate_cg:.2f}\n")
    
print("Finished writing summary file.")
