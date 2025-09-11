#!/bin/bash
#SBATCH --job-name=aggmatrix
#SBATCH --output=../slurm/matrix_chr_0.9_%j.out
#SBATCH --error=../slurm/matrix_chr_0.9_%j.err
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4   #4 threads specified
#SBATCH --mem=15GB 	        #per task by default
#SBATCH --time=6:00:00      #per task by default

module load anaconda3/2024.2
conda activate bioinformatics

# takes a list of CGmap files and assembles a consensus methylation matrix
# methylated sites that pass a read depth threshold (10 is default=count) and are present in a set proportion of samples are included in the matrix (0.8 is default=float)

# CGmap files are first iterated through to count sites and THEN iterated to through to retrieve values - ***this is memory intensive***!

output_dir="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/clock"

# if -S sample names are not provided as comma separated string OR path to txt with line separated sample labels, sample names are extracted from CGmap files. 
# List of CGmap files MUST CONTAIN FULL PATH!
# -count = outputs a count matrix with count of methylated cytosines and total observed cytosines
# -CG = only output CG sites

python3 -m bsbolt AggregateMatrix -F CGmaps_omit_mMarFla1.txt -S S_names.txt -min-coverage 10 -CG -t 4 -min-sample 0.9 -O $output_dir/marmot_chr_mMarFla1_0.9_165_freq.txt
