#!/bin/bash
#SBATCH --job-name=marmotindex
#SBATCH --error=index_mMarFla1.err
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G 	#per task
#SBATCH --time=6:00:00

module load anaconda3/2024.2
conda activate bioinformatics

# RRBS Index, MspI Cut Format, default 40bp Lower Fragment Bound, and 500bp Upper Fragment Bound - adjusting this to reflect previous ref genomes

python3 -m bsbolt Index -G mMarFla1.hap1.fna -DB /home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/indexed_chr_genome \
-rrbs -rrbs-cut-format C-CGG -rrbs-lower 50 -rrbs-upper 500