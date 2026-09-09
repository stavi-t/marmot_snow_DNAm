#!/bin/bash
#SBATCH --job-name=trim
#SBATCH --output=slurm/trim_%j.out
#SBATCH --error=slurm/trim_%j.err
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1GB 	        #per task by default with array (per cpu if multithreaded)
#SBATCH --time=4:00:00              #per task by default
#SBATCH --array=1-119

# samples names
SAMPLE_LIST="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/filenames_2.txt"

# Read the sample name corresponding to the array task ID + trim white space
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST | tr -d '\n')

# activate env with cutadapt
module load anaconda3/2024.2
conda activate bioinformatics

FASTQ1=(/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/raw/concat_new/${SAMPLE}_R1_001.fastq.gz)
FASTQ2=(/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/raw/concat_new/${SAMPLE}_R2_001.fastq.gz)
out_dir="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/01_trimmed"

echo "$SAMPLE"
echo "Read 1: $FASTQ1"
echo "Read 2: $FASTQ2"

# The NEBNext libraries for Illumina resemble TruSeq libraries and can be trimmed similar to TruSeq:
# Adaptor Read 1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA Adaptor Read 2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

cutadapt -q 20 --minimum-length 20 -a NNAGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A NNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT $FASTQ1 $FASTQ2 --cores=8 -o $out_dir/${SAMPLE}_R1_trim.fastq -p $out_dir/${SAMPLE}_R2_trim.fastq
