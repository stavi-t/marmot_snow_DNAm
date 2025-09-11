#!/bin/bash
#SBATCH --job-name=align119
#SBATCH --output=../slurm/align_chr_%j.out
#SBATCH --error=../slurm/align_chr_%j.err
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4GB 	        #per task by default with array (per cpu if multithreaded)
#SBATCH --time=24:00:00             #per task by default
#SBATCH --array=1-119

module load anaconda3/2024.2
conda activate bioinformatics

# align RR EMseq trimmed reads to new CHROMOSOME LEVEL REF GENOME: mMarFla1.hap1
SAMPLE_LIST="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/filenames_2.txt"

# Read the sample name corresponding to the array task ID + trim white space
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST | tr -d '\n')

# trimmed paired end reads from 01_trimmed original files for L1,L2 (*_trimmed == 70)
# trimmed paired end reads from 01_trimmed/concat_new files for L3,L4,L5 (*_trim == 119)
FASTQ1="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/01_trimmed/${SAMPLE}_R1_trim.fastq"
FASTQ2="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/01_trimmed/${SAMPLE}_R2_trim.fastq"
echo "Sample: $SAMPLE"

# Paired End Alignment Using Default Commands with new indexed chromosome level ref genome
python3 -m bsbolt Align -DB /home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/indexed_chr_genome -F1 $FASTQ1 -F2 $FASTQ2 -t 8 -O /home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/aligned_mMarFla1/${SAMPLE}