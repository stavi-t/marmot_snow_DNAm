#!/bin/bash
#SBATCH --job-name=spltalign2
#SBATCH --output=spltalign2_%A_%a.out
#SBATCH --error=spltalign2_%A_%a.err
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB 	#per task by default with array
#SBATCH --time=72:00:00     #per task by default
#SBATCH --array=1-10

module load anaconda3/2024.2
conda activate bioinformatics

trimmed_dir="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/01_trimmed"

# first run seqkit to split large .fastq files into 10 parts per read = 20 files total / sample
# seqkit split2 -p 10 -O $trimmed_dir/split_CO_MARM1-18_S34 -1 $trimmed_dir/CO_MARM1-18_S34_R1_trimmed.fastq -2 $trimmed_dir/CO_MARM1-18_S34_R2_trimmed.fastq
# seqkit split2 -p 10 -O $trimmed_dir/split_CO_MARM1-17_S33 -1 $trimmed_dir/CO_MARM1-17_S33_R1_trimmed.fastq -2 $trimmed_dir/CO_MARM1-17_S33_R2_trimmed.fastq

conda deactivate
conda activate pip-bsbolt

# Define the sample name
SAMPLE="CO_MARM1-18_S34"

# then set paths to paired, split files and outputs
READ1_SPLIT=$(ls $trimmed_dir/split_${SAMPLE}/${SAMPLE}_R1*.fastq | sed -n "${SLURM_ARRAY_TASK_ID}p")
READ2_SPLIT=$(ls $trimmed_dir/split_${SAMPLE}/${SAMPLE}_R2*.fastq | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo "Sample: $SAMPLE"
echo "Read 1 Parts: $READ1_SPLIT"
echo "Read 2 Parts: $READ2_SPLIT"

OUTPUT_BAM="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/aligned/${SAMPLE}_part_${SLURM_ARRAY_TASK_ID}_aligned"

# then align the 10 sections of each file in parallel with 4 threads 
# Paired End Alignment Using Default Commands with indexed genome
python3 -m bsbolt Align -DB /home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/indexed_genome -F1 $READ1_SPLIT -F2 $READ2_SPLIT -t 4 -O $OUTPUT_BAM