#!/bin/bash
#SBATCH --job-name=bammerge33
#SBATCH --output=slurm-out/merge_%j.out
#SBATCH --error=slurm-out/merge_%j.err
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB 	         # per task by default with array
#SBATCH --time=24:00:00      # per task by default

module load anaconda3/2024.2
conda activate bioinformatics

aligned_dir="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/aligned"

# sort bams by coordinate to merge 
sambamba sort -t 4 $aligned_dir/CO_MARM1-17_S33_part_1_aligned.bam
sambamba sort -t 4 $aligned_dir/CO_MARM1-17_S33_part_2_aligned.bam
sambamba sort -t 4 $aligned_dir/CO_MARM1-17_S33_part_3_aligned.bam
sambamba sort -t 4 $aligned_dir/CO_MARM1-17_S33_part_4_aligned.bam
sambamba sort -t 4 $aligned_dir/CO_MARM1-17_S33_part_5_aligned.bam
sambamba sort -t 4 $aligned_dir/CO_MARM1-17_S33_part_6_aligned.bam
sambamba sort -t 4 $aligned_dir/CO_MARM1-17_S33_part_7_aligned.bam
sambamba sort -t 4 $aligned_dir/CO_MARM1-17_S33_part_8_aligned.bam
sambamba sort -t 4 $aligned_dir/CO_MARM1-17_S33_part_9_aligned.bam
sambamba sort -t 4 $aligned_dir/CO_MARM1-17_S33_part_10_aligned.bam

# merge 10 split bams to one to call methylation
sambamba merge -t 4 --show-progress $aligned_dir/CO_MARM1-17_S33_aligned_merge.bam $aligned_dir/CO_MARM1-17_S33_part_1_aligned.sorted.bam \
$aligned_dir/CO_MARM1-17_S33_part_2_aligned.sorted.bam \
$aligned_dir/CO_MARM1-17_S33_part_3_aligned.sorted.bam \
$aligned_dir/CO_MARM1-17_S33_part_4_aligned.sorted.bam \
$aligned_dir/CO_MARM1-17_S33_part_5_aligned.sorted.bam \
$aligned_dir/CO_MARM1-17_S33_part_6_aligned.sorted.bam \
$aligned_dir/CO_MARM1-17_S33_part_7_aligned.sorted.bam \
$aligned_dir/CO_MARM1-17_S33_part_8_aligned.sorted.bam \
$aligned_dir/CO_MARM1-17_S33_part_9_aligned.sorted.bam \
$aligned_dir/CO_MARM1-17_S33_part_10_aligned.sorted.bam