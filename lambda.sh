#!/bin/bash
#SBATCH --job-name=lambdaconv
#SBATCH --output=../../slurm/lambda_%j.out
#SBATCH --error=../../slurm/lambda_%j.err
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G 			#per task by default in array
#SBATCH --time=8:00:00      #per task by default in array
#SBATCH --array=1-119

#must build lambda and pUC19 ref genomes in default dirr, specify location with -d

SAMPLE_LIST="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/filenames_2.txt"

# Read the sample name corresponding to the array task ID + trim white space
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST | tr -d '\n')

module load anaconda3/2024.2
conda activate bioinfo_py2

FASTQ1="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/01_trimmed/${SAMPLE}_R1_trim.fastq"
FASTQ2="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/01_trimmed/${SAMPLE}_R2_trim.fastq"
echo $SAMPLE

# set paths
lambda_genome="/home/vonholdt/VONHOLDT/stennen/ncbi_ref_genomes/lambda_NC_001416.1.fa"
index_path="/home/vonholdt/VONHOLDT/stennen/BSseeker2/bs_utils/reference_genomes"
wgbs_index="/home/vonholdt/VONHOLDT/stennen/BSseeker2/bs_utils/reference_genomes/lambda_NC_001416.1.fa_bowtie2"
out_dir="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/04_controls"

# align reads to lambda genome
/home/vonholdt/VONHOLDT/stennen/BSseeker2/bs_seeker2-align.py -1 $FASTQ1 -2 $FASTQ2 -g $lambda_genome --aligner=bowtie2 -d $index_path -o $out_dir/${SAMPLE}_lambda_aligned.bam

# call methylation at 10x in lambda genome
/home/vonholdt/VONHOLDT/stennen/BSseeker2/bs_seeker2-call_methylation.py -i $out_dir/${SAMPLE}_lambda_aligned.bam --CGmap=$out_dir/${SAMPLE}_lambda_10x.CGmap.gz -d $wgbs_index -r 10 
