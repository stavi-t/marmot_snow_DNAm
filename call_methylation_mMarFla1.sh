#!/bin/bash
#SBATCH --job-name=callmeth
#SBATCH --output=../slurm/call_chr_cov_%j.out
#SBATCH --error=../slurm/call_chr_cov_%j.err
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8    # 8 threads specified
#SBATCH --mem-per-cpu=4GB 	 # per task by default with array (per cpu if multithreaded)
#SBATCH --time=48:00:00      # per task by default
#SBATCH --array=1-119

module load anaconda3/2024.2
conda activate bioinformatics

# Methylation calls for WGBS and targeted bisulfite sequencing can be improved by removing PCR duplicate reads before calling methylation. 
# Marking duplicate alignments w/ RRBS data is not recommended as the sequencing reads will often share the same mapping coordinates due to enzymatic digestion
# Must specify a max and min read depth for calling: DEFAULT -min=10, -max=8000

# set paths
SAMPLE_LIST="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/filenames_2.txt"
aligned_dir="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/aligned_mMarFla1"
called_dir="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/called_mMarFla1/filter_10x-81x"
indexed_genome="/home/vonholdt/VONHOLDT/stennen/RMBL-MARMOTS/bsbolt/indexed_chr_genome"

# Read the sample name corresponding to the array task ID + trim white space
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST | tr -d '\n')

# Set the input BAM and sorted BAM file paths
BAM="$aligned_dir/${SAMPLE}.bam"
SORTED_BAM="$aligned_dir/${SAMPLE}.sorted.bam"

# # sort bams
# sambamba sort -t 12 -o $SORTED_BAM $BAM

# # index sorted bams
# sambamba index $SORTED_BAM

# Methylation Calling with 8 threads in verbose mode
python3 -m bsbolt CallMethylation -I $SORTED_BAM -O $called_dir/${SAMPLE} -DB $indexed_genome -t 8 -min 10 -max 81 -verbose > $called_dir/${SAMPLE}_stats.txt




