# Code used in Tennenbaum et al. 202X: 

"Climate change increases energetic demand, shaping epigenetic regulation of metabolism in a hibernating mammal"

## Metadata cleaning
0. Assemble sample metadata

## DNA methylation data
1. Shell scripts to analyze RR-EMseq data with MspI digest (NEB kit v1): trim, align, call methylation, aggregate matrix
2. Shell scripts to analyze RR-Emseq control DNA: lambda and pUC19 align, call methylation, aggregate matrix, summarize results
3. Python scripts to plot PCAs and analyze sequencing coverage
4. R script to run binomial GLMMs in lme4-breeding
5. R script to analyze model output, variance, and identify outlier CpGs
6. R script to visualize CpG methylation across candidate promoters
7. R script to analyze GO/KEGG functional enrichment from eggNOG annotations
