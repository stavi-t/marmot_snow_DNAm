# Code used in Tennenbaum et al. 202X: 

"Climate change increases energetic demand, shaping epigenetic regulation of metabolism in a hibernating mammal"

## Metadata cleaning
0. Assemble sample metadata

## DNA methylation data
1. Shell scripts to analyze RR-EMseq data (BSBOLT), (NEB kit v1 + MspI digest): make_index.sh, trim_read.sh, align_reads.sh, call_methylation.sh, aggregate_matrix.sh
2. Shell/python scripts to analyze control DNA (BSSEEKER2): lambda.sh, puc19.sh, lambda_summarize.py, puc19_summarize.py
3. Python scripts to analyze methylation variance and coverage: coverage_PCAs.py, coverage_violins.py
4. R script to run binomial GLMMs in lme4-breeding: 575k_snow.R
5. R script to analyze model output, variance, identify outlier CpGs: LME4Breed_out_v2.R
6. R script to visualize CpG methylation in candidate promoters: prom_characterization.R
7. R script to analyze GO/KEGG functional enrichment from eggNOG annotations: eggnog_2.R
