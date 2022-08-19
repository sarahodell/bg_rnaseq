# `scripts` Directory Contents


## RNASeq QC
- `qc_check.R`
- `deseq.R`
- `limma_voom.R`: normalization of gene counts from STAR
- `expression_matrix.R`
- `expression_plots.R`
- `ft_correlation.R`
- `read_check.R`
- `mapping_rates.py`
- `alignment_qc.R`

## *cis*-eQTL
- `ciseqtl_GridLMM.R`: running GridLMM on individual genes to identify cis-eqtl
- `cis_eqtl_plots.R`: creates a Manhattan plot of all the eQTL tests from `ciseqtl_GridLMM.R`
- `beta_fix.R`: make one txt file from results of `ciseqtl_GridLMM.R` for each chromosome, as well as the residuals after accounting for cis-variation. Used as input for running MegaLMM.
- `QTL_ciseQTL_overlap.R`: tests for positional overlap between previously identified agronomic QTL and cis-eQTL. Also looks to see if the effect sizes of the parental alleles are correlated for overlapping QTL.
- `var_explained.R`: calculates the amount of variation in gene expression is explained by cis-eQTL
- `check_h2.R`: looking at the heritability of cis-eQTL

## Running MegaLMM_
- `runMegaLMM.R`: running with gene counts for a single timepoint
- `runMegaLMM_test.R`: testing smaller subset of genes for troubleshooting
- `run_pheno_MegaLMM.R`: running with agronomic phenotypes included with gene counts
- `avg_posterior.R`: getting posterior means of parameter estimates from runMegaLMM.R, pairing correlated factors across chains
- `factor_groups.R`: identify genes and agronomic phenotypes that are grouped into the same factors by runMegaLMM.R, and how much variation is explained by each

## *trans*-eQTL
- `transeQTL_GridLMM.R`: running GridLMM on factors output from `run_MegaLMM.R` (F parameter) to identify trans-eqtl
- `trans_eqtl_plots.R`: creates a Manhattan plot of all the eQTL tests from `transeqtl_GridLMM.R`
