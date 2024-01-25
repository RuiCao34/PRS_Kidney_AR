Unfortunately, due to the consent limitation, the phenotype and genotype data in these two cohorts are not publicly available. To reproduce the PRS result, run the R file sequetially:
1. GWAS_dekaf_HLAmatching.R： parallel GWAS analysis on DeKAF cohort. Cox regression model: AR ~ recipient_age + recipient_gender + PRA_status + previous_tx_status + HLA_num_mismatch.
2. writeplink_dekaf_HLAmatching.R: clump the GWAS results.
3. PRS_dekaf_HLAmatching.R, PRS_gen03_HLAmatching.R: calculate whole-genome PRS scores in both DeKAF and GEN-03 cohorts.
4. PRS_dekaf_HLA_HLAmatching.R, PRS_gen03_HLA_HLAmatching.R: calculate HLA PRS scores in both DeKAF and GEN-03 cohorts.
5. PRS_valid_HLAmatching.R, PRS_valid_HLAmatching_dekafsample.R: validate the PRS in both DeKAF and GEN-03 cohorts.
