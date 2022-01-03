# Regulon analysis
This directory contains scripts used to perform the variant calling in scRNA- and bulk RNA-seq BAM files, and use the resulting VCF files to perform the regulon analysis. The analysis steps are as follows:
- 1) For the scRNA-seq BAM file, extract reads associated to different cells based on their sample affiliation, and build BAM files for each sample (in our case, we use the "grep" command in shell; you are welcome to try a different method if you know a faster one to do it)
- 2) For the sample-specific BAM files from the scRNA-seq and the bulk RNA-seq data, call for the variants in each file using GATK
- 3) Combine the VCF files from the different samples using ```compare_VCFs.R```
- 4) Find sample-specific mutations using ```sampleSpecificMutations.R```
- 5) Determine which mutations are non-synonymous using ```getCodonInfo_target_variants.R```
- 6) Calculate mutation percentages for top variants using ```find_variants_in_bam.R```
- 7) Calculate regressions for target variations using ```regression_mutPerc-vs-averageCounts.R```
