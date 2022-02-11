# Regulon analysis
This directory contains scripts used to perform the variant calling in scRNA- and bulk RNA-seq BAM files, and use the resulting VCF files to perform the regulon analysis. The analysis steps are as follows:
- 1) For the scRNA-seq BAM file, extract reads associated to different cells based on their sample affiliation, and build BAM files for each sample (in our case, we use the "grep" command in shell; )
- 2) For the sample-specific BAM files from the scRNA-seq and the bulk RNA-seq data, call for the variants in each file using `GATK`. Specifically, the following commands were used to call the variants from each sample:
```Shell
# BAM to grouped/sorted BAM
gatk AddOrReplaceReadGroups -I /path/to/sample.bam -O /path/to/sample.group.bam -SO coordinate --CREATE_INDEX true --RGID RNA --RGLB RNA --RGPL Illumina --RGPU Hiseq --RGSM sample --TMP_DIR .

# Mark duplicates
gatk MarkDuplicates -I /path/to/sample.group.bam -O /path/to/sample.marked.bam --METRICS_FILE /path/to/sample_reads.metrics --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT

# CIGAR reads
gatk SplitNCigarReads -R /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I /path/to/sample.marked.bam -O /path/to/sample.split.bam --tmp-dir .

# Variant calling
gatk HaplotypeCaller -R /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa -I /path/to/sample.split.bam -O /path/to/sample_variants.vcf

# Variant filtering
gatk VariantFiltration -R /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa -V /path/to/sample_variants.vcf -O /path/to/sample_variants_filt.vcf --window 35 --cluster 3 --filter-name FS --filter "FS > 30.0" --filter-name QD --filter "QD < 2.0"
```
- 3) Combine the VCF files from the different samples using ```combine_VCFs.R```
- 4) Find sample-specific mutations using ```sampleSpecificMutations.R```
- 5) Determine which mutations are non-synonymous using ```getCodonInfo_target_variants.R```
- 6) Calculate mutation percentages for top variants using ```find_variants_in_bam.R```
- 7) Calculate regressions for target variations using ```regression_mutPerc-vs-averageCounts.R```
