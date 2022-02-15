# library(VennDiagram)
library(VariantAnnotation)
library(stringr)

# Deal with command line arguments.
require(optparse)
option_list <- list(
  # Mandatory options
  make_option(c("-v", "--VCF"), action = "store", default=NA, type="character",
              help = "Comma-separated VCFs to compare"),
  make_option(c("-o", "--output"), action = "store", default=NA, type="character", dest="OutDir",
              help = "Output directory")
)
parser <- parse_args(OptionParser(option_list=option_list))
# check VCFs file are provided
if (!is.na(parser$VCF)) {
  VCFfiles <- strsplit(parser$VCF,split=",")[[1]]
}else {
  stop("VCF file not defined. See script usage (--help)")
}
# Output directory
if (!is.na(parser$OutDir)) {
  OutDir <- parser$OutDir
  dir.create(OutDir,showWarnings = F)
}else {
  stop("Output directory not defined. See script usage (--help)")
}
vcfNames <- sapply(basename(VCFfiles),function(name) gsub("_variants_filt.vcf","",name))

# Extract all variants in VCF file
vcf <- lapply(1:length(VCFfiles), function(i){
  # Read in VCF file
  vcfName <- vcfNames[i]
  VCF <- VCFfiles[i]
  compressVcf <- bgzip(VCF, tempfile())
  tab <- indexVcf(compressVcf)
  vcf <- readVcf(tab, "hg38")
  
  # Get coverage for reference and alternative nucleotides, including mutation rate (proportion of alternative reads out of total reads)
  mutNames <- names(rowRanges(vcf))
  mutations <- as.data.frame(t(sapply(mutNames, function(x) strsplit(x,":|_|/")[[1]])))
  colnames(mutations) <- c("chr","pos","REF","ALT")
  mutations[[paste("genotype",vcfName,sep="_")]] <- sapply(geno(vcf)$GT, function(mut) mut[[1]])
  mutations[[paste("#REF",vcfName,sep="_")]] <- sapply(geno(vcf)$AD, function(mut) mut[1])
  mutations[[paste("#ALT",vcfName,sep="_")]] <- sapply(geno(vcf)$AD, function(mut) mut[2])
  mutations[[paste("mutRate",vcfName,sep="_")]] <- 100*mutations[[paste("#ALT",vcfName,sep="_")]]/(mutations[[paste("#ALT",vcfName,sep="_")]]+mutations[[paste("#REF",vcfName,sep="_")]])
  return(mutations)
})
df <- Reduce(function(df1, df2){merge(df1, df2, all=T)}, vcf)

# Write out all variants found in all samples
write.table(df,file=file.path(OutDir,"VCF_annotations.csv"),
            col.names=T, row.names=F, quote=F, sep=",")
write.table(df[apply(df,1,function(row){all(!is.na(row))}),],file=file.path(OutDir,"VCF_annotations_all.csv"),
            col.names=T, row.names=F, quote=F, sep=",")
