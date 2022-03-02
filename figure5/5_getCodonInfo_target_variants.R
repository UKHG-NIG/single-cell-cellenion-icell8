library(VariantAnnotation)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Deal with command line arguments.
require(optparse)
option_list <- list(
  # Mandatory options
  make_option(c("-v", "--VCF"), action = "store", default=NA, type="character",
              help = "Comma-separated VCFs to compare"),
  make_option(c("-t", "--targetVars"), action = "store", default=NA, type="character",
              help = "Path to directory with sample-specific variants filtered out for low-depth variants"),
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
# check target variant path is set
if (!is.na(parser$targetVars)) {
  targetVars <- list.files(parser$targetVars, pattern="filtDepth", full.names=T)
}else {
  stop("Target variant path not defined. See script usage (--help)")
}
# Output directory
if (!is.na(parser$OutDir)) {
  OutDir <- parser$OutDir
  dir.create(OutDir,showWarnings = F)
}else {
  stop("Output directory not defined. See script usage (--help)")
}
vcfNames <- sapply(basename(VCFfiles),function(name) gsub("_variants_filt.vcf","",name))

# Extract all variants in VCF file in specific coordinates
vcf <- lapply(1:length(VCFfiles), function(i){
  vcfName <- vcfNames[i]
  VCF <- VCFfiles[i]
  compressVcf <- bgzip(VCF, tempfile())
  tab <- indexVcf(compressVcf)
  vcf <- readVcf(tab, "hg38")
  if(length(vcf)==0){
    print(paste("No mutation found for sample",vcfName))
    return(NULL)
  }
  ## Retrieve depth-filtered targets
  targetFile <- grep(strsplit(vcfName, "_")[[1]][1],targetVars, value=T)
  if(length(targetFile)==0){
    print(paste("No targets found for sample",vcfName))
    return(NULL)
  }
  targets <- read.table(targetFile, sep=",", header=T)
  targets <- sapply(1:nrow(targets), function(i) paste(targets$chr[i],":",targets$pos[i],"_",targets$REF[i],"/",targets$ALT[i],sep=""))
  
  vcf <- vcf[targets]
  ## Rename seqlevels with renameSeqlevesl()
  vcf <- renameSeqlevels(vcf, paste0("chr", seqlevels(vcf)))
  vcf <- renameSeqlevels(vcf, gsub("chrMT","chrM", seqlevels(vcf)))
  vcf <- keepSeqlevels(vcf, unique(seqnames(vcf)), pruning.mode = "coarse")
  ## Create data.frame with codon predictions
  mutations <- predictCoding(vcf, txdb, Hsapiens)
  if(length(mutations)==0){
    print(paste("No predicted codons for mutations in sample",vcfName))
    return(NULL)
  }
  mutations <- unique(mutations[,c("REF","ALT","varAllele","CDSLOC","PROTEINLOC","CONSEQUENCE","REFCODON","VARCODON","REFAA","VARAA")])
  
  ## Add genotype information
  mutations <- unique(merge(as.data.frame(mutations),geno(vcf)$GT, by="row.names"))
  rownames(mutations) <- mutations$Row.names
  mutations <- mutations[,grep("Row.names|end|width|CDSLOC.end|CDSLOC.width",colnames(mutations), invert=T)]
  colnames(mutations) <- c("chr","pos","strand","REF","ALT","varAllele","CDSLOC","PROTEINLOC","CONSEQUENCE","REFCODON","VARCODON","REFAA","VARAA","genotype")
  
  ## Add depth information
  tmp <- data.frame("REF"=sapply(geno(vcf)$AD,function(obj) obj[1]),
                    "ALT"=sapply(geno(vcf)$AD,function(obj) obj[2]))
  rownames(tmp) <- rownames(vcf)
  colnames(tmp) <- c("#REF","#ALT")
  mutations <- unique(merge(mutations, tmp, by="row.names"))
  rownames(mutations) <- mutations$Row.names
  mutations <- mutations[,-1]
  mutations$mutRate <- 100*mutations$`#ALT`/(mutations$`#ALT`+mutations$`#REF`)
  colnames(mutations)[5:ncol(mutations)] <- paste(colnames(mutations)[5:ncol(mutations)],vcfName,sep="_")
  return(mutations)
})
names(vcf) <- vcfNames
vcf <- vcf[lapply(vcf,length)>0] ## you can use sapply,rapply

# Reduce list of targets into data.frame
df <- Reduce(function(df1, df2){merge(df1, df2, all=T)}, vcf)
# Remove variations in control (GOE1303 and GOE1305) but not in disease patients
disease <- grep("1303|1305", grep("genotype*.*ICELL8",colnames(df), value=T), invert=T, value=T)
df <- df[rowSums(is.na(df[,disease]))<6,]
# Obtain samples where variants are expressed
## ICELL8 samples
df$samples_icell8 <- sapply(1:nrow(df), function(i){
  row <- df[i,]
  paste(gsub("genotype_|_ICELL8","",grep("genotype*.*ICELL8",colnames(row[,!is.na(row)]), value=T)),collapse=",")
})
## Bulk samples
df$samples_bulk <- sapply(1:nrow(df), function(i){
  row <- df[i,]
  paste(unique(gsub("genotype_|_bulk[123]","",grep("genotype*.*bulk*",colnames(row[,!is.na(row)]), value=T))),collapse=",")
})
# Remove variations in ICELL8 only and without bulk
df <- df[df$samples_bulk!="",]
df <- df[,grep("bulk[123]", colnames(df), invert = T)]
colnames(df) <- gsub("_ICELL8","",colnames(df))
# Annotate table
attributes <- "external_gene_name"
filters <- c("chromosome_name","start","end")
df$genes <- NA
for(i in 1:nrow(df)){
  values <- list(chromosome_name=gsub("chr","",gsub("chrM","chrMT",df$chr[i])),start=df$pos[i],end=df$pos[i])
  df$genes[i] <- paste(unique(getBM(attributes=attributes, filters=filters, values=values, mart=mart)),collapse=",")
}

write.table(df,file=file.path(OutDir, "VCF_annotations.csv"),
            col.names=T, row.names=F, quote=F, sep="\t")
