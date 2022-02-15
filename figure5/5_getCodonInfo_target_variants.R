library(VariantAnnotation)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(biomaRt)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
snps <- SNPlocs.Hsapiens.dbSNP151.GRCh38

mart <- useDataset("hsapiens_gene_ensembl", useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl"))

VCFfiles <- list.files(path="/zfstank/ngsdata/projects/icell8_cellenion_publication/figure5/new_target_variants", pattern = "_variants_filt.vcf$", full.names = T)
vcfNames <- sapply(VCFfiles, function(file) gsub("_variants_filt.vcf","",basename(file)))

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
  mutations$REF <- as.character(mutations$REF)
  mutations$ALT <- as.character(unlist(mutations$ALT))
  mutations$varAllele <- as.character(mutations$varAllele)
  mutations$CDSLOC <- start(mutations$CDSLOC)
  mutations$PROTEINLOC <- sapply(mutations$PROTEINLOC, function(loc) paste(as.character(loc),collapse=","))
  mutations$REFCODON <- as.character(mutations$REFCODON)
  mutations$VARCODON <- as.character(mutations$VARCODON)
  mutations$REFAA <- AMINO_ACID_CODE[unlist(as.character(mutations$REFAA))]
  mutations$VARAA <- AMINO_ACID_CODE[unlist(as.character(mutations$VARAA))]
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

## you can use sapply,rapply
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
  print(i)
  values <- list(chromosome_name=gsub("chr","",gsub("chrM","chrMT",df$chr[i])),start=df$pos[i],end=df$pos[i])
  df$genes[i] <- paste(unique(getBM(attributes=attributes, filters=filters, values=values, mart=mart)),collapse=",")
}
df$dbSNP <- NA
for(i in 1:nrow(df)){
  print(i)
  values <- GRanges(gsub("chr","",gsub("chrM","chrMT",df$chr[i])), IRanges(start=df$pos[i],end=df$pos[i]))
  value <- snpsByOverlaps(snps, values)
  if(length(value)>0)
    df$dbSNP[i] <- value$RefSNP_id
}

write.table(df,file="/zfstank/ngsdata/projects/icell8_cellenion_publication/figure5/new_target_variants/VCF_annotations.csv",
            col.names=T, row.names=F, quote=F, sep="\t")
