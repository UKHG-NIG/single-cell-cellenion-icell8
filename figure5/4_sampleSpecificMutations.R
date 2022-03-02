# Deal with command line arguments.
require(optparse)
option_list <- list(
  # Mandatory options
  make_option(c("-i", "--inFile"), action = "store", default=NA, type="character",
              help = "Comma-separated file with all variants in all samples"),
  make_option(c("-o", "--output"), action = "store", default=NA, type="character", dest="OutDir",
              help = "Output directory")
)
parser <- parse_args(OptionParser(option_list=option_list))
# check VCFs file are provided
if (!is.na(parser$inFile)) {
  inFile <- parser$inFile
}else {
  stop("Input file not defined. See script usage (--help)")
}
# Output directory
if (!is.na(parser$OutDir)) {
  OutDir <- parser$OutDir
  dir.create(OutDir,showWarnings = F)
}else {
  stop("Output directory not defined. See script usage (--help)")
}

tab <- read.csv(inFile, header = T)
samples <- c("Goe247","Goe486","Goe615","Goe800","Goe1309","Goe1360")

# Create table for icell8 variants
icell8_tab <- tab[,c(1:4,grep("icell8",colnames(tab), ignore.case = T))]
# Create table for icell8 variants
bulk_tab <- tab[,c(1:4,grep("bulk",colnames(tab)))]

sampleSpecific <- list()
for(sample in samples){ # Find mutations specific for each sample
  # Get columns with sample names
  notSamples <- paste(c("Goe1303","Goe1305",grep(sample,samples,invert=T, value=T)),collapse="|")
  sampleCols <- grep(sample,colnames(icell8_tab))
  nonSampleCols <- grep(notSamples,colnames(icell8_tab))
  mutRateCols <- grep(paste("mutRate",sample,sep="_"),colnames(icell8_tab))
  altDepthCols <- grep(paste("ALT",sample,sep="_"),colnames(icell8_tab))
  
  # Get all variations which are not NA in target sample
  icell8_sampleSpecific <- icell8_tab[apply(!is.na(icell8_tab[,sampleCols]),1,all),]
  # Remove mutations appearing in non-target sample
  icell8_sampleSpecific <- icell8_sampleSpecific[apply(is.na(icell8_sampleSpecific[,nonSampleCols]),1,all),]
  
  # Find mutations which are appearing in bulk as well
  bulk_sampleSpecific <- bulk_tab[rownames(icell8_sampleSpecific),]
  bulk_sampleSpecific <- bulk_sampleSpecific[apply(!is.na(bulk_sampleSpecific[,grep(sample,colnames(bulk_sampleSpecific))]),1,all),]
  
  tmp <- tab[rownames(bulk_sampleSpecific),]
  tmp <- tmp[order(tmp$chr, tmp$pos),]
  sampleSpecific[[sample]] <- tmp
  write.table(tmp,file.path(OutDir,paste("vars_sample",sample,".csv",sep="")), sep=",",row.names=F, col.names=T, quote = F)
}
# # Get depth statistics prior filtering
# statsDepth <- lapply(1:length(sampleSpecific), function(i){
#   obj <- sampleSpecific[[i]]
#   name <- names(sampleSpecific)[i]
#   apply(obj[,grep("genotype",grep(name, colnames(obj), value=T), value=T, invert=T)], 2,summary)
# })
# write.table(Reduce("cbind", statsDepth),file.path(OutDir,"depthStats_preFilt.csv"), sep=",",row.names=T, col.names=NA)

# Filter for depth
# 1) Alternative depth >= median alternative depth
# 2) Total depth >= (median reference depth + median alternative depth)
# 3) Mutation rate >= 10
sampleSpecific_postDepthFilt <- list()
for(i in 1:length(sampleSpecific)){
  obj <- sampleSpecific[[i]]
  sample <- names(sampleSpecific)[i]
  mediansRef <- statsDepth[[i]]["Median",grep("REF",colnames(statsDepth[[i]]))]
  mediansAlt <- statsDepth[[i]]["Median",grep("ALT",colnames(statsDepth[[i]]))]
  mediansMutRate <- statsDepth[[i]]["Median",grep("mutRate",colnames(statsDepth[[i]]))]
  obj <- obj[apply(obj[,names(mediansAlt)]>=mediansAlt,1,all),]
  obj <- obj[apply(obj[,names(mediansRef)]+obj[,names(mediansAlt)]>=mediansRef+mediansAlt,1,all),]
  obj <- obj[apply(obj[,names(mediansMutRate)]>=10,1,all),]
  sampleSpecific_postDepthFilt[[sample]] <- obj
  write.table(obj,file.path(OutDir,paste("vars_filtDepth_sample",sample,".csv",sep="")), sep=",",row.names=F, col.names=T, quote = F)
}

# # Get depth statistics post filtering
# statsDepth <- lapply(1:length(sampleSpecific_postDepthFilt), function(i){
#   obj <- sampleSpecific_postDepthFilt[[i]]
#   name <- names(sampleSpecific_postDepthFilt)[i]
#   apply(obj[,grep("genotype",grep(name, colnames(obj), value=T), value=T, invert=T)], 2,summary)
# })
# write.table(Reduce("cbind", statsDepth),file.path(OutDir,"depthStats_postFilt.csv"), sep=",",row.names=T, col.names=NA)

# # Get genotype statistics prior filtering
# statsDepth <- lapply(1:length(sampleSpecific_postDepthFilt), function(i){
#   obj <- sampleSpecific_postDepthFilt[[i]]
#   name <- names(sampleSpecific_postDepthFilt)[i]
#   as.data.frame(apply(obj[,grep("genotype",grep(name, colnames(obj), value=T), value=T)], 2,table))
# })
# write.table(Reduce("cbind", statsDepth),file.path(OutDir,"genoStats_preFilt.csv"), sep=",",row.names=T, col.names=NA)

# Filter for genotype (not 0/1)
for(i in 1:length(sampleSpecific_postDepthFilt)){
  obj <- sampleSpecific_postDepthFilt[[i]]
  sample <- names(sampleSpecific_postDepthFilt)[i]
  obj <- obj[apply(obj[,colnames(statsDepth[[i]])]!="0/1",1,all),]
  write.table(obj,file.path(OutDir,paste("vars_filtGeno_sample",sample,".csv",sep="")), sep=",",row.names=F, col.names=T, quote = F)
}
