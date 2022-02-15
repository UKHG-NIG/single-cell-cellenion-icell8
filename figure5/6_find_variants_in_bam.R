library(VariantAnnotation)
library(GenomicAlignments)
library(CogentDS)

# Deal with command line arguments.
require(optparse)
option_list <- list(
  # Mandatory options
  make_option(c("-r", "--Robj"), action = "store", default=NA, type="character",
              help = "R Object generated from CogentDS"),
  make_option(c("-b", "--bamFiles"), action = "store", default=NA, type="character",
              help = "BAM files containing all single-cell reads"),
  make_option(c("-v", "--VCF"), action = "store", default=NA, type="character",
              help = "Variant count file where to extract mutations from"),
  make_option(c("-c", "--coords"), action = "store", default=NA, type="character",
              help = "Coordinates to filter VCF for specific region"),
  make_option("--subSamples", action = "store", default=NA, type="character",
              help = "Samples used to plot tSNE"),
  make_option(c("-o", "--output"), action = "store", default=NA, type="character", dest="OutDir",
              help = "Output directory")
)
parser <- parse_args(OptionParser(option_list=option_list))
# check input file is provided
if (!is.na(parser$Robj)) {
  Robj <- parser$Robj
}else {
  stop("R object from CogentDS not defined. See script usage (--help)")
}
# check BAM file(s) is/are provided
if (!is.na(parser$bamFiles)) {
  bamFiles <- parser$bamFiles
  if(grepl(",",bamFiles)){
    bamFiles <- strsplit(bamFiles,split=",")[[1]]
  }
}else {
  stop("BAM files not defined. See script usage (--help)")
}
# check VCF file is provided
if (!is.na(parser$VCF)) {
  VCF <- parser$VCF
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
if(!is.na(parser$subSamples)) subSamples <- strsplit(parser$subSamples,",")[[1]]
print(parser)
  
# Load RDA object
load(Robj)
metadata <- CogentDS_data[["qc_data"]]$metadata
if(!is.na(parser$subSamples)) metadata <- CogentDS_data[["qc_data"]]$metadata[CogentDS_data[["qc_data"]]$metadata$Sample %in% subSamples,]
gm <- t(CogentDS_data[["pca_data"]]$pca_obj$x[, 1:CogentDS_data[["pca_data"]]$pc_count])
gm <- gm[,rownames(metadata)]
tsne_obj <- gm.tsne(gm = gm, metadata = metadata,
                    grouping_var = grouping_var, dims = 2, theta = 0.25,
                    max_iter = 1500, pca = FALSE, pca_center = FALSE, 
                    pca_scale = FALSE, perplexity = 20, plot = F)

# Extract all variants in VCF file in specific coordinates
compressVcf <- bgzip(VCF, tempfile())
tab <- indexVcf(compressVcf)
if(!is.na(parser$coords)){
  coords <- strsplit(parser$coords,",")[[1]]
  param = ScanVcfParam(
    which=GRanges(sapply(coords, function(coord) strsplit(coord,":")[[1]][1]),
                  IRanges(start=as.numeric(sapply(coords, function(coord) strsplit(strsplit(coord,":")[[1]][2],"-")[[1]][1])),
                          end=as.numeric(sapply(coords, function(coord) strsplit(strsplit(coord,":")[[1]][2],"-")[[1]][2])))))
  vcf <- readVcf(tab, "hg38", param)
}else{
  vcf <- readVcf(tab, "hg38")
}
for(mut in names(rowRanges(vcf))){
  metadata[[mut]] <- FALSE
  # metadata[[paste(mut,"REF",sep="_")]] <- 0
  metadata[[paste(mut,"ALT",sep="_")]] <- 0
  metadata[[paste(mut,"mutPerc",sep="_")]] <- 0
}

# Extract all reads associated with all VCF variants
param = ScanBamParam(which=rowRanges(vcf))
regions <- rowRanges(vcf)
for(i in 1:length(regions)){
  print(paste("Processing mutation",i,"of",length(regions)))
  region <- regions[i]
  if(length(bamFiles)>1){
    bam <- mergeBam(bamFiles, destination=tempfile(), region = region,
                    overwrite = T, header = character(), byQname = FALSE,
                    addRG = FALSE, compressLevel1 = FALSE, indexDestination = FALSE)
    bamIndex <- indexBam(bam)
  }else{
    bam <- bamFiles
  }
  galn <- readGappedReads(bam, param=param, use.names=T)
  mapRead <- mapToAlignments(rowRanges(vcf),galn)
  print(paste("Number of reads overlapping mutation is ",length(mapRead)))
  for(j in 1:length(mapRead)){
    map <- mapRead[j]
    mut <- rowRanges(vcf)[names(map)]
    read <- galn[map@seqnames@values]
    base <- as.character(qseq(read)[[1]][start(map)])
    barcode <- strsplit(names(read),"_")[[1]][2]
    if(!is.na(barcode)){
      metadata[metadata$Barcode==barcode,paste(names(mut),"ALT",sep="_")] <- metadata[metadata$Barcode==strsplit(names(read),"_")[[1]][2],paste(names(mut),"ALT",sep="_")]+1
      metadata[metadata$Barcode==barcode,names(mut)] <- TRUE
    }
  }
}
# Calculate mutation percentages
for(mut in names(rowRanges(vcf))){
  # metadata[,paste(mut,"mutPerc",sep="_")] <- metadata[,paste(mut,"ALT",sep="_")]/(metadata[,paste(mut,"ALT",sep="_")]+metadata[,paste(mut,"REF",sep="_")])*100
  metadata[is.na(metadata[,paste(mut,"mutPerc",sep="_")]),paste(mut,"mutPerc",sep="_")] <- 0
}
# Plot tSNE plots for each variant 
# width <- max_text_width(
#   names(rowRanges(vcf)), 
#   gp = gpar(fontsize = 12)
# )

mutPresenceCol <- structure(CogentDS:::add.alpha(colorRampPalette(takara_col(col = c("cellartis_red","light_gray")))(2), 
                                                 alpha = 1), names = c("TRUE","FALSE"))
for(coordinates in names(rowRanges(vcf))){
  if(sum(metadata[,coordinates])>0){
    pdf(file.path(OutDir,paste("Cells_withMut_",gsub("/","-",coordinates),".pdf", sep="")), width=12)
    print(reduction.plot(vis_obj = tsne_obj, metadata = metadata,
                         grouping_var = "Sample",
                         numeric_gradient_col=takara_col(col = c("light_gray", "cellartis_red")),
                         marker_size=0.5, alpha=1))
    print(reduction.plot(vis_obj = tsne_obj, metadata = metadata,
                         grouping_var = paste(coordinates,"ALT",sep="_"),
                         numeric_gradient_col=takara_col(col = c("light_gray", "cellartis_red")),
                         marker_size=0.5, alpha=1))
    print(reduction.plot(vis_obj = tsne_obj, metadata = metadata,
                         grouping_var = coordinates,col_override=mutPresenceCol,
                         marker_size=0.5, alpha=1))
    graphics.off()
  }
}

# Write out all cells containing at least one mutation
if(length(rowRanges(vcf))>1){
  df <- metadata[apply(metadata[,names(rowRanges(vcf))],1,any),c("Barcode","Sample",grep(paste(rowRanges(vcf),collapse="|"),colnames(metadata), value=T))]
}else{
  df <- metadata[metadata[,names(rowRanges(vcf))],c("Barcode","Sample",grep(paste(rowRanges(vcf),collapse="|"),colnames(metadata), value=T))]
}

write.table(df,file=file.path(OutDir,"Mutation_counts.csv"),
            col.names=T, row.names=F, quote=F, sep=",")

# Calculate mutation percentages per sample
mutPerc <- lapply(1:length(rowRanges(vcf)), function(i){
  mut <- names(rowRanges(vcf))[i]
  totalCells <- table(metadata$Sample)
  mutCells <- aggregate(x = metadata[[mut]],
            by = list(metadata$Sample),
            FUN = sum)
  percentages <- mutCells$x/totalCells*100
  df <- merge(merge(as.data.frame(totalCells), mutCells, by.x="Var1",by.y="Group.1"),as.data.frame(percentages), by="Var1")
  colnames(df) <- c("Sample","TotCells",paste("MutCells",mut,sep="-"),paste("MutantCellPerc",mut,sep="-"))
  df
})
mutPercDF <- Reduce(function(df1,df2) merge(df1,df2, by=c("Sample","TotCells")), mutPerc)
write.table(mutPercDF,file.path(OutDir,"MutationPerc_perSample.csv"),
            col.names=T, row.names=F, quote=F, sep=",")
