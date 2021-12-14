# This script is used to plot the normalised expressions of the top 100 markers (Figure 3b and Supp Fig 3)

library(DESeq2)
library(GGally)
library(ggpubr)

# Load CogentDS objects for composite and ICELL8 datasets
composite <- get(load("CogentDS_CellenONE/CogentDS.analysis.rda"))
icell8 <- get(load("CogentDS_ICELL8/CogentDS.analysis.rda"))

# Generate pseudocounts for icell8
icell8_samples <- unique(grep(pattern="_Ctrl", icell8$qc_data$metadata$Sample, invert=T, value=T))
icell8_pseudoGM <- sapply(icell8_samples, function(sample){
  indices <- icell8$raw_data$metadata$Sample==sample
  rowSums(as.matrix(icell8$raw_data$gm)[,indices])
})
dds_icell8 <- DESeqDataSetFromMatrix(icell8_pseudoGM, DataFrame(cond=colnames(icell8_pseudoGM)), ~ cond)
dds_icell8 <- estimateSizeFactors(dds_icell8)
icell8_norm <- as.data.frame(counts(dds_icell8, normalized=T))
colnames(icell8_norm) <- paste("ICELL8",toupper(colnames(icell8_norm)),sep="_")
rownames(icell8_norm) <- sapply(rownames(icell8_norm), function(gene) strsplit(x = gene,split = "_")[[1]][1])

# Generate pseudocounts for cellenion
cellenion_samples <- unique(grep(pattern="Goe247_wo__fix", composite$qc_data$metadata$Sample, invert=T, value=T))
cellenion_pseudoGM <- sapply(unique(cellenion_samples), function(sample){
  indices <- composite$raw_data$metadata$Sample==sample
  rowSums(as.matrix(composite$raw_data$gm)[,indices])
})
dds_cellenion <- DESeqDataSetFromMatrix(cellenion_pseudoGM, DataFrame(cond=colnames(cellenion_pseudoGM)), ~ cond)
dds_cellenion <- estimateSizeFactors(dds_cellenion)
cellenion_norm <- as.data.frame(counts(dds_cellenion, normalized=T))
colnames(cellenion_norm) <- paste("CellenONE",toupper(colnames(cellenion_norm)),sep="_")
rownames(cellenion_norm) <- sapply(rownames(cellenion_norm), function(gene) strsplit(x = gene,split = "_")[[1]][1])

# Load bulk data
pdata = read.table("bulk_table.csv",header=T,sep=',')
dds <- DESeqDataSetFromHTSeqCount(sampleTable = pdata,directory=".",design= ~group)
dds <- estimateSizeFactors(dds)
bulk_norm <- as.data.frame(counts(dds, normalized=T))

# Plot 
# Correlation panel
samples <- c("GOE1303","GOE1305","GOE1309","GOE1360","GOE247","GOE486","GOE615","GOE800")
for(sample in samples){
  # Load markers for each sample in each dataset
  df <- read.table(paste("top100_markers_",sample,".csv",sep=""),sep=",", header=T, row.names=1)
  markers_bulk <- sapply(rownames(df), function(gene) strsplit(gene, "_")[[1]][1])
  
  plotDF <- data.frame(icell8=icell8_norm[markers_bulk,grep(sample,colnames(icell8_norm))],
                       cellenion=cellenion_norm[markers_bulk,grep(sample,colnames(cellenion_norm))],
                       bulk_norm[markers_bulk,grep(sample,colnames(bulk_norm))])
  colnames(plotDF) <- gsub(paste(sample,"_",sep=""), "", colnames(plotDF))
  pdf(paste("Fig3b_scatterMatrix_",sample,".pdf",sep=""))
  print(ggpairs(plotDF, aes(alpha = 0.5),lower = list(continuous = "smooth"), axisLabels="show",
                diag = list(continuous = "blankDiag"),
                upper = list(continuous = wrap("cor", size = 7)))+
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major = element_blank(),strip.text.x=element_text(size=12),
                strip.text.y=element_text(size=12),
                axis.ticks = element_blank()))
  graphics.off()
}
