# This script is used to generate the Venn diagrams of differentially expressed markers (adjusted p-value <=0.05) for each sample in both novel and ICELL8 approaches (Fig3C and Supp. Fig4)

# Load library
library(VennDiagram)

# Set adjusted p-value parameter
padj <- 0.05

samples <- c("GOE247","GOE486","GOE615","GOE800","GOE1303","GOE1305","GOE1309","GOE1360")
for(sample in samples){
  # Upload tables for markers 
  icell8_markerTable <- read.table(paste("markers_icell8_",sample,".csv",sep=""), row.names = 1, header=T, sep=",")
  composite_markerTable <- read.table(paste("markers_composite_",sample,".csv",sep=""), row.names = 1, header=T, sep=",")
  
  # Find deregulated markers for each 
  icell8_markers <- rownames(icell8_markerTable)[!is.na(icell8_markerTable$padj) & icell8_markerTable$padj<=padj]
  composite_markers <- rownames(composite_markerTable)[!is.na(composite_markerTable$padj) & composite_markerTable$padj<=padj]
    
  # Plot venn diagrams
  vennList=list()
  vennList$composite <- composite_markers
  vennList$icell8 <- icell8_markers
  venn.diagram(x = vennList, col = "transparent",
               fill = c("blue", "yellow"),# fill = c("blue", "green","red"),
               alpha = 0.5,cex=2,cat.cex=1.5,
               category.names = c("CellenONE-ICELL8","ICELL8"),
               print.mode=c("raw","percent"),
               filename = paste("Fig3c_Venn_markers_",sample,".pdf",sep=""))
}
