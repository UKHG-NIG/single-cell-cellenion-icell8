# This scripts is used to generate tSNE plots for composite approach and ICELL8 platform
library(CogentDS)
FDR=0.05
FC=0.5

composite <- get(load("CogentDS_CellenONE_fig4/CogentDS.analysis.rda"))
icell8 <- get(load("CogentDS_ICELL8_fig4/CogentDS.analysis.rda"))

# Figure 4C: tSNE plot of composite approach with graph-based clusters ####
reduction.plot(vis_obj = composite[["vis_data"]], metadata = composite[["qc_data"]]$metadata, 
               grouping_var = "Graph-based_Clusters", marker_size=0.5, alpha=1, file="Fig4C_composite_clusters_tSNE.pdf")

# Figure 4D: tSNE plots for specific markers in data ####
reduction.plot(vis_obj = composite[["vis_data"]], gm=composite[["qc_data"]]$gm,
               metadata = composite[["qc_data"]]$metadata, 
               gene_highlight = "ENSG00000087245_MMP2", file="Fig4D_composite_MMP2_tSNE.pdf", marker_size=0.5, alpha=1)
reduction.plot(vis_obj = composite[["vis_data"]], gm=composite[["qc_data"]]$gm,
               metadata = composite[["qc_data"]]$metadata, 
               gene_highlight = "ENSG00000188486_H2AX", file="Fig4D_composite_H2AX_tSNE.pdf", marker_size=0.5, alpha=1)
reduction.plot(vis_obj = composite[["vis_data"]], gm=composite[["qc_data"]]$gm,
               metadata = composite[["qc_data"]]$metadata, 
               gene_highlight = "ENSG00000184752_NDUFA12", file="Fig4D_composite_NDUFA12_tSNE.pdf", marker_size=0.5, alpha=1)

# Figure 4E: heatmap of top 30 markers based on their ranks ####
# Obtain markers for each cluster in the composite approach
composite_markers <- de.markers(gm = composite$qc_data$gm,
                                metadata =  composite$qc$metadata,
                                min_fc = 0,
                                grouping_var = "Graph-based_Clusters",
                                transform_type = composite$qc_data$params$gm_log_base)
# Create tables of markers tested for each cluster
for(i in 1:length(composite_markers)){
  sample <- names(composite_markers)[i]
  obj <- composite_markers[[i]]
  colnames(obj)[colnames(obj)=="pval"] <- "pvalue"
  colnames(obj)[colnames(obj)=="log_fc"] <- "log2FoldChange"
  obj$pvalue[obj$pvalue==0] <- .Machine$double.xmin
  obj$rank <- -log10(obj$pvalue)*obj$log2FoldChange
  tmpRes <- obj[order(obj$rank), ]
  gene_id <- sapply(rownames(tmpRes), function(gene) strsplit(x = gene, split="_")[[1]][1])
  Gene_Name <- sapply(rownames(tmpRes), function(gene) strsplit(x = gene, split="_")[[1]][2])
  obj <- as.data.frame(cbind(gene_id, Gene_Name, tmpRes))
  composite_markers[[i]] <- obj
  write.table(composite_markers[[i]], file = paste("markers_composite_clusters_",sample,".csv", sep=""),
              sep=",", quote = F, row.names=T, col.names = NA)
}
# Retrieve top 30 markers with smallest adjusted p-value in each cluster
heatmapGenes <- Reduce("union",lapply(composite_markers, function(obj){
  tmpGenes <- rownames(obj[order(obj$padj, decreasing = F),][1:30,])
  sapply(tmpGenes, function(gene) strsplit(gene,"_")[[1]][2])
}))
# Retrieve deregulated markers (adjusted p-value <=0.05 and absolute log2 fold-change >1) for each cluster
deGenes <- Reduce("union",lapply(composite_markers, function(obj){
  obj <- obj[abs(obj$log2FoldChange)>FC & !is.na(obj$padj) & obj$padj<=FDR,]
  obj$Gene_Name
}))
# Obtain deregulated markers also in the top30 list
deTop30 <- intersect(heatmapGenes, deGenes)
# Get ranks for each of the clusters
composite_matTopMarkers <-  matrix(0, nrow=length(deTop30), ncol=length(composite_markers), dimnames=list(deTop30,names(composite_markers)))
for(i in 1:length(composite_markers)){
  obj <- composite_markers[[colnames(composite_matTopMarkers)[i]]]
  composite_matTopMarkers[,i] <- obj$rank[match(rownames(composite_matTopMarkers), obj$Gene_Name)]
}
# Plot ranks of target genes
plotHeight <- if(round(nrow(composite_matTopMarkers)/4)<7){7}else{round(nrow(composite_matTopMarkers)/4)}
pdf("Fig4E_heatmap_composite_clusters_ranks_top30MarkersPerCluster.pdf", height = plotHeight, width = 9)
print(Heatmap(scale(composite_matTopMarkers), name = "ranks", column_title = "Cluster", row_title = "Gene",
              row_names_max_width = max_text_width(rownames(composite_matTopMarkers), gp = gpar(fontsize = 12)),
              column_names_max_height = max_text_width(colnames(composite_matTopMarkers), gp = gpar(fontsize = 12))))
graphics.off()

# Supp Figure 4a: find top 30 markers for each ICELL8 sample, and generate a heatmap similar to 4E
# Obtain markers for each sample in the ICELL8 dataset
icell8_markers <- de.markers(gm = icell8$qc_data$gm,
                                metadata =  icell8$qc$metadata,
                                min_fc = 0,
                                grouping_var = "Sample",
                                transform_type = icell8$qc_data$params$gm_log_base)
# Create tables of markers tested for each sample
for(i in 1:length(icell8_markers)){
  sample <- names(icell8_markers)[i]
  obj <- icell8_markers[[i]]
  colnames(obj)[colnames(obj)=="pval"] <- "pvalue"
  colnames(obj)[colnames(obj)=="log_fc"] <- "log2FoldChange"
  obj$pvalue[obj$pvalue==0] <- .Machine$double.xmin
  obj$rank <- -log10(obj$pvalue)*obj$log2FoldChange
  tmpRes <- obj[order(obj$rank), ]
  gene_id <- sapply(rownames(tmpRes), function(gene) strsplit(x = gene, split="_")[[1]][1])
  Gene_Name <- sapply(rownames(tmpRes), function(gene) strsplit(x = gene, split="_")[[1]][2])
  obj <- as.data.frame(cbind(gene_id, Gene_Name, tmpRes))
  icell8_markers[[i]] <- obj
  write.table(icell8_markers[[i]], file = paste("markers_icell8_samples_",sample,".csv", sep=""),
              sep=",", quote = F, row.names=T, col.names = NA)
}
# Retrieve top 30 markers with smallest adjusted p-value in each sample
heatmapGenes <- Reduce("union",lapply(icell8_markers, function(obj){
  tmpGenes <- rownames(obj[order(obj$padj, decreasing = F),][1:30,])
  sapply(tmpGenes, function(gene) strsplit(gene,"_")[[1]][2])
}))
# Retrieve deregulated markers (adjusted p-value <=0.05 and absolute log2 fold-change >1) for each sample
deGenes <- Reduce("union",lapply(icell8_markers, function(obj){
  obj <- obj[abs(obj$log2FoldChange)>FC & !is.na(obj$padj) & obj$padj<=FDR,]
  obj$Gene_Name
}))
# Obtain deregulated markers also in the top30 list
deTop30 <- intersect(heatmapGenes, deGenes)
# Get ranks for each of the samples
icell8_matTopMarkers <-  matrix(0, nrow=length(deTop30), ncol=length(icell8_markers), dimnames=list(deTop30,names(icell8_markers)))
for(i in 1:length(icell8_markers)){
  obj <- icell8_markers[[colnames(icell8_matTopMarkers)[i]]]
  icell8_matTopMarkers[,i] <- obj$rank[match(rownames(icell8_matTopMarkers), obj$Gene_Name)]
}
# Plot ranks of target genes
plotHeight <- if(round(nrow(icell8_matTopMarkers)/4)<7){7}else{round(nrow(icell8_matTopMarkers)/4)}
pdf("SuppFig4a_heatmap_icell8_samples_ranks_top30MarkersPerCluster.pdf", height = plotHeight, width = 9)
print(Heatmap(scale(icell8_matTopMarkers), name = "ranks", column_title = "Sample", row_title = "Gene",
              row_names_max_width = max_text_width(rownames(icell8_matTopMarkers), gp = gpar(fontsize = 12)),
              column_names_max_height = max_text_width(colnames(icell8_matTopMarkers), gp = gpar(fontsize = 12))))
graphics.off()

# Supp Figure 4b: Venn diagram of top 30 CellenONE and ICELL8 markers
require(VennDiagram)
# Plot Venn diagram overlapping markers from both datasets
pdf("SuppFig4b_Venn.pdf")
grid.draw(venn.diagram(x = list("CellenONE-ICELL8"=rownames(composite_matTopMarkers),
                                "ICELL8"=rownames(icell8_matTopMarkers)), col = "transparent",
                       fontfamily ="sans",cat.fontfamily="sans",
                       fill = c("blue", "yellow"),# fill = c("blue", "green","red"),
                       alpha = 0.5, cat.cex=1.5, cex=1.5, height=100, width=100,filename=NULL,
                       print.mode=c("raw","percent")))
graphics.off()
