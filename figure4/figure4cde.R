# This scripts is used to generate tSNE plots for composite approach and ICELL8 platform
library(CogentDS)

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
composite_markers <- de.markers(gm = composite$qc_data$gm,
                  metadata =  composite$qc$metadata,
                  min_fc = 0,
                  grouping_var = "Graph-based_Clusters",
                  transform_type = composite$qc_data$params$gm_log_base)
for(i in 1:length(res)){
  colnames(composite_markers[[i]])[colnames(composite_markers[[i]])=="pval"] <- "pvalue"
  colnames(composite_markers[[i]])[colnames(composite_markers[[i]])=="log_fc"] <- "log2FoldChange"
  composite_markers[[i]]$rank <- -log10(composite_markers[[i]]$pvalue)*res[[i]]$log2FoldChange
  tmpRes <- composite_markers[[i]][order(composite_markers[[i]]$rank), ]
  gene_id <- sapply(rownames(tmpRes), function(gene) strsplit(x = gene, split="_")[[1]][1])
  Gene_Name <- sapply(rownames(tmpRes), function(gene) strsplit(x = gene, split="_")[[1]][2])
  composite_markers[[i]] <- as.data.frame(cbind(gene_id, Gene_Name, tmpRes))
}
heatmapGenes <- Reduce("union",lapply(composite_markers, function(obj){
  tmpGenes <- rownames(obj[order(obj$padj, decreasing = F),1:3][1:30,])
  sapply(tmpGenes, function(gene) strsplit(gene,"_")[[1]][2])
}))
mat <- as.matrix(vennValues[,grep("rank",colnames(vennValues), value=T)])
colnames(mat) <- gsub("rank_","",colnames(mat))
rownames(mat) <- vennValues$Gene_Name
mat[is.na(mat)] <- 0
matTopMarkers <- mat[rownames(mat) %in% heatmapGenes,]
plotHeight <- if(round(nrow(matTopMarkers)/4)<7){7}else{round(nrow(matTopMarkers)/4)}
pdf("Fig4E_heatmap_composite_clusters_ranks_top30MarkersPerCluster.pdf", height = plotHeight, width = 9)
print(Heatmap(scale(matTopMarkers), name = "ranks", column_title = "Cluster", row_title = "Gene",
              row_names_max_width = max_text_width(rownames(matTopMarkers), gp = gpar(fontsize = 12)),
              column_names_max_height = max_text_width(colnames(matTopMarkers), gp = gpar(fontsize = 12))))
graphics.off()

# Supp Figure 4a: find top 30 markers for each ICELL8 sample, and generate a heatmap similar to 4E
icell8_markers <- de.markers(gm = icell8$qc_data$gm,
                  metadata =  icell8$qc$metadata,
                  min_fc = 0,
                  grouping_var = "Sample",
                  transform_type = icell8$qc_data$params$gm_log_base)

for(i in 1:length(res)){
  colnames(icell8_markers[[i]])[colnames(icell8_markers[[i]])=="pval"] <- "pvalue"
  colnames(icell8_markers[[i]])[colnames(icell8_markers[[i]])=="log_fc"] <- "log2FoldChange"
  icell8_markers[[i]]$rank <- -log10(icell8_markers[[i]]$pvalue)*res[[i]]$log2FoldChange
  tmpRes <- icell8_markers[[i]][order(icell8_markers[[i]]$rank), ]
  gene_id <- sapply(rownames(tmpRes), function(gene) strsplit(x = gene, split="_")[[1]][1])
  Gene_Name <- sapply(rownames(tmpRes), function(gene) strsplit(x = gene, split="_")[[1]][2])
  res[[i]] <- as.data.frame(cbind(gene_id, Gene_Name, tmpRes))
}
heatmapGenes <- Reduce("union",lapply(res, function(obj){
  tmpGenes <- rownames(obj[order(obj$padj, decreasing = F),1:3][1:30,])
  sapply(tmpGenes, function(gene) strsplit(gene,"_")[[1]][2])
}))
mat <- as.matrix(vennValues[,grep("rank",colnames(vennValues), value=T)])
colnames(mat) <- gsub("rank_","",colnames(mat))
rownames(mat) <- vennValues$Gene_Name
mat[is.na(mat)] <- 0
matTopMarkers <- mat[rownames(mat) %in% heatmapGenes,]
plotHeight <- if(round(nrow(matTopMarkers)/4)<7){7}else{round(nrow(matTopMarkers)/4)}
pdf("Fig4E_heatmap_composite_clusters_ranks_top30MarkersPerCluster.pdf", height = plotHeight, width = 9)
print(Heatmap(scale(matTopMarkers), name = "ranks", column_title = "Cluster", row_title = "Gene",
              row_names_max_width = max_text_width(rownames(matTopMarkers), gp = gpar(fontsize = 12)),
              column_names_max_height = max_text_width(colnames(matTopMarkers), gp = gpar(fontsize = 12))))
graphics.off()

# Supp Figure 4b: Venn diagram of top 30 CellenONE and ICELL8 markers
