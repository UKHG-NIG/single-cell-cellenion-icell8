library(CogentDS)
library(ComplexHeatmap)
library(circlize)
# Average normalised counts per sample
load("/zfstank/ngsdata/projects/icell8_cellenion_publication/figure4/CogentDS_samples_icell8_abslowcov5000/CogentDS.analysis.rda")
samples <- unique(CogentDS_data$qc_data$metadata$Sample)

sampleMat <- sapply(samples, function(sample){
  barcodes <- rownames(CogentDS_data$qc_data$metadata[CogentDS_data$qc_data$metadata$Sample %in% sample,])
  gm <- CogentDS_data$qc_data$gm[,barcodes]
  rowMeans(gm)
})
markers <- de.markers(gm=CogentDS_data$qc_data$gm, metadata=CogentDS_data$qc_data$metadata,
                      grouping_var="Sample",
                      transform_type=CogentDS_data$qc_data$params$transform_type)
markers <- lapply(markers, function(obj) rownames(obj[obj$padj<=0.05,]))

# Load mutant percentages
mutPercs <- read.table("/zfstank/ngsdata/projects/icell8_cellenion_publication/figure5/new_target_variants/MutationPerc_perSample.csv", sep=",", header=T)
mutationGenes <- read.table("/zfstank/ngsdata/projects/icell8_cellenion_publication/figure5/final_targets.csv", header=T, sep=",", comment.char="")
mutationGenes <- mutationGenes[,c("samples_icell8","genes","chr","pos")]

mutationTests <- list()
markersTmp <- Reduce("union", markers)
tmpSampleMat <- sampleMat[markersTmp,]
for(i in 1:nrow(mutationGenes)){#c(6,9)){
  mutation <- gsub("MutCells.","",grep(mutationGenes$pos[i],colnames(mutPercs), value=T)[1])
  print(paste("Run F-test for mutation",mutation))
  gene <- mutationGenes$genes[i]
  geneMat <- sapply(samples, function(sample){
    barcodes <- rownames(CogentDS_data$qc_data$metadata[CogentDS_data$qc_data$metadata$Sample %in% sample,])
    gm <- CogentDS_data$qc_data$gm[grep(gene, rownames(CogentDS_data$qc_data$gm)),barcodes]
    mean(gm)
  })
  df <- data.frame(t(tmpSampleMat),
                   mutPercs=mutPercs[match(colnames(tmpSampleMat),toupper(mutPercs$Sample)),grep(mutation,colnames(mutPercs))[2]],
                   geneAverage=geneMat[match(colnames(tmpSampleMat),names(geneMat))])
  
  ftests <- lapply(grep("ENS",colnames(df)), function(j){
    model <- lm(get(colnames(df)[j]) ~ mutPercs+geneAverage, data = df)
    tmp <- summary(model)$fstatistic
    pval <- pf(summary(model)$fstatistic[1],summary(model)$fstatistic[2],summary(model)$fstatistic[3])[[1]]
    data.frame(F=tmp[[1]], df1=tmp[[2]],df2=tmp[[3]], pvalueF=pval,pvalMutPerc=summary(model)$coefficients["mutPercs",4])
  })
  ftests <- as.data.frame(Reduce("rbind",ftests))
  rownames(ftests) <- grep("ENS",colnames(df), value=T)
  ftests$padjF <- p.adjust(ftests$pvalueF, method = "BH")
  ftests$padjMutPerc <- p.adjust(ftests$pvalMutPerc, method = "BH")
  ftests <- ftests[order(ftests$pvalMutPerc, decreasing = F),]
  mutationTests[[mutation]] <- ftests
}
OutDir <- "/zfstank/ngsdata/projects/icell8_cellenion_publication/figure5/new_target_variants/mutation_regressions"
dir.create(OutDir,showWarnings = F)
for(i in 1:length(mutationTests)){
  obj <- mutationTests[[i]]
  name <- names(mutationTests)[i]
  write.table(obj, file.path(OutDir, paste("F-test_sample",name,".csv",sep="")),
              row.names = T, col.names=NA, sep=",", quote=F)
}

# Heatmap of top regressed genes ####
samples <- as.factor(CogentDS_data$qc_data$metadata$Sample)
col <- structure(as.vector(takara_col(length(levels(samples)),alpha = 0.5)), names = levels(samples))
for(i in 1:length(mutationTests)){
  # Merge gene p-values from model test to gene matrix
  gm <- CogentDS_data$qc_data$gm[markersTmp,]
  gm <- merge(gm,mutationTests[[i]][,c("pvalMutPerc","padjMutPerc")], by="row.names")
  rownames(gm) <- gm[,1]
  gm <- gm[order(gm$pvalMutPerc, decreasing=F),]
  plotMutPerc <- gm[,grep("MutPerc",colnames(gm))]
  gm <- t(gm[,grep("_",colnames(gm))])
  # Merge Sample information and mutation percetage to gene matrix
  metadata <- CogentDS_data$qc_data$metadata[,c("Barcode","Sample")]
  gm <- merge(gm, metadata, all=T, by="row.names")
  rownames(gm) <- gm$Row.names
  gm <- merge(gm, mutPercs[,c(1,grep(names(mutationTests)[i],colnames(mutPercs))[2])], by="Sample", all.x=T)
  rownames(gm) <- gm$Row.names
  gm <- gm[order(gm[,ncol(gm)], rownames(gm), decreasing = F),]
  
  plotMutPerc <- plotMutPerc[1:30,]
  gm <- gm[,colnames(gm) %in% c("Sample",rownames(plotMutPerc),colnames(gm[ncol(gm)]))]
  rownames(plotMutPerc) <- sapply(rownames(plotMutPerc), function(gene) strsplit(gene,"_")[[1]][2])
  colnames(gm)[grep("ENS",colnames(gm))] <- sapply(colnames(gm)[grep("ENS",colnames(gm))], function(gene) strsplit(gene,"_")[[1]][2])
  
  mat <- t(gm[,-c(1,ncol(gm))])
  dend1 <- cluster_between_groups(mat,gm$Sample)
  
  fa <- order(gm[,ncol(gm)],gm$Sample, rownames(gm), decreasing=F)
  
  row_ha = rowAnnotation(MutPerc=anno_barplot(gm[,ncol(gm)], ylim=c(0,max(gm[,ncol(gm)]+0.1)), gp=gpar(fill=col)), width=unit(3,"cm"))#,Sample=gm$Sample)
  pvalue = plotMutPerc$pvalMutPerc
  padj = plotMutPerc$padjMutPerc
  is_sig = padj < 0.05
  pch = rep("*", length(is_sig))
  pch[!is_sig] = NA
  
  # color mapping for -log10(pvalue)
  pvalue_col_fun = colorRamp2(c(0, 2, 3), c("green", "white", "brown4")) 
  col_ha = HeatmapAnnotation(
    pvalue = anno_simple(-log10(pvalue), pt_gp = gpar(col = "white"), col = pvalue_col_fun, pch = pch),
    annotation_name_side = "left")
  ht <- Heatmap(scale(t(as.matrix(mat))), name = "NormExp",
                show_row_names=F,show_column_names=T,height = unit(200, "mm"),
                # cluster_rows = dend1,
                cluster_columns = F,
                column_names_side = "top",
                # row_dend_reorder = TRUE,cluster_rows = dend1,
                row_order = fa, row_split= gm$Sample,
                border = TRUE,row_gap = unit(0.5, "mm"),
                row_names_max_width = max_text_width(colnames(mat), gp = gpar(fontsize = 12)),
                column_names_max_height = max_text_width(rownames(mat), gp = gpar(fontsize = 12)),#,
                left_annotation = row_ha, top_annotation = col_ha)
  # now we generate two legends, one for the p-value
  # see how we define the legend for pvalue
  lgd = Legend(title = "p-value", col_fun = pvalue_col_fun, at = c(0, 1, 2, 3),
               labels = c("1", "0.1", "0.01", "0.001"))
  # and one for the significant p-values
  lgd_sig = Legend(pch = "*", type = "points", labels = "padj < 0.05")
  # these two self-defined legends are added to the plot by `annotation_legend_list`
  
  pdf(file.path(OutDir,paste("Heatmap_regressions_mutation_",names(mutationTests[i]),".pdf",sep="")), height = 10, width = 9)
  draw(ht, annotation_legend_list = list(lgd, lgd_sig))
  graphics.off()
  }