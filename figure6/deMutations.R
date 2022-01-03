require(circlize)
require(CogentDS)
require(ReactomePA)
require(Rsamtools)
require(GenomicAlignments)
require(ComplexHeatmap)
require(clusterProfiler)
require(DOSE)
require(org.Hs.eg.db)
require(GOSemSim)
require(WebGestaltR)
require(enrichplot)
setwd("/zfstank/ngsdata/projects/icell8_cellenion_publication/")
targets <- read.table("figure5/final_targets.csv",sep=",", comment.char = "", header=T)
load("figure4/CogentDS_samples_icell8_abslowcov5000/CogentDS.analysis.rda")
bamFiles <- list.files("figure5/new_target_variants", pattern="ICELL8_variants.bam", full.names = T)
markers <- list()
for(i in 1:nrow(targets)){
  sample <- targets$samples_icell8[i]
  chr <- gsub("M","MT",targets$chr[i])
  pos <- targets$pos[i]
  param = ScanBamParam(which=GRanges(seqnames=chr,
                                     IRanges(start=pos,end=pos),
                                     strand=targets$strand[i]))
  galn <- readGappedReads(grep(sample,bamFiles, value=T), param=param, use.names=T)
  cells <- sapply(names(galn), function(read) strsplit(read, "_")[[1]][2])
  cells <- cells[!is.na(cells)]
  mut <- paste(chr,":",pos,"_",sample,sep="")
  CogentDS_data$qc_data$metadata[[mut]] <- NA
  CogentDS_data$qc_data$metadata[[mut]][grep(sample, CogentDS_data$qc_data$metadata$Sample, ignore.case = T)] <- "not_mutated"
  CogentDS_data$qc_data$metadata[[mut]][CogentDS_data$qc_data$metadata$Barcode %in% cells] <- "mutated"

  metadata <- CogentDS_data$qc_data$metadata[!is.na(CogentDS_data$qc_data$metadata[[mut]]),]
  gm <- CogentDS_data$qc_data$gm[,rownames(metadata)]
  markers[[mut]] <- de.markers(gm=gm,min_fc = 0.25,
                               metadata=metadata,
                               grouping_var=mut,
                               transform_type = CogentDS_data$qc_data$params$transform_type)
  markers[[mut]] <- markers[[mut]]$mutated
  colnames(markers[[mut]]) <- c("baseMeanMutated","baseMeanNotMutated","log2FoldChange","pvalue","padj")
  markers[[mut]]$rank <- -log10(markers[[mut]]$pvalue)*markers[[mut]]$log2FoldChange
  markers[[mut]] <- markers[[mut]][order(markers[[mut]]$rank),]
  colnames(markers[[mut]]) <- paste(colnames(markers[[mut]]),sample,sep="_")
  markers[[mut]]$ENSEMBL <- sapply(rownames(markers[[mut]]), function(gene) strsplit(gene,"_")[[1]][1])
  markers[[mut]]$SYMBOL <- sapply(rownames(markers[[mut]]), function(gene) strsplit(gene,"_")[[1]][2])
}
for(i in 1:length(markers)){
  obj <- markers[[i]]
  name <- names(markers)[i]
  write.table(obj, file.path("figure6",paste("statResults_mut",name,".csv", sep="")), row.names = F, sep=",")
}
df <- Reduce(function(x,y) merge(x = x, y = y, all = TRUE, by=c("ENSEMBL","SYMBOL")), markers)

# Heatmap of DE genes in each test ####
FDR <- 0.05
mutations <- names(markers)
for(mut in mutations){
  name <- mut
  sample <- strsplit(mut, split="_")[[1]][2]
  obj <- markers[[mut]]
  genes <- rownames(obj[obj[[paste("padj",sample,sep="_")]]<=FDR,])
  if(length(genes)>0){
    # genes <- rownames(obj[obj[[paste("padj",sample,sep="_")]]<=FDR,])
    metadata <- CogentDS_data$qc_data$metadata[!is.na(CogentDS_data$qc_data$metadata[[mut]]),]
    gm <- expm1(CogentDS_data$qc_data$gm[genes,rownames(metadata)])
    rownames(gm) <- sapply(rownames(gm), function(gene) strsplit(gene, "_")[[1]][2])
    
    # library(dendextend)
    plotHeight <- if(round(nrow(gm)/4)<7){7}else{round(nrow(gm)/4)}
    
    cellMut <- rainbow( length(unique(metadata[[mut]]) ))
    names(cellMut) <- as.character(unique(metadata[[mut]]))
    
    lfc=obj[genes, grep("log2FoldChange",colnames(obj))]
    names(lfc) <- obj[genes,"SYMBOL"]
    lfc_col = colorRamp2(c(min(lfc),0, max(lfc)), c("red","white", "darkgreen"))
    padj=-log10(obj[genes, grep("padj",colnames(obj))])
    names(padj) <- obj[genes,"SYMBOL"]
    padj_col = colorRamp2(c(0, max(padj)), c("white", "darkgreen"))
    ranks=data.frame(rank=obj[genes, grep("rank",colnames(obj))])
    rownames(ranks) <- sapply(genes, function(gene) strsplit(gene,"_")[[1]][2])
    
    row_ha <- rowAnnotation(log2FC=lfc, "-log10(padj)"=padj,
                            col=list(log2FC=lfc_col,"-log10(padj)"=padj_col))
    
    col_ha = HeatmapAnnotation(df = data.frame(Mutation=metadata[[mut]]),
                               col = list(mut=cellMut))
    
    # pdf(file.path("figure6",paste("Heatmap_cells_",sample,".pdf",sep="")), height = plotHeight, width = 9)
    # print(Heatmap(scale(gm), name = "NormExp", column_title = "Cell", row_title = "Gene",
    #               height = unit(150*nrow(gm)/31, "mm"),show_column_names=F, show_row_names=T,
    #               top_annotation = col_ha, left_annotation=row_ha,column_split=metadata[[mut]],
    #               row_names_max_width = max_text_width(rownames(gm), gp = gpar(fontsize = 12)),
    #               column_names_max_height = max_text_width(colnames(gm), gp = gpar(fontsize = 12))))
    # graphics.off()
    
    # plotObj <- obj[genes,grep("baseMean",colnames(obj))]
    # colnames(plotObj) <- gsub("_","",gsub(sample,"",gsub("baseMean","",colnames(plotObj))))
    # rownames(plotObj) <- sapply(rownames(plotObj), function(gene) strsplit(gene,"_")[[1]][2])
    pdf(file.path("figure6",paste("Heatmap_sample_",sample,".pdf",sep="")), height = plotHeight, width = 9)
    print(Heatmap(as.matrix(ranks), show_row_names=T, show_column_names=F, left_annotation=row_ha))
    # print(Heatmap(plotObj, name = "BaseMeans", column_title = "Condition",
    # height = unit(150*nrow(gm)/31, "mm"),show_column_names=T, show_row_names=T,
    # left_annotation=row_ha,
    # row_names_max_width = max_text_width(rownames(plotObj), gp = gpar(fontsize = 12)),
    # column_names_max_height = max_text_width(colnames(plotObj), gp = gpar(fontsize = 12))))
    graphics.off()
  }
}

# Enrichment analysis ####
mutations <- c("6:33320019_Goe615","MT:15452_Goe1360")
for(mut in mutations){
  name <- mut
  sample <- strsplit(mut, split="_")[[1]][2]
  obj <- markers[[mut]]
  genes_ranked <- data.frame(ID=obj$ENSEMBL, rank=obj[[grep("rank", colnames(obj))]])

  de_genes <- obj[obj[[grep("padj", colnames(obj))]]<=0.05,]
  de_up <- de_genes$ENSEMBL[de_genes[[grep("log2FoldChange", colnames(de_genes))]]>0]
  de_down <- de_genes$ENSEMBL[de_genes[[grep("log2FoldChange", colnames(de_genes))]]<0]
  universe <- obj$ENSEMBL
  
  terms <- c("geneontology_Biological_Process_noRedundant","geneontology_Cellular_Component_noRedundant","geneontology_Molecular_Function_noRedundant","pathway_KEGG","pathway_Reactome")
  
  dir.create(file.path("figure6",paste("Webgestalt_upReg_sample_",sample,sep="")))
  WebGestaltR(
    enrichMethod = "ORA",
    organism = "hsapiens",
    enrichDatabase = terms,
    interestGene=de_up,interestGeneType="ensembl_gene_id",
    referenceGene = universe, referenceGeneType="ensembl_gene_id",
    sigMethod="fdr", fdrThr=FDR,minNum=5,outputDirectory=file.path("figure6",paste("Webgestalt_upReg_sample_",sample,sep=""))
  )
  dir.create(file.path("figure6",paste("Webgestalt_downReg_sample_",sample,sep="")))
  WebGestaltR(
    enrichMethod = "ORA",
    organism = "hsapiens",
    enrichDatabase = terms,
    interestGene=de_up,interestGeneType="ensembl_gene_id",
    referenceGene = universe, referenceGeneType="ensembl_gene_id",
    sigMethod="fdr", fdrThr=FDR,minNum=5,outputDirectory=file.path("figure6",paste("Webgestalt_downReg_sample_",sample,sep=""))
  )
  dir.create(file.path("figure6",paste("Webgestalt_GSEA_sample_",sample,sep="")))
  WebGestaltR(
    enrichMethod = "GSEA",
    organism = "hsapiens",
    enrichDatabase = terms,
    interestGene=genes_ranked,
    interestGeneType="ensembl_gene_id",
    sigMethod="fdr", fdrThr=FDR,minNum=5,outputDirectory=file.path("figure6",paste("Webgestalt_GSEA_sample_",sample,sep=""))
  )
}

# GO
for(mut in mutations){
  name <- mut
  obj <- markers[[mut]]
  genes_ranked <- obj[[grep("rank", colnames(obj))]]
  names(genes_ranked) <- obj$ENSEMBL
  genes_ranked <- genes_ranked[!is.na(names(genes_ranked))]
  genes_ranked <- sort(genes_ranked, decreasing = T)

  de_genes <- obj[obj[[grep("padj", colnames(obj))]]<=0.05,]
  de_up <- de_genes$ENSEMBL[de_genes[[grep("log2FoldChange", colnames(de_genes))]]>0]
  de_down <- de_genes$ENSEMBL[de_genes[[grep("log2FoldChange", colnames(de_genes))]]<0]
  universe <- obj$ENSEMBL

  # Plot Gene-concept network
  gsea <- gseGO(geneList      = genes_ranked,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENSEMBL",
                ont           = "ALL",
                minGSSize     = 10,
                maxGSSize     = 500,
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                verbose       = FALSE)
  gsea <- setReadable(gsea, 'org.Hs.eg.db', 'ENSEMBL')

  ora_up <- enrichGO(gene        = de_up,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENSEMBL",
                ont           = "ALL",
                minGSSize     = 10,
                maxGSSize     = 500,
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                universe      = TRUE)
  ora_up <- setReadable(ora_up, 'org.Hs.eg.db', 'ENSEMBL')
  ora_down <- enrichGO(gene        = de_down,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENSEMBL",
                ont           = "ALL",
                minGSSize     = 10,
                maxGSSize     = 500,
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",
                universe      = TRUE)
  ora_down <- setReadable(ora_down, 'org.Hs.eg.db', 'ENSEMBL')

  if(nrow(gsea)>0){
    p1 <- cnetplot(gsea, foldChange=genes_ranked, cex_label_category = 2,
                   cex_label_gene = 1.5)
    pdf(file.path("figure6",paste("GO_GSEA_mut_",name,"_geneConceptMap.pdf", sep="")), width=40, height=20)
    print(p1)
    graphics.off()
    p2 <- heatplot(gsea, foldChange=genes_ranked)
    pdf(file.path("figure6",paste("GO_GSEA_mut_",name,"_heatmap.pdf", sep="")), width=55, height=10)
    print(p2)
    graphics.off()
  }
  if(nrow(gsea@result)>0) write.table(gsea@result, file.path("figure6",paste("GO_GSEA_mut_",name,".csv", sep="")), sep=",",row.names=F)

  if(nrow(ora_up)>0){
    p1 <- cnetplot(ora_up, cex_label_category = 2,
                   cex_label_gene = 1.5)
    pdf(file.path("figure6",paste("GO_ORA_upReg_mut_",name,"_geneConceptMap.pdf", sep="")), width=40, height=20)
    print(p1)
    graphics.off()
    p2 <- heatplot(ora_up)
    pdf(file.path("figure6",paste("GO_ORA_upReg_mut_",name,"_heatmap.pdf", sep="")), width=55, height=10)
    print(p2)
    graphics.off()
  }
  if(nrow(ora_up@result)>0) write.table(ora_up@result, file.path("figure6",paste("GO_ORA_upReg_mut_",name,".csv", sep="")), sep=",",row.names=F)

  if(nrow(ora_down)>0){
    p1 <- cnetplot(ora_down, cex_label_category = 2,
                   cex_label_gene = 1.5)
    pdf(file.path("figure6",paste("GO_ORA_downReg_mut_",name,"_geneConceptMap.pdf", sep="")), width=40, height=20)
    print(p1)
    graphics.off()
    p2 <- heatplot(ora_down)
    pdf(file.path("figure6",paste("GO_ORA_downReg_mut_",name,"_heatmap.pdf", sep="")), width=55, height=10)
    print(p2)
    graphics.off()
  }
  if(nrow(ora_down@result)>0) write.table(ora_down@result, file.path("figure6",paste("GO_ORA_downReg_mut_",name,".csv", sep="")), sep=",",row.names=F)
}

# KEGG
for(mut in mutations){
  name <- mut
  obj <- markers[[mut]]
  # genes_ranked <- obj[[grep("rank", colnames(obj))]]
  # names(genes_ranked) <- unlist(as.list(org.Hs.egENSEMBL2EG)[obj$ENSEMBL])
  # genes_ranked <- genes_ranked[!is.na(names(genes_ranked))]
  # genes_ranked <- sort(genes_ranked, decreasing = T)

  de_genes <- obj[obj[[grep("padj", colnames(obj))]]<=0.05,]
  de_up <- unlist(as.list(org.Hs.egENSEMBL2EG)[de_genes$ENSEMBL[de_genes[[grep("log2FoldChange", colnames(de_genes))]]>0]])
  de_down <- unlist(as.list(org.Hs.egENSEMBL2EG)[de_genes$ENSEMBL[de_genes[[grep("log2FoldChange", colnames(de_genes))]]<0]])
  universe <- unlist(as.list(org.Hs.egENSEMBL2EG)[obj$ENSEMBL])

  # Plot Gene-concept network
  # gsea <- gseKEGG(geneList     = genes_ranked,
  #               organism = "hsa",
  #               keyType      = "kegg",
  #               minGSSize    = 10,
  #               maxGSSize    = 500,
  #               pvalueCutoff = 0.05,
  #               pAdjustMethod = "BH",
  #               verbose      = FALSE)
  # gsea <- setReadable(gsea, 'org.Hs.eg.db', 'ENTREZID')

  ora_up <- enrichKEGG(gene     = de_up,
                       organism = "hsa",
                       keyType      = "kegg",
                       minGSSize    = 10,
                       maxGSSize    = 500,
                       pvalueCutoff = 0.05,
                       universe = universe,
                       pAdjustMethod = "BH")
  ora_down <- enrichKEGG(gene     = de_down,
                         organism = "hsa",
                         keyType      = "kegg",
                         minGSSize    = 10,
                         maxGSSize    = 500,
                         pvalueCutoff = 0.05,
                         universe = universe,
                         pAdjustMethod = "BH")

  # if(nrow(gsea)>0){
  #   p1 <- cnetplot(gsea, foldChange=genes_ranked, cex_label_category = 2,
  #                  cex_label_gene = 1.5)
  #   pdf(file.path("figure6",paste("KEGG_GSEA_mut_",name,"_geneConceptMap.pdf", sep="")), width=40, height=20)
  #   print(p1)
  #   graphics.off()
  #   p2 <- heatplot(gsea, foldChange=genes_ranked)
  #   pdf(file.path("figure6",paste("KEGG_GSEA_mut_",name,"_heatmap.pdf", sep="")), width=55, height=10)
  #   print(p2)
  #   graphics.off()
  # }
  # if(nrow(gsea@result)>0) write.table(gsea@result, file.path("figure6",paste("KEGG_GSEA_mut_",name,".csv", sep="")), sep=",", row.names=F)

  if(!is.null(ora_up) && nrow(ora_up)>0){
    ora_up <- setReadable(ora_up, 'org.Hs.eg.db', 'ENTREZID')
    p1 <- cnetplot(ora_up, cex_label_category = 2, cex_label_gene = 1.5)
    pdf(file.path("figure6",paste("KEGG_ORA_upReg_mut_",name,"_geneConceptMap.pdf", sep="")), width=40, height=20)
    print(p1)
    graphics.off()
    p2 <- heatplot(ora_up)
    pdf(file.path("figure6",paste("KEGG_ORA_upReg_mut_",name,"_heatmap.pdf", sep="")), width=55, height=10)
    print(p2)
    graphics.off()
  }
  if(!is.null(ora_up) && nrow(ora_up@result)>0) write.table(ora_up@result, file.path("figure6",paste("KEGG_ORA_upReg_mut_",name,".csv", sep="")), sep=",", row.names=F)

  if(!is.null(ora_down) && nrow(ora_down)>0){
    ora_down <- setReadable(ora_down, 'org.Hs.eg.db', 'ENTREZID')
    p1 <- cnetplot(ora_down, cex_label_category = 2, cex_label_gene = 1.5)
    pdf(file.path("figure6",paste("KEGG_ORA_downReg_mut_",name,"_geneConceptMap.pdf", sep="")), width=40, height=20)
    print(p1)
    graphics.off()
    p2 <- heatplot(ora_down)
    pdf(file.path("figure6",paste("KEGG_ORA_downReg_mut_",name,"_heatmap.pdf", sep="")), width=55, height=10)
    print(p2)
    graphics.off()
  }
  if(!is.null(ora_down) && nrow(ora_down@result)>0) write.table(ora_down@result, file.path("figure6",paste("KEGG_ORA_downReg_mut_",name,".csv", sep="")), sep=",", row.names=F)
}

# REACTOME
for(mut in mutations){
  name <- mut
  obj <- markers[[mut]]
  # genes_ranked <- obj[[grep("rank", colnames(obj))]]
  # names(genes_ranked) <- unlist(as.list(org.Hs.egENSEMBL2EG)[obj$ENSEMBL])
  # genes_ranked <- genes_ranked[!is.na(names(genes_ranked))]
  # genes_ranked <- sort(genes_ranked, decreasing = T)

  de_genes <- obj[obj[[grep("padj", colnames(obj))]]<=0.05,]
  de_up <- unlist(as.list(org.Hs.egENSEMBL2EG)[de_genes$ENSEMBL[de_genes[[grep("log2FoldChange", colnames(de_genes))]]>0]])
  de_down <- unlist(as.list(org.Hs.egENSEMBL2EG)[de_genes$ENSEMBL[de_genes[[grep("log2FoldChange", colnames(de_genes))]]<0]])
  universe <- unlist(as.list(org.Hs.egENSEMBL2EG)[obj$ENSEMBL])

  # # Plot Gene-concept network
  # gsea <- gsePathway(geneList     = genes_ranked,
  #                 organism = "human",
  #                 minGSSize    = 10,
  #                 maxGSSize    = 500,
  #                 pvalueCutoff = 0.05,
  #                 pAdjustMethod = "BH",
  #                 verbose      = FALSE)
  # gsea <- setReadable(gsea, 'org.Hs.eg.db', 'ENTREZID')

  ora_up <- enrichPathway(gene     = de_up,
                    organism = "human",
                    minGSSize    = 10,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.05,
                    universe = universe,
                    pAdjustMethod = "BH")
  ora_down <- enrichPathway(gene     = de_down,
                    organism = "human",
                    minGSSize    = 10,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.05,
                    universe = universe,
                    pAdjustMethod = "BH")

  # if(nrow(gsea)>0){
  #   ora_up <- setReadable(ora_up, 'org.Hs.eg.db', 'ENTREZID')
  #   p1 <- cnetplot(gsea, foldChange=genes_ranked, cex_label_category = 2,
  #                  cex_label_gene = 1.5)
  #   pdf(file.path("figure6",paste("Reactome_GSEA_mut_",name,"_geneConceptMap.pdf", sep="")), width=40, height=20)
  #   print(p1)
  #   graphics.off()
  #   p2 <- heatplot(gsea, foldChange=genes_ranked)
  #   pdf(file.path("figure6",paste("Reactome_GSEA_mut_",name,"_heatmap.pdf", sep="")), width=55, height=10)
  #   print(p2)
  #   graphics.off()
  # }
  # if(nrow(gsea@result)>0) write.table(gsea@result, file.path("figure6",paste("Reactome_GSEA_mut_",name,".csv", sep="")), sep=",", row.names=F)

  if(!is.null(ora_up) && nrow(ora_up)>0){
    ora_up <- setReadable(ora_up, 'org.Hs.eg.db', 'ENTREZID')
    p1 <- cnetplot(ora_up, cex_label_category = 2, cex_label_gene = 1.5)
    pdf(file.path("figure6",paste("Reactome_ORA_upReg_mut_",name,"_geneConceptMap.pdf", sep="")), width=40, height=20)
    print(p1)
    graphics.off()
    p2 <- heatplot(ora_up)
    pdf(file.path("figure6",paste("Reactome_ORA_upReg_mut_",name,"_heatmap.pdf", sep="")), width=55, height=10)
    print(p2)
    graphics.off()
  }
  if(!is.null(ora_up) && nrow(ora_up@result)>0) write.table(ora_up@result, file.path("figure6",paste("Reactome_ORA_upReg_mut_",name,".csv", sep="")), sep=",", row.names=F)
  
  if(!is.null(ora_down) && nrow(ora_down)>0){
    ora_down <- setReadable(ora_down, 'org.Hs.eg.db', 'ENTREZID')
    p1 <- cnetplot(ora_down, cex_label_category = 2, cex_label_gene = 1.5)
    pdf(file.path("figure6",paste("Reactome_ORA_downReg_mut_",name,"_geneConceptMap.pdf", sep="")), width=40, height=20)
    print(p1)
    graphics.off()
    p2 <- heatplot(ora_down)
    pdf(file.path("figure6",paste("Reactome_ORA_downReg_mut_",name,"_heatmap.pdf", sep="")), width=55, height=10)
    print(p2)
    graphics.off()
  }
  if(!is.null(ora_down) && nrow(ora_down@result)>0) write.table(ora_down@result, file.path("figure6",paste("Reactome_ORA_downReg_mut_",name,".csv", sep="")), sep=",", row.names=F)
}

save(CogentDS_data,file = "figure6/CogentDS_samples_icell8_abslowcov5000_withMuts.rda")
