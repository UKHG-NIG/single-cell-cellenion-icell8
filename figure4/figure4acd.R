# This scripts is used to generate tSNE plots for composite approach and ICELL8 platform
library(CogentDS)

# Figure 4A: compare tSNE of composite approach with icell8 platform ####
composite <- get(load("CogentDS_CellenONE/CogentDS.analysis.rda"))
icell8 <- get(load("CogentDS_ICELL8/CogentDS.analysis.rda"))

# Composite
composite_metadata <- composite[["qc_data"]]$metadata
composite_metadata <- composite_metadata[composite_metadata$Sample!="Goe247_wo__fix",]
composite_pca <- gm.pca(gm_qc = as.matrix(composite[["qc_data"]]$gm[,rownames(composite_metadata)]),
                        gm_raw = as.matrix(composite[["raw_data"]]$gm[,rownames(composite_metadata)]), metadata = composite_metadata,
                        grouping_var = "Sample", filt_method = "top_var",
                        thresh_cut = 0, top_genes = 2000,
                        transform_type = "none",
                        plot = FALSE, var_weight = FALSE)
composite_tsne <- gm.tsne(gm = t(composite_pca$pca_obj$x[,1:composite_pca$pc_count]), metadata = composite_metadata,
                    grouping_var = "Sample", dims = 2, theta = 0.25,
                    max_iter = 1500, pca = FALSE, pca_center = FALSE, 
                    pca_scale = FALSE, perplexity = 20, plot = F)
reduction.plot(vis_obj = composite_tsne, metadata = composite_metadata, 
               grouping_var = "Sample", marker_size=0.5, alpha=1, file="Fig4A_composite_samples_tSNE.pdf")

# icell8
icell8_metadata <- icell8[["qc_data"]]$metadata
icell8_metadata <- icell8_metadata[!grepl("Ctrl",icell8_metadata$Sample),]
icell8_pca <- gm.pca(gm_qc = as.matrix(icell8[["qc_data"]]$gm[,rownames(icell8_metadata)]),
                        gm_raw = as.matrix(icell8[["raw_data"]]$gm[,rownames(icell8_metadata)]), metadata = icell8_metadata,
                        grouping_var = "Sample", filt_method = "top_var",
                        thresh_cut = 0, top_genes = 2000,
                        transform_type = "none",
                        plot = FALSE, var_weight = FALSE)
icell8_tsne <- gm.tsne(gm = t(icell8_pca$pca_obj$x[,1:icell8_pca$pc_count]), metadata = icell8_metadata,
                          grouping_var = "Sample", dims = 2, theta = 0.25,
                          max_iter = 1500, pca = FALSE, pca_center = FALSE, 
                          pca_scale = FALSE, perplexity = 20, plot = F)
reduction.plot(vis_obj = icell8_tsne, metadata = icell8_metadata, 
               grouping_var = "Sample", marker_size=0.5, alpha=1,file="Fig4A_ICELL8_samples_tSNE.pdf")

# Figure 4C: tSNE plot of composite approach with graph-based clusters ####
# TODO: re-cluster cells with graph-based method and generate tSNE plot
CogentDS.cluster(m = CogentDS_data[["pca_data"]]$pca_obj$x[, 
            1:CogentDS_data[["pca_data"]]$pc_count], metadata = CogentDS_data[["qc_data"]]$metadata, 
            method = "Graph-based")
reduction.plot(vis_obj = composite_tsne, metadata = composite_metadata, 
               grouping_var = "Graph-based_Clusters", marker_size=0.5, alpha=1, file="Fig4C_composite_clusters_tSNE.pdf")

# Figure 4D: tSNE plots for specific markers in data ####
reduction.plot(vis_obj = composite_tsne, gm=composite[["qc_data"]]$gm,
               metadata = composite_metadata, 
               gene_highlight = "ENSG00000087245_MMP2", file="Fig4D_composite_MMP2_tSNE.pdf", marker_size=0.5, alpha=1)
reduction.plot(vis_obj = composite_tsne, gm=composite[["qc_data"]]$gm,
               metadata = composite_metadata, 
               gene_highlight = "ENSG00000188486_H2AX", file="Fig4D_composite_H2AX_tSNE.pdf", marker_size=0.5, alpha=1)
reduction.plot(vis_obj = composite_tsne, gm=composite[["qc_data"]]$gm,
               metadata = composite_metadata, 
               gene_highlight = "ENSG00000184752_NDUFA12", file="Fig4D_composite_NDUFA12_tSNE.pdf", marker_size=0.5, alpha=1)

# Figure 4E: heatmap of top 30 markers based on their ranks ####

# Supp Figure 4a: find top 30 markers for each ICELL8 sample, and generate a heatmap similar to 4E

# Supp Figure 4b: Venn diagram of top 30 CellenONE and ICELL8 markers
