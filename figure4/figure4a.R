library(CogentDS)
library(umap)

# Figure 4A: compare tSNE of cellenion with icell8 ####
load("figure4/CogentDS_samples_cellenion_abslowcov5000/CogentDS.analysis.rda")
Cogent_cellenion <- CogentDS_data
load("figure4/CogentDS_samples_icell8_abslowcov5000/CogentDS.analysis.rda")
Cogent_icell8 <- CogentDS_data

# CellenONE
gm = t(Cogent_cellenion[["pca_data"]]$pca_obj$x[, 1:Cogent_cellenion[["pca_data"]]$pc_count])
metadata = Cogent_cellenion[["qc_data"]]$metadata
tsne_obj <- gm.tsne(gm = gm, metadata = metadata,
                    grouping_var = grouping_var, dims = 2, theta = 0.25,
                    max_iter = 1500, pca = FALSE, pca_center = FALSE, 
                    pca_scale = FALSE, perplexity = 20, plot = F)
reduction.plot(vis_obj = tsne_obj, metadata = metadata, 
               grouping_var = "Sample", marker_size=0.5, alpha=1, file="figure4/Fig4A_CellenONE_tSNE.pdf")

# icell8
gm = t(Cogent_icell8[["pca_data"]]$pca_obj$x[, 1:Cogent_icell8[["pca_data"]]$pc_count])
metadata = Cogent_icell8[["qc_data"]]$metadata
tsne_obj <- gm.tsne(gm = gm, metadata = metadata,
                    grouping_var = grouping_var, dims = 2, theta = 0.25,
                    max_iter = 1500, pca = FALSE, pca_center = FALSE, 
                    pca_scale = FALSE, perplexity = 20, plot = F)
reduction.plot(vis_obj = tsne_obj, metadata = metadata, 
               grouping_var = "Sample", marker_size=0.5, alpha=1,file="figure4/Fig4A_ICELL8_tSNE.pdf")
