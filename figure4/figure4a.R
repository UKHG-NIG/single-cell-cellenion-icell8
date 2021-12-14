library(CogentDS)
library(umap)

setwd("/zfstank/ngsdata/projects/icell8_cellenion_publication")

# # Get UMAP for iCell8 samples
# gm_loc<-"/zfstank/ngsdata/Hiseq/210216_K00295_0247_AHKVTJBBXY/Aligned/Wollnik_Cogent/merge_chips/merged_gm"
# metadata_loc<-"/zfstank/ngsdata/Hiseq/210216_K00295_0247_AHKVTJBBXY/Aligned/Wollnik_Cogent/merge_chips/merged_stats"
# gene_info_loc<-"/zfstank/ngsdata/Hiseq/210316_K00295_0279_BHKVJNBBXY/Aligned/Wollnik/analysis_L4/gene_info.csv"
# icell8_samples <- c("GOE247","GOE486","GOE615","GOE800","GOE1303","GOE1305","GOE1309","GOE1360")
# standard_analysis(gm_loc = gm_loc, metadata_loc = metadata_loc,
#                   gene_info_loc=gene_info_loc, cluster_analysis=T, grouping_var="Sample", cell_type=icell8_samples,
#                   mito_frac = list(0,1), intergenic_frac = list(0,1), auto_outlier_removal=F,
#                   qc_cell_abslowcov = 5000,
#                   qc_cells_abslowgenecount = 200,
#                   parallel_cores=100, report=T, CogentAP_data = T, output_dir = file.path("figure4","CogentDS_samples_icell8_abslowcov5000"))
# 
# standard_analysis(gm_loc = gm_loc, metadata_loc = metadata_loc,names=c("chip1","chip2","chip3"),
#                   gene_info_loc=gene_info_loc, cluster_analysis=T, grouping_var="Group", cell_type=icell8_samples,
#                   mito_frac = list(0,1), intergenic_frac = list(0,1), auto_outlier_removal=F,
#                   qc_cell_abslowcov = 5000,
#                   qc_cells_abslowgenecount = 200,
#                   parallel_cores=100, report=T, CogentAP_data = T, output_dir = file.path("figure4","CogentDS_chips_icell8_abslowcov5000"))
# 
# # Get UMAP for cellenion samples
# gm<-"/zfstank/ngsdata/Hiseq/210316_K00295_0279_BHKVJNBBXY/Aligned/Wollnik/merge_lanes/analysis_merged_genematrix.csv"
# metadata<-"/zfstank/ngsdata/Hiseq/210316_K00295_0279_BHKVJNBBXY/Aligned/Wollnik/merge_lanes/analysis_merged_stats.csv"
# gene_info_loc<-"/zfstank/ngsdata/Hiseq/210316_K00295_0279_BHKVJNBBXY/Aligned/Wollnik/analysis_L4/gene_info.csv"
# cellenion_samples <- c("Goe247", "Goe486", "Goe615", "Goe800", "Goe1303", "Goe1305", "Goe1309","Goe1360")
# standard_analysis(gm = gm, metadata = metadata,
#                   gene_info_loc=gene_info_loc, cluster_analysis=T, grouping_var="Sample", cell_type=cellenion_samples,
#                   mito_frac = list(0,1), intergenic_frac = list(0,1), auto_outlier_removal=F,
#                   qc_cell_abslowcov = 5000,
#                   qc_cells_abslowgenecount = 200,
#                   parallel_cores=100, report=T, CogentAP_data = T, output_dir = file.path("figure4","CogentDS_samples_cellenion_abslowcov5000"))
# standard_analysis(gm = gm, metadata = metadata,
#                   gene_info_loc=gene_info_loc, cluster_analysis=T, grouping_var="Sample", cell_type=cellenion_samples,
#                   mito_frac = list(0,1), intergenic_frac = list(0,1), auto_outlier_removal=T,
#                   qc_cell_abslowcov = 5000,
#                   qc_cells_abslowgenecount = 200,
#                   parallel_cores=100, report=T, CogentAP_data = T, output_dir = file.path("figure4","CogentDS_samples_cellenion_abslowcov5000_noOutlierRemoval"))

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

# gm = Cogent_cellenion[["pca_data"]]$pca_obj$x[,1:Cogent_cellenion[["pca_data"]]$pc_count]
# umap_obj <- gm.umap(gm = gm, metadata = metadata, grouping_var = grouping_var, plot = F)
# reduction.plot(vis_obj = umap_obj, metadata = metadata, 
#                grouping_var = grouping_var, file = "figure4/Fig4A_CellenONE_UMAP.png", marker_size=0.5)

# icell8
gm = t(Cogent_icell8[["pca_data"]]$pca_obj$x[, 1:Cogent_icell8[["pca_data"]]$pc_count])
metadata = Cogent_icell8[["qc_data"]]$metadata
tsne_obj <- gm.tsne(gm = gm, metadata = metadata,
                    grouping_var = grouping_var, dims = 2, theta = 0.25,
                    max_iter = 1500, pca = FALSE, pca_center = FALSE, 
                    pca_scale = FALSE, perplexity = 20, plot = F)
reduction.plot(vis_obj = tsne_obj, metadata = metadata, 
               grouping_var = "Sample", marker_size=0.5, alpha=1,file="figure4/Fig4A_ICELL8_tSNE.pdf")

# gm = Cogent_icell8[["pca_data"]]$pca_obj$x[,1:Cogent_icell8[["pca_data"]]$pc_count]
# umap_obj <- gm.umap(gm = gm, metadata = metadata, grouping_var = grouping_var, plot = F)
# reduction.plot(vis_obj = umap_obj, metadata = metadata, 
#                grouping_var = grouping_var, file = "figure4/Fig4A_ICELL8_UMAP.png", marker_size=0.5)
