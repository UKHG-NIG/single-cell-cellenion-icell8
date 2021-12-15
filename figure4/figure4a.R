require(CogentDS)

# ICELL8 dataset
standard_analysis(gm_loc = "GSM5411464_geneMatrix_files.txt", metadata_loc = "GSM5411464_metadata_files.txt", gene_info_loc="GSE179204_gene_info.csv",
                  cluster_analysis=T, grouping_var="Sample", names=c("chip1","chip2","chip3"), report=T, CogentAP_data = T, output_dir = "CogentDS_ICELL8_fig4",
                  qc_cell_abslowcov=5000, qc_cells_abslowgenecount = 200, mito_frac=c(0,1), intergenic_frac=c(0,1), auto_outlier_removal=F,
                  cell_type=c("GOE247","GOE486","GOE615","GOE800","GOE1303","GOE1305","GOE1309","GOE1360"), reduction_type="tSNE")

# CellenONE/ICELL8 dataset
standard_analysis(gm = "GSM5411467_CellenONE_geneMatrix.csv", metadata = "GSM5411467_CellenONE_metadata.csv", gene_info_loc="GSM5411467_gene_info.csv",
                  cluster_analysis=T, grouping_var="Sample", report=T, CogentAP_data = T, output_dir = "CogentDS_CellenONE_fig4",
                  qc_cell_abslowcov=5000, qc_cells_abslowgenecount = 200 mito_frac=c(0,1), intergenic_frac=c(0,1), auto_outlier_removal=F,
                  cell_type=c("Goe247", "Goe486", "Goe615", "Goe800", "Goe1303", "Goe1305", "Goe1309","Goe1360"), reduction_type="tSNE")
