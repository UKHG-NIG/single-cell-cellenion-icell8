require(CogentDS)

# ICELL8 dataset
## Create text file referring to gene matrix files
fileConn<-file("GSM5411464_geneMatrix_files.txt")
writeLines(c("GSM5411464_ICELL8_geneMatrix_chip1.csv","GSM5411465_ICELL8_geneMatrix_chip2.csv","GSM5411466_ICELL8_geneMatrix_chip3.csv"), fileConn)
close(fileConn)
## Create text file referring to metadata files
fileConn<-file("GSM5411464_metadata_files.txt")
writeLines(c("GSM5411464_ICELL8_metadata_chip1.csv","GSM5411465_ICELL8_metadata_chip2.csv","GSM5411466_ICELL8_metadata_chip3.csv"), fileConn)
close(fileConn)

standard_analysis(gm_loc = "GSM5411464_geneMatrix_files.txt", metadata_loc = "GSM5411464_metadata_files.txt", gene_info_loc="GSE179204_gene_info.csv",
                  cluster_analysis=T, grouping_var="Sample", names=c("chip1","chip2","chip3"), report=T, CogentAP_data = T, output_dir = "CogentDS_ICELL8",
                  cell_type=c("GOE1303","GOE1305","GOE1309","GOE1360","GOE247","GOE486","GOE615","GOE800"), mito_frac=c(0,1), intergenic_frac=c(0,1))

# CellenONE/ICELL8 dataset
## Default analysis for samples
standard_analysis(gm = "GSM5411467_CellenONE_geneMatrix.csv", metadata = "GSM5411467_CellenONE_metadata.csv", gene_info_loc="GSM5411467_gene_info.csv",
                  cluster_analysis=T, grouping_var="Sample", report=T, CogentAP_data = T, output_dir = "CogentDS_CellenONE",
                  cell_type=c("Goe1303","Goe1305","Goe1309","Goe1360","Goe247","Goe486","Goe615","Goe800"), mito_frac=c(0,1), intergenic_frac=c(0,1))
