require(CogentDS)

# ICELL8 dataset

# Default analysis for samples
standard_analysis(gm = "analysis_L4_old/analysis_L4_genematrix.csv", metadata = "analysis_L4_old/analysis_L4_new_stats.csv",
                  gene_info_loc=geneInfo, cluster_analysis=T, grouping_var="Sample",
                  parallel_cores=100, report=T, CogentAP_data = T, output_dir = file.path(OutDir,"CogentDS"))
