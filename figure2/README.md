Figure 2: QC metrics of both approaches (2C) and number of genes expressed (2D)

To generate the plot for figure 2, run the commands in the following order:

1) Download single-cell and bulk data from NCBI
Rscript Download_NCBI.R

2) Generate the CogentDS reports for the single-cell datasets
Rscript generate_CogentDS_reports.R

3) Generate figure 2C (bar plots for scRNA-seq metrics from CogentDS report)
Rscript figure2c.R

4) Generate figure 2D (Venn and violin plots for number of expressed genes in cells in each approach)
Rscript figure2d.R