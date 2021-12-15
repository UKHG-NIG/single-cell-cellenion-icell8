Workflow to generate figure 3 and Supplementary Figures 3 and 4 are as follows:

1) Obtain the markers from each sample in each approach, and write them out
Rscript getSampleMarkers.R

2) Generate the correlation plots for the top 100 ranked markers in both ICELL8 and composite approach (figure 3A)
Rscript figure3a.R

3) Get the gene counts from the bulk RNA-seq data
Rscript Download_NCBI_bulk.R

4) Generate the correlation plots for the top 100 ranked markers from (2) for pseudocounts of the scRNA-seq approaches and compare them to the bulk data (figure 3b and Supplementary Figure 3)
Rscript figure3b.R

5) Generate the Venn diagrams for the markers found in both scRNA-seq approaches (figure 3c and Supplementary Figure 4)
Rscript figure3c.R
