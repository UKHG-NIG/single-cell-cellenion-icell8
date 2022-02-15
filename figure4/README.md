# Analysis of transcriptional heterogeneity of dermal fibroblast derived from patients using the novel approach (Figure 4, Supplementary Figure 4)

1. Generate the CogentDS reports for the ICELL8 and novel approach datasets (Figure 4a)
```
Rscript figure4a.R
```

2. Plot a Circo plot displaying the top 100 markers found for each sample, and the extent to which they overlap between the different samples (Figure 4B)
```
Rscript figure4b.R
```

3. Using the CogentDS objects, plot the following plots:
    - t-SNE with unsupervised cluster (Figure 4C)
    - the gene expression gradients for genes MMP2, H2AX and NDUFA12 in all cells (Figure 4D)
    - heatmap for the top 30 ranked markers in each cluster in the dataset from the novel approach (Figure 4E)
    - heatmap for the top 30 ranked markers in each cluster in the dataset from the ICELL8 (Supp. Fig. 4A)
    - a venn diagram overlapping the markers from the two heatmaps above (Supp. Fig. 4B)
```
Rscript figure4cde.R
```
