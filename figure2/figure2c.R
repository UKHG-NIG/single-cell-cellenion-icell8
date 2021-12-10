# This script generates bar plots of the QC metrics from the two scRNA-seq approaches (Fig2c)

# Load supporting libraries
require(ggpubr)
require(reshape2)

# Load datasets
## Composite approach
composite <- get(load("CogentDS_CellenONE/CogentDS.analysis.rda"))
## ICELL8 approach
icell8 <- get(load("CogentDS_ICELL8/CogentDS.analysis.rda"))

# Extract data for columns of interest and sum them up for all cells
columns <- c("Barcoded_Reads","Trimmed_Reads", "Unmapped_Reads", "Mapped_Reads",
             "Uniquely_Mapped_Reads","Multimapped_Reads", "Exon_Reads",
             "Ambiguous_Exon_Reads", "Intron_Reads", "Ambiguous_Intron_Reads",
             "Intergenic_Reads", "Mitochondrial_Reads", "Ribosomal_Reads")
df <- t(data.frame(composite=colSums(composite$raw_data$metadata[,columns]),
                 icell8=colSums(icell8$raw_data$metadata[,columns])))
# Calculate percentages for each QC metric
df <- df/df[,1]*100
# Melt data.frame to allow easier bar plotting
df <- melt(df)
colnames(df) <- c("approach","metrics", "value")
df$approach <- as.character(df$approach)
df$metrics <- as.character(df$metrics)
# Plot bars
ggbarplot(df, x = "metrics", y = "value",# ylim=c(0,100),
          fill = "approach", ylab="%Reads",
          color = "approach", palette = c("#999999", "#E69F00"),
          label = F, position = position_dodge(0.9),
          x.text.angle = 45           # Rotate vertically x axis texts
)+scale_y_continuous(breaks = get_breaks(from=0,to=100, by=10))
