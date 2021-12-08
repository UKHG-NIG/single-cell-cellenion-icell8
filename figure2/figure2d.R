# Generate Venn diagram and violin plots comparing number of genes detected in each approach
setwd("/path/to/InputDir")

# Load supporting libraries
require(VennDiagram)
require(reshape2)

# Load datasets
## Composite approach
composite <- get(load("CogentDS_CellenONE/CogentDS.analysis.rda"))
## ICELL8 approach
icell8 <- get(load("CogentDS_ICELL8/CogentDS.analysis.rda"))

# Plot Venn diagram for quality-control (QC) genes between the approaches
grid.draw(venn.diagram(x = list("CellenONE-ICELL8"=rownames(composite$qc_data$gm),
                                "ICELL8"=rownames(icell8$qc_data$gm)), col = "transparent",
                       fill = c("blue", "yellow"),# fill = c("blue", "green","red"),
                       alpha = 0.5, cat.cex=1.5, cex=1.5, height=100, width=100,filename=NULL,
                       print.mode=c("raw","percent")))

# Plot violin plot for genes detected per cell
# Raw reads counted for different platforms 
## Icell8 
icell8_mat <- icell8$qc_data$metadata[,c("Sample","No_of_Genes")]
icell8_mat$platform="ICELL8"
icell8_mat <- icell8_mat[!(icell8_mat$Sample %in% c("Pos_Ctrl","Neg_Ctrl")),]
## Cellenion 
cellenion_mat <- composite$qc_data$metadata[,c("Sample","No_of_Genes")]
cellenion_mat$platform="CellenONE"
cellenion_mat <- cellenion_mat[!(cellenion_mat$Sample %in% "Goe247_wo__fix"),]
# Plot 
plotDF <- melt(rbind(icell8_mat, cellenion_mat))
ggplot(plotDF[plotDF$variable=="No_of_Genes",], aes(factor(platform), value, fill=platform))+
  geom_violin(alpha=0.5) + scale_y_log10() +geom_boxplot(width = 0.2, fill="white")+theme_bw()+
  scale_fill_manual(values = c("blue", "yellow"))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  labs(title="Number of genes detected in cells (read count >0)", x="Approach", y="Genes detected per cell")

