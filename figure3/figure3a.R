# This script is used to correlate the ranking scores of the single-cell approaches (Figure 3A)
library(ggpubr)

# Generate correlation plot for top 100 markers in both 
samples <- c("GOE1303","GOE1305","GOE1309","GOE1360","GOE247","GOE486","GOE615","GOE800")
for(sample in samples){
  # Load markers for each sample in each dataset
  df <- read.table(paste("top100_markers_",sample,".csv",sep=""),sep=",", header=T, row.names=1)
  plotDF <- data.frame(icell8=log(df[[grep("icell8",grep(sample,colnames(df), ignore.case=T, value=T), value=T)]]),
                       cellenion=log(df[[grep("cellenion",grep(sample,colnames(df), ignore.case=T, value=T), value=T)]]))
  
  pdf(paste("Fig3a_Correlation_geneExpression_top100Icell8plusCellenion_",sample,".pdf",sep=""))
  print(ggscatter(plotDF, x = "icell8",xlab="",
                  title="",
                  y = "cellenion", ylab="",
                  add = "reg.line",
                  conf.int = TRUE)+
          geom_abline(intercept = 0, slope = 1)+
          stat_cor(label.y = 4.4) +
          stat_regline_equation(label.y = 4.2))
  graphics.off()
}
