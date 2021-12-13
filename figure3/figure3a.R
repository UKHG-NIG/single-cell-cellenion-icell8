# This script is used to correlate the ranking scores of the single-cell approaches (Figure 3A)
library(ggpubr)

# Generate correlation plot for top 100 markers in both 
samples <- c("GOE1303","GOE1305","GOE1309","GOE1360","GOE247","GOE486","GOE615","GOE800")
for(sample in samples){
  # Load markers for each sample in each dataset
  markersComposite <- read.table(paste("markers_composite_",sample,".csv",sep=""), sep = ",", header=T, row.names = 1)
  markersICELL8 <- read.table(paste("markers_icell8_",sample,".csv",sep=""), sep = ",", header=T, row.names = 1)
  colnames(markersICELL8) <- paste("icell8",colnames(markersICELL8),sep="_")
  colnames(markersComposite) <- paste("cellenion",colnames(markersComposite),sep="_")
  
  # Merge marker tables together
  df <- merge(markersICELL8, markersComposite, by="row.names",all=T)
  df <- df[complete.cases(df), ]
  # Sort markers by sum of absolute rank scores in both approaches
  df <- df[order(abs(df$icell8_rank)+abs(df$cellenion_rank), decreasing=TRUE),]
  # Select top 100 markers to plot
  df <- df[1:100,]
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
