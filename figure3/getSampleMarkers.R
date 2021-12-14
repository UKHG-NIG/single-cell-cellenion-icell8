# This script finds the markers in each of the approaches, and prints them out into CSV files
# Set wokring directory to input directory where CogentDS directoeries are
setwd("/path/to/InputDir")

# Load libraries
library(CogentDS)

# Load RDA objects from CogentDS analyses
composite <- get(load("CogentDS_CellenONE/CogentDS.analysis.rda"))
icell8 <- get(load("CogentDS_ICELL8/CogentDS.analysis.rda"))

# Get markers for ICELL8 data
icell8_samples <- grep(pattern="_Ctrl", icell8$qc_data$metadata$Sample, invert=T)
icell8_gm <- icell8$qc_data$gm[,icell8_samples]
icell8_metadata <- icell8$qc_data$metadata[icell8_samples,]
icell8_markers <- de.markers(gm = icell8_gm,
                             metadata = icell8_metadata,
                             min_fc = 0,
                             grouping_var = "Sample",
                             transform_type = icell8$qc_data$params$gm_log_base)
icell8_markers <- lapply(icell8_markers, function(obj){
  obj$rank <- -log10(obj$pval)*obj$log_fc
  obj[order(obj$rank, decreasing = T),]
})
names(icell8_markers) <- toupper(names(icell8_markers))

# Get markers for composite approach data
composite_samples <- grep(pattern="Goe247_wo__fix", composite$qc_data$metadata$Sample, invert=T)
composite_gm <- composite$qc_data$gm[,composite_samples]
composite_metadata <- composite$qc_data$metadata[composite_samples,]
composite_markers <- de.markers(gm = composite_gm,
                             metadata = composite_metadata,
                             min_fc = 0,
                             grouping_var = "Sample",
                             transform_type = composite$qc_data$params$gm_log_base)
composite_markers <- lapply(composite_markers, function(obj){
  obj$rank <- -log10(obj$pval)*obj$log_fc
  obj[order(obj$rank, decreasing = T),]
})
names(composite_markers) <- toupper(names(composite_markers))

# Print out markers for each sample
for(sample in names(icell8_markers)){
  # ICELL8 markers
  write.table(icell8_markers[[sample]],
              file = paste("markers_icell8_",sample,".csv", sep=""),
              sep=",", quote = F, row.names=T, col.names = NA)
  
  # Composite markers
  write.table(composite_markers[[sample]],
              file = paste("markers_composite_",sample,".csv", sep=""),
              sep=",", quote = F, row.names=T, col.names = NA)
  
  # Merge marker tables together
  markersICELL8 <- icell8_markers[[sample]]
  markersComposite <- composite_markers[[sample]]
  colnames(markersICELL8) <- paste("icell8",colnames(markersICELL8),sep="_")
  colnames(markersComposite) <- paste("cellenion",colnames(markersComposite),sep="_")
  df <- merge(markersICELL8, markersComposite, by="row.names",all=T)
  df <- df[complete.cases(df), ]
  # Sort markers by sum of absolute rank scores in both approaches
  df <- df[order(abs(df$icell8_rank)+abs(df$cellenion_rank), decreasing=TRUE),]
  # Select top 100 markers to plot
  df <- df[1:100,]
  write.table(df[,-1],sep=",",
              file = paste("top100_markers_",sample,".csv", sep=""),
              quote = F, row.names=df$Row.names, col.names = NA)
}
