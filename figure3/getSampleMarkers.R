# This script finds the markers in each of the approaches, and prints them out into CSV files
# Set wokring directory to input directory where CogentDS directoeries are
setwd("/path/to/InputDir")

# Load libraries
library(CogentDS)

# Load RDA objects from CogentDS analyses
Cogent_cellenion <- get(load("CogentDS_CellenONE/CogentDS.analysis.rda"))
Cogent_icell8 <- get(load("CogentDS_ICELL8/CogentDS.analysis.rda"))

# Get markers for icell8
icell8_samples <- grep(pattern="_Ctrl", Cogent_icell8$qc_data$metadata$Sample, invert=T)
icell8_gm <- Cogent_icell8$qc_data$gm[,icell8_samples]
icell8_metadata <- Cogent_icell8$qc_data$metadata[icell8_samples,]
icell8_markers <- de.markers(gm = icell8_gm,
                             metadata = icell8_metadata,
                             min_fc = 0,
                             grouping_var = "Sample",
                             transform_type = Cogent_icell8$qc_data$params$gm_log_base)
icell8_markers <- lapply(icell8_markers, function(obj){
  obj$rank <- -log10(obj$pval)*obj$log_fc
  obj[order(obj$rank, decreasing = T),]
})
names(icell8_markers) <- toupper(names(icell8_markers))

# Get markers for cellenion
cellenion_samples <- grep(pattern="Goe247_wo__fix", Cogent_cellenion$qc_data$metadata$Sample, invert=T)
cellenion_gm <- Cogent_cellenion$qc_data$gm[,cellenion_samples]
cellenion_metadata <- Cogent_cellenion$qc_data$metadata[cellenion_samples,]
cellenion_markers <- de.markers(gm = cellenion_gm,
                             metadata = cellenion_metadata,
                             min_fc = 0,
                             grouping_var = "Sample",
                             transform_type = Cogent_cellenion$qc_data$params$gm_log_base)
cellenion_markers <- lapply(cellenion_markers, function(obj){
  obj$rank <- -log10(obj$pval)*obj$log_fc
  obj[order(obj$rank, decreasing = T),]
})
names(cellenion_markers) <- toupper(names(cellenion_markers))

# Print out markers for each sample
for(sample in names(icell8_markers)){
  # ICELL8 markers
  write.table(icell8_markers[[sample]],
              file = paste("markers_icell8_",sample,".csv", sep=""),
              sep=",", quote = F, row.names=T, col.names = NA)
  
  # CellenIon markers
  write.table(cellenion_markers[[sample]],
              file = paste("markers_composite_",sample,".csv", sep=""),
              sep=",", quote = F, row.names=T, col.names = NA)
}
