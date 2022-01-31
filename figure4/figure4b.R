require(circlize)

# Retrieve files with top 100 markers 
topMarkersFiles <- list.files(pattern = "top100_markers*")
markersLists <- lapply(topMarkersFiles, function(file){
  tmp <- read.csv(file, sep =",", row.names = 1, stringsAsFactors = FALSE, header = T)
  rownames(tmp)
})
names(markersLists) <- gsub(".csv","",gsub("top100_markers_","",topMarkersFiles))

# Generate data.frame with number of top100 markers intersecting in each pair of samples
df <- as.data.frame(t(combn(names(markersLists),2)))
colnames(df) <- c("from","to")
df$value <- 0
for(i in 1:nrow(df)){
  sample1 <- df$from[i]
  sample2 <- df$to[i]
  df$value[i] <- length(intersect(markersLists[[sample1]],markersLists[[sample2]]))
}

# Generate Chord diagram
col_vec <- c(GOE800="#FFD01580",
             GOE247="#5C5D6080",
             GOE486="#24A53B80",
             GOE615="#FF610380",
             GOE1303="#DF1A2280",
             GOE1305="#1C82E080",
             GOE1309="#6E008280",
             GOE1360="#005A7980")

jpeg("figure4C_circulize_samples.jpeg", res = 150, width = 600, height = 600)
chordDiagramFromDataFrame(df, annotationTrack = c("name", "grid"), grid.col = col_vec)
graphics.off()
