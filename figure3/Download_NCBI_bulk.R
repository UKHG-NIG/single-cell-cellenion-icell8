# Load libraries required to run the analyses
require(GEOquery)
require(RCurl)

# Set directory to which to download the files

# Download gene count files related to bulk dataset
datasets <- 19:42
for(d in datasets){
  dataset <- paste("GSM54115",d,sep="")
  # gsm <- getGEO(dataset)
  ## Gene counts
  # url <- gsm@header$supplementary_file_1
  # download.file(url, basename(url))
  gunzip(basename(url), overwrite=T, remove=T)
}

# Create sample2group CSV file 
df <- as.data.frame(Reduce("rbind",lapply(list.files(pattern = "_geneCounts.txt"), function(file){
  group <- strsplit(basename(file), split = "_")[[1]][4]
  name <- paste(strsplit(basename(file), split = "_")[[1]][c(4,3)], collapse ="-")
  c(name, file, group)
})))
colnames(df) <- c("sampleName","fileName","group")
write.table(df, "bulk_table.csv", sep=",", row.names = F, col.names = T)
