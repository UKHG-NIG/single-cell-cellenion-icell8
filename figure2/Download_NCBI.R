# This script is used to download the GEO datasets from R
# Load libraries required to run the analyses
require(GEOquery)
require(RCurl)

# Set directory to which to download the files

# Download datasets from ICELL8
datasets <- c("GSM5411464","GSM5411465","GSM5411466")
for(data in datasets){
  gsm <- getGEO(dataset)
  ## Gene matrix
  url <- gsm@header$supplementary_file_1
  download.file(url, basename(url))
  gunzip(basename(url), overwrite=T, remove=T)
  ## Metadata
  url <- gsm@header$supplementary_file_2
  download.file(url, basename(url))
  gunzip(basename(url), overwrite=T, remove=T)
}
# Download geneInfo
gse <- getGEO("GSE179204", GSEMatrix = F)
url <- gsm@header$supplementary_file_1
download.file(url, basename(url))
gunzip(basename(url), overwrite=T, remove=T)

# Download datasets from ICELL8
gsm <- getGEO("GSM5411467")
## Gene matrix
url <- gsm@header$supplementary_file_1
download.file(url, basename(url))
gunzip(basename(url), overwrite=T, remove=T)
## Metadata
url <- gsm@header$supplementary_file_2
download.file(url, basename(url))
gunzip(basename(url), overwrite=T, remove=T)
## geneInfo
url <- gsm@header$supplementary_file_3
download.file(url, basename(url))
gunzip(basename(url), overwrite=T, remove=T)
