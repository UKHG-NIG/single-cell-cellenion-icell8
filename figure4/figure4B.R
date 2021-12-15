require(circlize)


goe247 <- read.csv("/home/maren/Dokumente/Julia_singleCell/documents/manuscript/markers/top100_markers_GOE247.txt", sep ="\t",
                   stringsAsFactors = FALSE, header = FALSE)
goe486 <- read.csv("/home/maren/Dokumente/Julia_singleCell/documents/manuscript/markers/top100_markers_GOE486.txt", sep ="\t",
                   stringsAsFactors = FALSE, header = FALSE)
goe615 <- read.csv("/home/maren/Dokumente/Julia_singleCell/documents/manuscript/markers/top100_markers_GOE615.txt", sep ="\t",
                   stringsAsFactors = FALSE, header = FALSE)
goe800 <- read.csv("/home/maren/Dokumente/Julia_singleCell/documents/manuscript/markers/top100_markers_GOE800.txt", sep ="\t",
                   stringsAsFactors = FALSE, header = FALSE)
goe1303 <- read.csv("/home/maren/Dokumente/Julia_singleCell/documents/manuscript/markers/top100_markers_GOE1303.txt", sep ="\t",
                    stringsAsFactors = FALSE, header = FALSE)
goe1305 <- read.csv("/home/maren/Dokumente/Julia_singleCell/documents/manuscript/markers/top100_markers_GOE1305.txt", sep ="\t",
                    stringsAsFactors = FALSE, header = FALSE)
goe1309 <- read.csv("/home/maren/Dokumente/Julia_singleCell/documents/manuscript/markers/top100_markers_GOE1309.txt", sep ="\t",
                    stringsAsFactors = FALSE, header = FALSE)
goe1360 <- read.csv("/home/maren/Dokumente/Julia_singleCell/documents/manuscript/markers/top100_markers_GOE1360.txt", sep ="\t",
                    stringsAsFactors = FALSE, header = FALSE)


df <- data.frame(Goe800=0,
                 Goe247=length(intersect(goe247$V1, goe800$V1)),
                 Goe486=length(intersect(goe486$V1, goe800$V1)),
                 Goe615=length(intersect(goe615$V1, goe800$V1)),
                 Goe1303=length(intersect(goe1303$V1, goe800$V1)),
                 Goe1305=length(intersect(goe1305$V1, goe800$V1)),
                 Goe1309=length(intersect(goe1309$V1, goe800$V1)),
                 Goe1360=length(intersect(goe1360$V1, goe800$V1)))

df[2,] <- c(length(intersect(goe247$V1, goe800$V1)), 0,
            length(intersect(goe247$V1, goe486$V1)),length(intersect(goe247$V1, goe615$V1)),
            length(intersect(goe247$V1, goe1303$V1)), length(intersect(goe247$V1, goe1305$V1)),
            length(intersect(goe247$V1, goe1309$V1)), length(intersect(goe247$V1, goe1360$V1)))

df[3,] <- c(length(intersect(goe486$V1, goe800$V1)), length(intersect(goe486$V1, goe247$V1)),
            0,length(intersect(goe486$V1, goe615$V1)),
            length(intersect(goe486$V1, goe1303$V1)), length(intersect(goe486$V1, goe1305$V1)),
            length(intersect(goe486$V1, goe1309$V1)), length(intersect(goe486$V1, goe1360$V1)))

df[4,] <- c(length(intersect(goe615$V1, goe800$V1)), length(intersect(goe615$V1, goe247$V1)),
            length(intersect(goe615$V1, goe486$V1)), 0,
            length(intersect(goe615$V1, goe1303$V1)), length(intersect(goe615$V1, goe1305$V1)),
            length(intersect(goe615$V1, goe1309$V1)), length(intersect(goe615$V1, goe1360$V1)))

df[5,] <- c(length(intersect(goe1303$V1, goe800$V1)), length(intersect(goe1303$V1, goe247$V1)),
            length(intersect(goe1303$V1, goe486$V1)), length(intersect(goe1303$V1, goe615$V1)), 
            0, length(intersect(goe1303$V1, goe1305$V1)),
            length(intersect(goe1303$V1, goe1309$V1)), length(intersect(goe1303$V1, goe1360$V1)))

df[6,] <- c(length(intersect(goe1305$V1, goe800$V1)), length(intersect(goe1305$V1, goe247$V1)),
            length(intersect(goe1305$V1, goe486$V1)), length(intersect(goe1305$V1, goe615$V1)), 
            length(intersect(goe1305$V1, goe1303$V1)), 0,
            length(intersect(goe1305$V1, goe1309$V1)), length(intersect(goe1305$V1, goe1360$V1)))

df[7,] <- c(length(intersect(goe1309$V1, goe800$V1)), length(intersect(goe1309$V1, goe247$V1)),
            length(intersect(goe1309$V1, goe486$V1)), length(intersect(goe1309$V1, goe615$V1)), 
            length(intersect(goe1309$V1, goe1303$V1)), length(intersect(goe1309$V1, goe1305$V1)),
            0, length(intersect(goe1309$V1, goe1360$V1)))

df[8,] <- c(length(intersect(goe1360$V1, goe800$V1)), length(intersect(goe1360$V1, goe247$V1)),
            length(intersect(goe1360$V1, goe486$V1)), length(intersect(goe1360$V1, goe615$V1)), 
            length(intersect(goe1360$V1, goe1303$V1)), length(intersect(goe1360$V1, goe1305$V1)),
            length(intersect(goe1360$V1, goe1309$V1)), 0)
rownames(df) <- colnames(df)
df

df2 = data.frame(from = rep(rownames(df), times = ncol(df)),
                 to = rep(colnames(df), each = nrow(df)),
                 value = as.vector(as.matrix(df)),
                 stringsAsFactors = FALSE)
df2


col_vec <- c("#FFD01580", "#5C5D6080", "#24A53B80", "#FF610380", "#DF1A2280", "#1C82E080", "#6E008280", "#005A7980")

names(col_vec) <- rownames(df)
jpeg("/zfstank/ngsdata/projects/icell8_cellenion_publication/figure4C_circulize_samples.jpeg", res = 150, width = 600, height = 600)
#svg("/zfstank/ngsdata/projects/icell8_cellenion_publication/figure4C_circulize_samples.svg")
chordDiagramFromDataFrame(df2, annotationTrack = c("name", "grid"), grid.col = col_vec)
dev.off()