

gNumDir <- "/storage/PGO/scripts/deletions_allGenes/gnumbers/"
gNumFiles <- list.files(gNumDir)

lins <- gsub("nr.txt|.txt", "", gNumFiles)
lins <- gsub("\\-.*", "", lins)
lins <- sort(unique(lins))
lins <- lins[lins != "LNA" & lins != "README.md"]

for(i in seq_along(lins)){
        lin <- lins[i]
        linFiles <- gNumFiles[grep(lin, gNumFiles)]
        print(linFiles)
        gNumsLin <- c()
        for(j in seq_along(linFiles)){
                linFile <- linFiles[j]
                gNumFile <- read.table(paste0(gNumDir, linFile),sep = "\t")
                gNumFile <- as.character(gNumFile[, 1])
                gNumsLin <- c(gNumsLin, gNumFile)
        }
        write.table(gNumsLin, 
                    file = sprintf("/home/guisana/scripts/mapsPipe/data/%snr.txt", lin),
                    sep = "\t",
                    quote = F,
                    row.names = F,
                    col.names = F)
        print(sprintf("%snr.txt saved at /home/guisana/scripts/mapsPipe/data/", lin))
}
