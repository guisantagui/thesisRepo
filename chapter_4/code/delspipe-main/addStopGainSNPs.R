# Read btw snps file and add genes that have stop gain mutations in each lineage to 
# the lists of genes to remove from models. 
if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)

parser <- arg_parser("This script adds stop gain SNPs present in all the isolates within each lineage to the toRem files based on gene deletions.")

parser <- add_argument(parser = parser,
                       arg = c("--toRemDir", 
                               "--btwSNPFile"),
                       help = c("Directory where toRem files are stored",
                                "Path to the btwSNPs file"),
                       default = list("--toRemDir" = "/storage/PGO/results/mtb/deletions_analysis/toRemFromLinMods/",
                                      "--btwSNPFile" = "/storage/PGO/data/mtb/wrk_dataset/mutations/positions_btw_DR/btwlin_positions.annot"),
                       flag = c(F, F))

parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################
toRemFileDir <- parsed$toRemDir
snpFile <- parsed$btwSNPFile

# Load data, parse and get stop gain dataframe and lists per lineage. 
#########################################################################################################
print(sprintf("Parsing %s...", basename(snpFile)))
if(length(grep("snpeff", basename(snpFile))) >= 1){
        btwSNPs <- read.table(snpFile, 
                              sep = '\t', 
                              header = F,
                              row.names = NULL)
        colnames(btwSNPs)[2:11] <- c("POS", "DOT", "REF", "ALT", "LIN", "INFO", "CHANGE", "TYPE", "GENE", "AA")

        btwSNPs <- btwSNPs[, c("POS", "REF", "ALT", "LIN", "TYPE", "GENE", "AA")]
        btwSNPs$AA <- gsub("p.", "", btwSNPs$AA, fixed = T)
        btwSNPs$LOCUS <- sapply(btwSNPs$GENE, 
                                function(x) strsplit(as.character(x), "_")[[1]][length(strsplit(as.character(x), "_")[[1]])])
        geneNames <- c()
        for(i in 1:nrow(btwSNPs)){
                gene <- as.character(btwSNPs$GENE[i])
                locus <- btwSNPs$LOCUS[i]
                nme <- gsub(paste0("_", locus), "", gene)
                geneNames <- c(geneNames, nme)
        }
        btwSNPs$NAME <- geneNames
        btwSNPs$LIN[grep("a4", btwSNPs$LIN)] <- gsub("a4", "A4", btwSNPs$LIN[grep("a4", btwSNPs$LIN)])
        btwSNPs <- btwSNPs[, colnames(btwSNPs) != "GENE"]
        colnames(btwSNPs)[colnames(btwSNPs) == "LOCUS"] <- "GENE"
        stopGainSNPs <- btwSNPs[btwSNPs$TYPE == "stop_gained", ]
}else{
        btwSNPs <- read.table(snpFile, 
                              sep = '\t', 
                              header = T)
        btwSNPs$TYPE <- factor(gsub(" ", "", btwSNPs$TYPE))
        stopGainSNPs <- btwSNPs[btwSNPs$TYPE == "stopgain", ]
}

                      






lins <- unique(unlist(sapply(btwSNPs$LIN, function(x) strsplit(as.character(x), "_"))))

stopGList <- list()
for(i in 1:length(lins)){
        lin <- lins[i]
        linStopG <- gsub(" ", "", as.character(stopGainSNPs$GENE[grep(lin, stopGainSNPs$LIN)]))
        stopGList[[i]] <- linStopG
}
names(stopGList) <- lins
print(stopGList)
# Filter to remove lineages that don't have any stop gain snp. 
stopGList <- stopGList[unlist(lapply(stopGList, function(x) length(x) > 0))]
#print(stopGList)


# Read toRem files and add the stop-gained SNPs
#########################################################################################################
toRemFiles <- list.files(toRemFileDir)[grep(".txt", list.files(toRemFileDir))]

for(i in 1:length(toRemFiles)){
        fileName <- toRemFiles[i]
        linName <- gsub(".txt", "", gsub("toRem", "", fileName))
        if(linName %in% names(stopGList)){
                toRemGenes  <- as.character(read.table(paste(toRemFileDir, fileName, sep = ""), 
                                            sep = "\t", 
                                            header = F)$V1)
                stopGLin <- stopGList[[match(linName, names(stopGList))]]
                toRemGenes <- c(toRemGenes, 
                                stopGLin)
                toRemGenes <- unique(toRemGenes)
                write.table(toRemGenes, 
                            file = paste(toRemFileDir, fileName, sep = ""),
                            sep = "\t",
                            quote = F,
                            row.names = F,
                            col.names = F)
                print(paste("Genes with stop gain SNPs in ", 
                            linName,
                            " (",
                            paste(stopGLin, collapse = ", "), 
                            ") added to ",
                             fileName, 
                             sep = ""))
        }
}