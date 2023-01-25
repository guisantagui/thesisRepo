# This script takes as input all the deletionsDF csv files created in the by the 
# deleted_epitopes.py script and builds the allDels.csv file, with all the deletions
# in all strains with multiple identifiers

setwd("/storage/PGO/scripts/deletionsAnalysis/")
if(!require(KEGGREST)) BiocManager::install("KEGGREST", update = F)
library(KEGGREST)

delDFsPath <- "/storage/PGO/results/mtb/deletions_allGenes/"
dictDir <- "/storage/PGO/data/mtb/annotations/annotDict/"

allDelsOutDir <- "/storage/PGO/results/mtb/deletions_allGenes/"

fileNames <- list.files(delDFsPath)
fileNames <- fileNames[grep("deletionsDF", fileNames)]

allDels <- data.frame()
for(i in 1:length(fileNames)){
        delsLin <- read.csv(paste(delDFsPath, fileNames[i], sep = "/"))
        rownames(delsLin) <- delsLin[, 1]
        delsLin <- delsLin[, 2:ncol(delsLin)]
        lin <- gsub("nr", "", strsplit(fileNames[i], "_|\\.|-")[[1]][2])
        colnames(delsLin) <- paste(colnames(delsLin), lin, sep = "_")
        if(nrow(allDels) == 0){
                allDels <- delsLin
        }else{
                allDels <- cbind.data.frame(allDels, delsLin)
        }
}
print("allDels object generated")

keggGenes <- gsub("mtu:", "", names(keggList("mtu")))
# This function retrieves a set of identifiers from KEGG, given a vector of Rv codes
getGeneDBIDs <- function(RvCodes){
        ort <- c()
        ncbi_gID <- c()
        ncbi_pID <- c()
        uniProt <- c()
        for(i in 1:length(RvCodes)){
                #print(i)##Edit
                #print(RvCodes[i])##Edit
                if(RvCodes[i] %in% keggGenes){
                kGet <- keggGet(paste("mtu", RvCodes[i], sep = ":"))[[1]]
                        a <- kGet$ORTHOLOGY[1]
                        dbs <- kGet$DBLINKS
                        if(!is.null(a)){
                                ort <- c(ort, a)
                        }else{
                                ort <- c(ort, NA)
                        }
                        gID <- dbs[grep("GeneID", dbs)]
                        pID <- dbs[grep("ProteinID", dbs)]
                        up <- dbs[grep("UniProt", dbs)]
                        if(length(gID) != 0){
                                ncbi_gID <- c(ncbi_gID, gsub("NCBI-GeneID: ", "", gID))
                        }else{
                                ncbi_gID <- c(ncbi_gID, NA)
                        }
                        if(length(pID) != 0){
                                ncbi_pID <- c(ncbi_pID, gsub("NCBI-ProteinID: ", "", pID))
                        }else{
                                ncbi_pID <- c(ncbi_pID, NA)
                        }
                        if(length(up) != 0){
                                uniProt <- c(uniProt, gsub("UniProt: ", "", up))
                        }else{
                                uniProt <- c(uniProt, NA)
                        }
                }else{
                        ort <- c(ort, NA)
                        ncbi_gID <- c(ncbi_gID, NA)
                        ncbi_pID <- c(ncbi_pID, NA)
                        uniProt <- c(uniProt, NA)
                }
                
        }
        annotDF <- data.frame(ncbi_geneID = ncbi_gID,
                              ncbi_proteinID = ncbi_pID,
                              uniProt = uniProt,
                              orthology = ort)
        
        rownames(annotDF) <- RvCodes
        return(annotDF)
}

annotDF <- getGeneDBIDs(rownames(allDels))

# Load gene dictionary and add product to the annotation DF
h37RvGenesDict <- read.csv(paste(dictDir, "allGenesMTBC_dict.csv", sep = ""))

annotDF$product <- h37RvGenesDict$product[match(rownames(annotDF), h37RvGenesDict$locus_tag)]

ECnum <- sub("\\].*", "", sub(".*\\[", "", as.character(annotDF$orthology))) 

ECnum[!grepl("EC:", ECnum)] <- NA

annotDF$EC_number <- ECnum 

allDels <- cbind.data.frame(annotDF, allDels)

allDels <- cbind.data.frame(rownames(allDels), allDels)
colnames(allDels)[1] <- "gene"

write.csv(allDels, file = paste(allDelsOutDir, "allDels.csv", sep = ""), row.names = F)

print(paste("allDels.csv created at", paste(allDelsOutDir, "allDels.csv", sep = ""), sep = " "))