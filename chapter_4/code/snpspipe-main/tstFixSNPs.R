if(!require(rtracklayer)) install.packages("rtracklayer", repos='http://cran.us.r-project.org')
library(rtracklayer)


allIntraSNPsDir <- "/storage/PGO/results/mtb/SNPs/"
outDir <- "/storage/PGO/results/mtb/SNPs/"
annotDir <- "/storage/PGO/data/mtb/annotations/"
btwSNPFile = "/storage/PGO/data/mtb/wrk_dataset/mutations/positions_btw_DR/btw.snps.snpeff.annot"

#
# Functions needed in script: 
##########################################################################################################

# Translates 3 amino acid substitution code to 1 aminoacid.  
three2one <- function(change){
        aaDict <- data.frame(AA = c("alanine",
                                    "arginine",
                                    "asparagine",
                                    "aspartic acid",
                                    "asparagine or aspartic acid",
                                    "cysteine",
                                    "glutamic acid",
                                    "glutamine",
                                    "glutamine or glutamic acid",
                                    "glycine",
                                    "histidine",
                                    "isoleucine",
                                    "leucine",
                                    "lysine",
                                    "methionine",
                                    "phenylalanine",
                                    "proline",
                                    "serine",
                                    "threonine",
                                    "tryptophan",
                                    "tyrosine",
                                    "valine"), 
                             three = c("ala",
                                       "arg",
                                       "asn",
                                       "asp",
                                       "asx",
                                       "cys",
                                       "glu",
                                       "gln",
                                       "glx",
                                       "gly",
                                       "his",
                                       "ile",
                                       "leu",
                                       "lys",
                                       "met",
                                       "phe",
                                       "pro",
                                       "ser",
                                       "thr",
                                       "trp",
                                       "tyr",
                                       "val"),
                             one = c("A",
                                     "R",
                                     "N",
                                     "D",
                                     "B",
                                     "C",
                                     "E",
                                     "Q",
                                     "Z",
                                     "G",
                                     "H",
                                     "I",
                                     "L",
                                     "K",
                                     "M",
                                     "F",
                                     "P",
                                     "S",
                                     "T",
                                     "W",
                                     "Y",
                                     "V"))
        if(grepl("*", change, fixed = T)){
                frst <- tolower(substr(change, 1, 3))
                scnd <- "*"
                num <- substr(change, 4, nchar(change)-1)
                out <- paste(aaDict$one[match(frst, aaDict$three)],
                             num,
                             scnd,
                             sep = "")
        }else if(is.na(change)){
                out <- change
        }else if(nchar(change) >= 7){
                frst <- tolower(substr(change, 1, 3))
                scnd <- tolower(substr(change, nchar(change)-2, nchar(change)))
                num <- substr(change, 4, nchar(change)-3)
                out <- paste(aaDict$one[match(frst, aaDict$three)],
                             num,
                             aaDict$one[match(scnd, aaDict$three)],
                             sep = "")
        }else{
                out <- change
        }
        return(out)
}

# Given a certain SNP position located in an intergenic region outputs a string indicating left and right 
# locus tags of corresponding genes.
getIntrGens <- function(position){
        position <- as.numeric(as.character(position))
        left <- annot$locus_tag[max(which(annot$end < position))-1]
        right <- annot$locus_tag[min(which(annot$start > position))]
        out <- paste(left, right, sep = "-")
        return(out)
}


annot <- as.data.frame(readGFF(paste(annotDir, "genes.gff", sep = "")))

# Given the gene string of SNP eff including gene name and locus tag isolates locus tag.
# If there are more than one (intergenic regions) collapses all of them in a single string
getGeneLocus <- function(gene){
        if(!is.na(gene)){
                gene <- gsub(" ", "", gene)
                gene <- gsub("-", "_", gene)
                geneSplit <- strsplit(gene, "_")[[1]]
                geneSplit <- unique(geneSplit[grep("rv", tolower(geneSplit))])
                if(length(geneSplit) > 1){
                        locus <- paste(geneSplit, collapse = "-")
                }else{
                        locus <- geneSplit[1]
                }
        }else{
                locus <- NA
        }
        return(locus)
}

# Given the gene string of SNP eff including gene name and locus tag isolates gene name.
# If there are more than one (intergenic regions) collapses all of them in a single string
getGeneName <- function(gene){
        if(!is.na(gene)){
                gene <- gsub(" ", "", gene)
                gene <- gsub("-", "_", gene)
                geneSplit <- strsplit(gene, "_")[[1]]
                geneSplit <- unique(geneSplit[-grep("rv", tolower(geneSplit))])
                if(length(geneSplit) > 1){
                        name <- paste(geneSplit, collapse = "_")
                }else{
                        name <- geneSplit[1]
                }
        }else{
                name <- NA
        }
        return(name)
}

# btw SNPs file. Better to use the SnpEff one...
print(sprintf("Parsing %s file...", basename(btwSNPFile)))
if(length(grep("snpeff", basename(btwSNPFile))) >= 1){
        print("Ja!")
        allBtwSNPs <- read.table(btwSNPFile, 
                                 sep = '\t', 
                                 fill = T,
                                 header = F,
                                 row.names = NULL)
        colnames(allBtwSNPs)[2:11] <- c("POS", "DOT", "REF", "ALT", "LIN", "INFO", "CHANGE", "ANN", "GENE", "AA")

        allBtwSNPs <- allBtwSNPs[, c("POS", "REF", "ALT", "LIN", "ANN", "GENE", "AA")]
        allBtwSNPs$ANN <- as.character(allBtwSNPs$ANN)
        allBtwSNPs$AA <- gsub("p.", "", allBtwSNPs$AA, fixed = T)
        #allBtwSNPs$LOCUS <- sapply(allBtwSNPs$GENE, function(x) strsplit(as.character(x), "_")[[1]][length(strsplit(as.character(x), "_")[[1]])])
        allBtwSNPs$LOCUS <- sapply(allBtwSNPs$GENE, getGeneLocus)
        geneNames <- c()
        for(i in 1:nrow(allBtwSNPs)){
                gene <- as.character(allBtwSNPs$GENE[i])
                locus <- allBtwSNPs$LOCUS[i]
                nme <- gsub(paste0("_", locus), "", gene)
                geneNames <- c(geneNames, nme)
        }
        #allBtwSNPs$NAME <- geneNames
        allBtwSNPs$NAME <- sapply(allBtwSNPs$GENE, getGeneName)
        allBtwSNPs$AA_one <- sapply(allBtwSNPs$AA, three2one)
        allBtwSNPs$LIN[grep("a4", allBtwSNPs$LIN)] <- gsub("a4", "A4", allBtwSNPs$LIN[grep("a4", allBtwSNPs$LIN)])
        #allBtwSNPs$LOCUS[allBtwSNPs$ANN == "intergenic_region"] <- sapply(allBtwSNPs$POS[allBtwSNPs$ANN == "intergenic_region"], getIntrGens)
}else{
        allBtwSNPs <- read.table(btwSNPFile, 
                                 sep = '\t', 
                                 header = T)
        allBtwSNPs$TYPE <- gsub(" ", "", as.character(btwSNPs$TYPE))
        colnames(allBtwSNPs) <- c("POS", "REF", "ALT", "ANN", "LOCUS", "AA_one", "LIN")
        allBtwSNPs$ANN[allBtwSNPs$ANN == "nonsynonymous"] <- gsub("nonsynonymous", 
                                                                  "missense_variant", 
                                                                  allBtwSNPs$ANN[allBtwSNPs$ANN == "nonsynonymous"])
        allBtwSNPs$ANN[allBtwSNPs$ANN == "synonymous"] <- gsub("synonymous", 
                                                               "synonymous_variant", 
                                                               allBtwSNPs$ANN[allBtwSNPs$ANN == "synonymous"])
        allBtwSNPs$ANN[allBtwSNPs$ANN == "stopgain"] <- gsub("stopgain", 
                                                             "stop_gained", 
                                                             allBtwSNPs$ANN[allBtwSNPs$ANN == "stopgain"])
        allBtwSNPs$ANN[allBtwSNPs$ANN == "stoplost"] <- gsub("stoplost", 
                                                             "stop_lost&splice_region_variant", 
                                                             allBtwSNPs$ANN[allBtwSNPs$ANN == "stoplost"])
        allBtwSNPs$NAME <- annot$gene[match(allBtwSNPs, annot$locus_tag)]
}

allBtwSNPs <- allBtwSNPs[, c("POS", "REF", "ALT", "LIN", "ANN", "LOCUS", "NAME", "AA_one")]

allBtwSNPs$propSNP <- rep(1, nrow(allBtwSNPs))
allBtwSNPs$TOTAL <- rep(NA, nrow(allBtwSNPs))

# Convert all btwSNPs dataframe in another dataframe where in the columns "LIN" there is only a single 
# lineage-->this dataframe will have more rows. 
allBtwSNPs_new <- data.frame(matrix(ncol = 10, 
                                    nrow = 0, 
                                    dimnames = list(NULL, 
                                                    c("POS", 
                                                      "REF", 
                                                      "ALT", 
                                                      "LIN", 
                                                      "ANN", 
                                                      "LOCUS", 
                                                      "NAME", 
                                                      "AA_one",
                                                      "propSNP",
                                                      "TOTAL"))))

for(i in 1:nrow(allBtwSNPs)){
        inLins <- strsplit(as.character(allBtwSNPs$LIN[i]), "_")[[1]]
        rowNew <- allBtwSNPs[i, ]
        for(j in seq_along(inLins)){
                lin <- inLins[j]
                rowNew$LIN <- lin
                allBtwSNPs_new <- rbind.data.frame(allBtwSNPs_new, rowNew)
        }
}

print(sprintf("%s file parsed.", basename(btwSNPFile)))

print("Reading allIntraSNPs.csv...")
allIntraSNPs <- read.csv(paste0(allIntraSNPsDir, "allIntraSNPs.csv"), row.names = 1)

allIntraSNPs$NAME <- sapply(allIntraSNPs$GENE, getGeneName)
allIntraSNPs$LOCUS <- sapply(allIntraSNPs$GENE, getGeneLocus)


intergenLocs <- sapply(allIntraSNPs$POS[allIntraSNPs$ANN == "intergenic_region"], # Run once, after that comment...
                                                                      getIntrGens)
                                                                      
print(intergenLocs)

print(which(is.na(intergenLocs)))

print(allIntraSNPs$LOCUS[allIntraSNPs$ANN == "intergenic_region"])

# allIntraSNPs$LOCUS[allIntraSNPs$ANN == "intergenic_region"] <- intergenLocs

print("Du!")
                                                                     
write.csv(allIntraSNPs, file = sprintf("%sallIntraSNPs_2.csv", allIntraSNPsDir))
print("All intra SNP files parsed!")

# Mix intraSNPs and btwSNPs in a single dataframe and export it. If flags present export too filtered
# nonsynonymous SNPs set and stop gain SNPs .txt files. 
##########################################################################################################

allBtwSNPs_new$LIN <- as.character(allBtwSNPs_new$LIN)
allIntraSNPs$LIN <- as.character(allIntraSNPs$LIN)

# Remove columns of allIntraSNPs that are not present in allBtwSNPs_new.
allIntraSNPs <- allIntraSNPs[,  c("POS", 
                                  "REF", 
                                  "ALT", 
                                  "LIN", 
                                  "ANN", 
                                  "LOCUS", 
                                  "NAME",
                                  "AA_one", 
                                  "propSNP", 
                                  "TOTAL")]


allSNPs <- rbind.data.frame(allBtwSNPs_new, allIntraSNPs)
allSNPs$snpCode <- paste(allSNPs$LOCUS, paste0(allSNPs$REF, allSNPs$POS, allSNPs$ALT), sep = ".")

# Remove rows empty in lineages (which are 4, caused because of 4 SNPs in the intergenic region of 
# Rv3924c that don't have altered nucleotide. These SNPs are very infrequent, so should not
# be a problem. Reorder mixed dataframe by lineage.  

allSNPs <- allSNPs[allSNPs$LIN != "", ]
allSNPs <- allSNPs[order(as.character(allSNPs$LIN)), ]
rownames(allSNPs) <- 1:nrow(allSNPs)

write.csv(allSNPs, file = sprintf("%sallSNPs.csv", outDir))
print(sprintf("allSNPs.csv saved at %s.", outDir))
                                                                      







# Build a matrix where rows are lineages and columns unique SNPs and each value corresponds to the 
# proportion of isolates within each lineage carrying the SNP
##########################################################################################################

lineages <- c("A1", 
              "A2", 
              "A3", 
              "A4", 
              "L1", 
              "L2", 
              "L3", 
              "L4", 
              "L5", 
              "L6", 
              "L7", 
              "L8", 
              "L9")

allUniqSNPs <- unique(allSNPs$snpCode)
                        
allSNPsMat <- matrix(rep(0, length(lineages)*length(allUniqSNPs)), 
                     nrow = length(lineages), 
                     ncol = length(allUniqSNPs),
                     dimnames = list(lineages, allUniqSNPs))

print(rownames(allSNPsMat))
print("Building all SNP matrix...")
                
for(i in 1:ncol(allSNPsMat)){
        snp <- colnames(allSNPsMat)[i]
        #print(snp)
        snpSubMat <- allSNPs[allSNPs$snpCode == snp, c("LIN", "propSNP")]
        #print(snpSubMat)
        #print(snpSubMat$LIN)
        #print(match(snpSubMat$LIN, rownames(allSNPsMat)))
        #print(allSNPsMat[match(snpSubMat$LIN, rownames(allSNPsMat)), snp])
        allSNPsMat[match(snpSubMat$LIN, rownames(allSNPsMat)), snp] <- as.numeric(snpSubMat$propSNP)
}

write.csv(allSNPsMat, file = sprintf("%sallSNPsMat_2.csv", outDir))

print(sprintf("allSNPsMat_2.csv saved at %s.", outDir))                                                                     
                                                                      
                                                                      