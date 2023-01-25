##########################################################################################################
#                                                                                                        #
# Parses intra SNP files into a single matrix stating the proportion of strains within each lineage with #
# the altered nucleotide. Generate too a matrix where each column is s SNP and each row is a lineage.    #
# Cell values represent proportion of the isolates within the lineage with the given SNP.                #
#                                                                                                        #
##########################################################################################################
if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(rtracklayer)) install.packages("rtracklayer", repos='http://cran.us.r-project.org')
library(rtracklayer)
# Create parser, which accepts arguments from terminal and pass to the script
parser <- arg_parser("This script generates binarized versions of input deletion file, with the option of filtering the enzymatic genes and exporting a binarized set of them.")

parser <- add_argument(parser = parser,
                       arg = c("--btwSNPFile", 
                               "--intraSNPDir", 
                               "--annotDir",
                               "--output",
                               "--thrshldStpGa"),
                       help = c("Path to btw SNPs file",
                                "Path to directory where intra SNPs file are allocated", 
                                "Path to annotation directory",
                                "Output directory",
                                "Threshold of proportion of presence of a Stop Gain SNP to be written in Stop Gain .txt file (if --filtStopGa flag is present)"),
                       flag = c(F, F, F, F, F),
                       default = list("--btwSNPFile" = "/storage/PGO/data/mtb/wrk_dataset/mutations/positions_btw_DR/btw.snps.snpeff.annot",
                                      "--intraSNPDir" = "/storage/PGO/data/mtb/wrk_dataset/mutations/positions_intra_DR/",
                                      "--annotDir" = "/storage/PGO/data/mtb/annotations/",
                                      "--output" = "/storage/PGO/results/mtb/SNPs/",
                                      "--thrshldStpGa" = 0.85))
                                      
parser <- add_argument(parser = parser, 
                       arg = c("--filtNonSyn",
                               "--filtStopGa"),
                       help = c("Filters SNPs to get a dataframe of just non-synonymous SNPs",
                                "Filters SNPs to get .txt files of stop gain SNPs present within a lineage at above a given proportion, defined by --thrshldStpGa"),
                       flag = c(T, T))
                                      
parsed <- parse_args(parser)


#
# Directory stuff
##########################################################################################################
outDir <- parsed$output
intraSNPDir <- parsed$intraSNPDir
btwSNPFile <- parsed$btwSNPFile
annotDir <- parsed$annotDir

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

#
# Load and format data
##########################################################################################################

# Load annotation file and parse all the different annotations in a single .txt file 
grntAnnot <- as.data.frame(readGFF(paste(annotDir, "genes.gff", sep = "")))
ncbiAnnot <- as.data.frame(readGFF(paste(annotDir, "h37RvNCBIAnnot.gff3", sep = "")))
myc1Annot <- as.data.frame(readGFF(paste(annotDir, "Mycobacterium_tuberculosis_H37Rv_gff_v1.gff", sep = "")))
myc2Annot <- as.data.frame(readGFF(paste(annotDir, "Mycobacterium_tuberculosis_H37Rv_gff_v2.gff", sep = "")))

# Assign grntAnnot to annot, as it is the general one that we'll use in this script
annot <- grntAnnot


# Create file txt gene position file for running provean. 
genes <- grntAnnot

genes4txt <- genes[!is.na(genes$locus_tag), c("start", "end", "strand", "locus_tag")]

ncbiAnnot <- ncbiAnnot[1:nrow(ncbiAnnot) %% 2 == 1 & !is.na(ncbiAnnot$locus_tag), ]

ncbiAnnot2Add <- ncbiAnnot[!ncbiAnnot$locus_tag %in% genes4txt$locus_tag, c("start", "end", "strand", "locus_tag")]

genes4txt <- rbind.data.frame(genes4txt, ncbiAnnot2Add)

genes4txt <- genes4txt[order(genes4txt$start), ]

genes4txt$strand <- as.factor(sapply(genes4txt$strand, function(x) if(x == "+") x <- "F" else x <- "R"))

write.table(genes4txt, 
            file = sprintf("%sallGenesMTBC.txt", outDir), 
            sep = "\t", 
            quote = F, 
            row.names = F, 
            col.names = F)

print(sprintf("allGenesMTBC.txt saved at %s.", outDir))


# btw SNPs file. Better to use the SnpEff one...
print(sprintf("Parsing %s file...", basename(btwSNPFile)))
if(length(grep("snpeff", basename(btwSNPFile))) >= 1){
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
        #geneNames <- c()
        #for(i in 1:nrow(allBtwSNPs)){
        #        gene <- as.character(allBtwSNPs$GENE[i])
        #        locus <- allBtwSNPs$LOCUS[i]
        #        nme <- gsub(paste0("_", locus), "", gene)
        #        geneNames <- c(geneNames, nme)
        #}
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

# intraSNP files: parse them in a single dataframe with proportions within lineage with each SNP. Resaves the 
# file when finishes running a whole lineage, so in case that the running is interrupted something is saved, 
# and next time runs from what was done from beforehand.
                                                 
allIntraSNPs <- data.frame(matrix(ncol = 12, 
                                  nrow = 0, 
                                  dimnames = list(NULL, 
                                                  c("POS", 
                                                    "REF", 
                                                    "ALT", 
                                                    "LIN", 
                                                    "ANN", 
                                                    "GENE", 
                                                    "LOCUS", 
                                                    "NAME", 
                                                    "AA", 
                                                    "AA_one",
                                                    "propSNP",
                                                    "TOTAL"))))

intraSNPFiles <- list.files(intraSNPDir)
intraSNPFiles <- intraSNPFiles[grep(".annot", intraSNPFiles, fixed = T)]


outFiles <- list.files(outDir)

if("allIntraSNPs.csv" %in% outFiles){
        allIntraSNPs <- read.csv(sprintf("%sallIntraSNPs.csv", outDir), row.names = 1)
}

alreadyDoneLins <- unique(allIntraSNPs$LIN)

for(file in intraSNPFiles){
        lin <- strsplit(file, "_")[[1]][1]
        if(!lin %in% alreadyDoneLins){
                print(sprintf("Parsing %s intra SNPs file...", lin))
                tab <- read.table(sprintf("%s%s", intraSNPDir, file), sep = "\t", header = T, row.names = NULL)
                if("row.names" %in% colnames(tab)){
                        colNames <- colnames(tab)[2:ncol(tab)]
                        tab <- tab[, 1:(ncol(tab) - 1)]
                        colnames(tab) <- colNames
                        tab$lin <- rep(lin, nrow(tab))
                }
                colnames(tab)[grep("linage", colnames(tab))] <- "lin"
                for(i in seq_along(tab$counter)){
                        countTab <- data.frame(t(sapply(sapply(strsplit(as.character(tab$counter[i]), 
                                                                        "), (", 
                                                                        fixed = T)[[1]], 
                                                               function(x) gsub("[[:punct:]]", "", x)), 
                                                        function(x) strsplit(x, " ")[[1]])))
                        tabRow <- tab[i, ]
                        if(ncol(countTab) == 2){
                                colnames(countTab) <- c("nt", "counts")
                                rownames(countTab) <- 1:nrow(countTab)
                                countTab$counts <- as.numeric(as.character(countTab$counts))
                                snpsInPos <- strsplit(as.character(tab$HOMO_SNP[i]), ",")[[1]]
                                annoInPos <- strsplit(as.character(tab$annotation[i]), ",")[[1]]
                                aaChInPos <- strsplit(as.character(tab$aa.change[i]), ",")[[1]]
                                for(j in seq_along(snpsInPos)){
                                        snp <- snpsInPos[j]
                                        aaChan <- gsub("p.", "", aaChInPos[j], fixed = T)
                                        gene <- tabRow$Gene
                                        annoSNP <- annoInPos[j]
                                        #if(!is.na(gene)){
                                        #        geneSplit <- strsplit(as.character(gene), "_")[[1]]
                                        #        if(length(geneSplit) == 1){
                                        #                geneName <- geneSplit[1]
                                        #                geneLocus <- geneName
                                        #        }else{
                                        #                geneLocus <- geneSplit[length(geneSplit)]
                                        #                geneName <- gsub(paste0("_", geneLocus), "", gene)
                                        #        }
                                        ##}else if(annoSNP == "intergenic_region" & !is.na(annoSNP)){ # Edit
                                        #}else if(annoSNP == "intergenic_region"){                    # Edit --> to get the format gene1-gene2 in intergen regs' locus col in all cases.
                                        #        tabRow$Gene <- getIntrGens(as.numeric(as.character(tabRow$Position)))
                                        #        geneName <- tabRow$Gene
                                        #        geneLocus <- tabRow$Gene
                                        #}else{
                                        #        geneName <- NA
                                        #        geneLocus <- NA
                                        #}
                                        geneName <- getGeneName(gene)
                                        geneLocus <- getGeneLocus(gene)
                                        rowNewDF <- tabRow
                                        rowNewDF$lin <- lin ## Edit. 
                                        rowNewDF$HOMO_SNP <- snp
                                        rowNewDF$annotation <- annoSNP
                                        rowNewDF$aa.change <- aaChan
                                        propSNP <- countTab$counts[match(snp, countTab$nt)]/sum(countTab$counts)
                                        totLin <- sum(countTab$counts)
                                        rowNewDF$propSNP <- propSNP
                                        rowNewDF$totLin <- totLin
                                        rowNewDF$AA_one <- three2one(aaChan)
                                        rowNewDF$LOCUS <- geneLocus
                                        rowNewDF$NAME <- geneName
                                        rowNewDF <- rowNewDF[, c("Position", 
                                                                 "REF", 
                                                                 "HOMO_SNP", 
                                                                 "lin", 
                                                                 "annotation", 
                                                                 "Gene", 
                                                                 "LOCUS", 
                                                                 "NAME", 
                                                                 "aa.change", 
                                                                 "AA_one", 
                                                                 "propSNP", 
                                                                 "totLin")]
                                        colnames(rowNewDF) <- c("POS", 
                                                                "REF", 
                                                                "ALT", 
                                                                "LIN", 
                                                                "ANN", 
                                                                "GENE", 
                                                                "LOCUS", 
                                                                "NAME", 
                                                                "AA", 
                                                                "AA_one", 
                                                                "propSNP", 
                                                                "TOTAL")
                                        allIntraSNPs <- rbind.data.frame(allIntraSNPs, rowNewDF)
                                }
                        }
                }
                write.csv(allIntraSNPs, file = sprintf("%sallIntraSNPs.csv", outDir))
                print(sprintf("%s intra SNPs parsed and added to allIntraSNPs.csv.", lin))
        }else{
                print(sprintf("%s intra SNP file already parsed.", lin))
        }
}

allIntraSNPs$LOCUS <- sapply(allIntraSNPs$GENE, getGeneLocus)
allIntraSNPs$NAME <- sapply(allIntraSNPs$GENE, getGeneName)

#allIntraSNPs$LOCUS[allIntraSNPs$ANN == "intergenic_region"] <- sapply(allIntraSNPs$POS[allIntraSNPs$ANN == "intergenic_region"], # Run once, after that comment...
#                                                                      getIntrGens)

#allIntraSNPs$LOCUS[is.na(allIntraSNPs$LOCUS)] <- allIntraSNPs$NAME[is.na(allIntraSNPs$LOCUS)]
write.csv(allIntraSNPs, file = sprintf("%sallIntraSNPs.csv", outDir))                                                            # Run once, after that comment...

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
allSNPs$snpCode <- paste0(allSNPs$REF, allSNPs$POS, allSNPs$ALT)

# Remove rows empty in lineages (which are 4, caused because of 4 SNPs in the intergenic region of 
# Rv3924c that don't have altered nucleotide. These SNPs are very infrequent, so should not
# be a problem. Reorder mixed dataframe by lineage.  

allSNPs <- allSNPs[allSNPs$LIN != "", ]
allSNPs <- allSNPs[order(as.character(allSNPs$LIN)), ]
rownames(allSNPs) <- 1:nrow(allSNPs)

print(head(allSNPs))

write.csv(allSNPs, file = sprintf("%sallSNPs.csv", outDir))
print(sprintf("allSNPs.csv saved at %s.", outDir))

# Filter non synonymous, or missense, for PROVEAN and/or SIFT
if(parsed$filtNonSyn == T){
        print("Filtering missense SNPs...")
        allSNPsNonSyn <- allSNPs[allSNPs$ANN == "missense_variant" & !is.na(allSNPs$ANN), ]
        write.csv(allSNPsNonSyn, file = sprintf("%sallSNPsNonSyn.csv", outDir))
        print(sprintf("allSNPsNonSyn.csv saved at %s.", outDir))
}

# Generate .txt files with lists of genes that have stop gained SNPs at more than a given proportion 
# within a lineage
if(parsed$filtStopGa == T){
        stopGaOutDir <- paste0(outDir, "stopGained/")
        if(!dir.exists(stopGaOutDir)){
                dir.create(stopGaOutDir)
        }
        print(sprintf("Filtering stop gained SNPs... Threshold for considering stop gain SNP in lineage: %s%%.", 
                      as.character(round(parsed$thrshldStpGa*100, digits = 2))))
        allSNPsStopGa <- allSNPs[allSNPs$ANN == "stop_gained" & allSNPs$propSNP >= parsed$thrshldStpGa, ]
        lins <- unique(as.character(allSNPsStopGa$LIN))
        lins <- lins[!is.na(lins)]
        for(l in lins){
                stopGa_genes <- allSNPsStopGa[allSNPsStopGa$LIN == l, c("LOCUS", "AA_one", "propSNP")]
                stopGa_genes <- stopGa_genes[!is.na(stopGa_genes$LOCUS), ]
                write.table(stopGa_genes,
                            quote = F,
                            sep = "\t",
                            row.names = F,
                            col.names = T,  
                            file = sprintf("%s%s_stopGain.txt", stopGaOutDir, l))
                print(sprintf("%s_stopGain.txt saved at %s.", l, stopGaOutDir))
        }
}
# Build a matrix where rows are lineages and columns unique SNPs and each value corresponds to the 
# proportion of isolates within each lineage carrying the SNP
##########################################################################################################

#lineages <- sapply(intraSNPFiles, function(x) strsplit(x, "_")[[1]][1])

#allUniqSNPs <- unique(allSNPs$snpCode)
                        
#allSNPsMat <- matrix(rep(0, length(lineages)*length(allUniqSNPs)), 
#                     nrow = length(lineages), 
#                     ncol = length(allUniqSNPs),
#                     dimnames = list(lineages, allUniqSNPs))

#print(rownames(allSNPsMat))
#print("Building all SNP matrix...")
                
#for(i in 1:ncol(allSNPsMat)){
#        snp <- colnames(allSNPsMat)[i]
#        print(snp)
##        snpSubMat <- allSNPs[allSNPs$snpCode == snp, c("LIN", "propSNP")]
#        print(snpSubMat)
##        print(snpSubMat$LIN)
#        print(match(snpSubMat$LIN, rownames(allSNPsMat)))
#        print(allSNPsMat[match(snpSubMat$LIN, rownames(allSNPsMat)), snp])
#        allSNPsMat[match(snpSubMat$LIN, rownames(allSNPsMat)), snp] <- as.numeric(snpSubMat$propSNP)
#}

#write.csv(allSNPsMat, file = sprintf("%sallSNPsMat.csv", outDir))

#print(sprintf("allSNPsMat.csv saved at %s.", outDir))