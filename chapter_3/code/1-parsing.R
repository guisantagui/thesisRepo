####################################################################################################
####################################################################################################
#                                                                                                  #
# Parse the Liquid Chromatography-Mass Spectrometry metabolomic data to a format appropriate for   #
# "metabolomics" and "NormalizeMets" R packages. Add to the peaks data the added metabolites in a  #
# newest set where more metabolites where identified and unifies the nomenclature of both.         #
#                                                                                                  #
####################################################################################################
####################################################################################################

if(!require(readxl)) install.packages("readxl")
library(readxl)
if(!require(dplyr)) install.packages("dplyr")
library(dplyr)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)

# Directory stuff
####################################################################################################
rootDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/PA_paper_scripts/" ## Change this 
                                                                                      ## path 
                                                                                      ## according 
                                                                                      ## user's 
                                                                                      ## direcory
                                                                                      ## location.
dataDir <- paste0(rootDir, "data/")
outDir <- paste0(rootDir, "results/")
plotDir <- paste0(rootDir, "plots/")

# Create directories if necessary
if(!dir.exists(plotDir)){
        dir.create(plotDir)
}
if(!dir.exists(outDir)){
        dir.create(outDir)
}

# Functions
####################################################################################################

# KEGGREST functions are giving problems, so let's build custom functions that do the same
keggFind_cust <- function(query, database, option){
        if(!require(RCurl)) install.packages("RCurl")
        library(RCurl)
        URL <- sprintf("https://rest.kegg.jp/find/%s/%s/%s", database, query, option)
        tab <- getURL(URL)
        if(nchar(tab) > 1){
                outDF <- read.table(URL, quote="", sep="\t")
                out <- as.character(outDF[, 2])
                names(out) <- outDF[, 1]
        }else{
                out <- NULL
        }
        
        return(out)
}

keggLink_cust <- function(target, source){
        if(!require(RCurl)) install.packages("RCurl")
        library(RCurl)
        URL <- sprintf("https://rest.kegg.jp/link/%s/%s", target, source)
        tab <- getURL(URL)
        if(nchar(tab) > 1){
                outDF <- read.table(URL, quote="", sep="\t")
                out <- as.character(outDF[, 2])
                names(out) <- outDF[, 1]
        }else{
                out <- NULL
        }
        
        return(out)
}

keggList_cust <- function(database){
        if(!require(RCurl)) install.packages("RCurl")
        library(RCurl)
        URL <- sprintf("https://rest.kegg.jp/list/%s", database)
        tab <- getURL(URL)
        if(nchar(tab) > 1){
                outDF <- read.table(URL, quote="", sep="\t")
                out <- as.character(outDF[, 2])
                names(out) <- outDF[, 1]
        }else{
                out <- NULL
        }
        
        return(out)
}
####################################################################################################
#                                                                                                  #
# Load the data.                                                                                   #
#                                                                                                  #
####################################################################################################

# Dictionary of metabolite names between different metabolite sets
dictionary <- read.csv(file = paste0(dataDir, "dictionary.csv"), 
                       row.names = 1, 
                       stringsAsFactors = F)

# Dictionary of enzymatic genes present in each strain. 
dictEnzyms <- read.csv(file = paste0(dataDir, "dictEnzymes.csv"), 
                       row.names = 1,
                       stringsAsFactors = F)

####################################################################################################
# Peaks from old set, which has several putative compounds that need to be handled. 
####################################################################################################

#Load Names
sampleNames <- read_excel(paste0(dataDir, "metData/sampleName1.xlsx"))
colnames(sampleNames) <- paste(sampleNames[1, ], sampleNames[2, ], sep = '_')

#Load data
negativePeaks <- read_excel(paste0(dataDir, "metData/ClinicalNeg-2.xlsx"))
colnames(negativePeaks)[6:89] <- colnames(sampleNames)
negativePeaks$negativePeaks = 1
negativePeaks$mode <- "negative"

positivePeaks1 <- read_excel(paste0(dataDir, "metData/ClinicalPos-4.1.xlsx"))
colnames(positivePeaks1)[6:89] <- colnames(sampleNames)
positivePeaks1$positivePeaks1 = 1
positivePeaks1$mode <- "positive1"

positivePeaks2 <- read_excel(paste0(dataDir, "metData/ClinicalPos-4.2.xlsx"))
colnames(positivePeaks2)[6:89] <- colnames(sampleNames)
positivePeaks2$positivePeaks2 = 1
positivePeaks2$mode <- "positive2"

metabolitePeaks_old <- merge(negativePeaks, positivePeaks1, all = T)
metabolitePeaks_old <- merge(metabolitePeaks_old, positivePeaks2, all = T)
metabolitePeaks_old[is.na(metabolitePeaks_old)] <- 0

#Drop variables that we don't want
drops <- c("CAS ID", "negativePeaks", "positivePeaks1", "positivePeaks2")
metabolitePeaks_old <- metabolitePeaks_old[, !(names(metabolitePeaks_old) %in% drops)]


#Convert to format the data
sampleNames <- t(sampleNames)
sampleNames <- sampleNames[, c(1, 2)]
colnames(sampleNames) <- c("Group", "Sample")
sampleNames <- data.frame(sampleNames[, c("Sample", "Group")])

rownames(metabolitePeaks_old) <- metabolitePeaks_old[, 1]
metFormat <- data.frame(t(metabolitePeaks_old[, -c(1, 2, 3, 4, 89)]))
metFormat <- merge(sampleNames, metFormat, by = "row.names",
                                   all.x = T, sort = F)
metFormat <- metFormat[, -c(1, 2)]
row.names(metFormat) <- row.names(sampleNames)
colnames(metFormat)[2:ncol(metFormat)] <- metabolitePeaks_old$`Compound Name`

#Replace missing values/zeros
metFormat[metFormat == 0] <- NA


# Change metabolite names to consensus. 
colnames(metFormat)[2:ncol(metFormat)] <- as.character(dictionary$Consensus[match(colnames(metFormat)[2:ncol(metFormat)], 
                                                                                  dictionary$Old.Data.Names)])

# Remove the columns which new name is NA, as in the dictionary are the compounds that are redundant 
# (have same values)
metFormat <- metFormat[, !is.na(colnames(metFormat))]


# Check if we can rescue any of the putative compounds based on the molecular weights of the 
# peaks and the enzymatic content of our clinical isolates. 
####################################################################################################

# Get the molecular weight of the metabolites in metFormat data frame
metPckgFrmt_oldNames <- dictionary$Old.Data.Names[match(colnames(metFormat)[2:ncol(metFormat)],
                                                                 dictionary$Consensus)]

molWght <- metabolitePeaks_old$Mass[match(metPckgFrmt_oldNames, metabolitePeaks_old$`Compound Name`)]
names(molWght) <- colnames(metFormat)[2:ncol(metFormat)]

# Filter putative compounds (the ones with the "_?" tag)
putMolWght <- molWght[grep("_?", names(molWght), fixed = T)]

# Obtain all the existing compounds with these molecular weights from KEGG
#possCpds <- sapply(putMolWght, keggFind, database = "compound", option = "exact_mass")
possCpds <- sapply(putMolWght, keggFind_cust, database = "compound", option = "exact_mass") # Edit
possCpds <- possCpds[sapply(possCpds, function(x) length(x) != 0)]

# Obtain all the enzymes that produce or consume each possible candidate compound (with same) 
# molecular weight that our peak). Examine what of these enzymes are in prokka annotation. 
# Count the number of enzymes related to each candidate compound appears in the list of
# our annotated enzymes. 
enzListsPutComp <- list()
for(i in seq_along(possCpds)){
        enzsPossComps <- list()
        for(j in seq_along(possCpds[[i]])){
                #enzs <- keggLink(target = "enzyme", source = gsub("cpd:", 
                #                                                  "", 
                #                                                  names(possCpds[[i]])[j]))
                enzs <- keggLink_cust(target = "enzyme", source = gsub("cpd:", 
                                                                       "", 
                                                                       names(possCpds[[i]])[j]))
                enzs <- gsub("ec:", 
                             "", 
                             enzs)
                inAnnot <- enzs %in% dictEnzyms$ECnums
                names(inAnnot) <- enzs
                inAnnot <- c(inAnnot, sum(inAnnot))
                names(inAnnot)[length(inAnnot)] <- "sums"
                enzsPossComps[[j]] <- inAnnot
        }
        names(enzsPossComps) <- names(possCpds[[i]])
        enzListsPutComp[[i]] <- enzsPossComps
}
names(enzListsPutComp) <- names(possCpds)

# Identify, for each of the ambiguous compounds, the candidates that we have genetic evidence of 
# their presence.  
# Keep the candidates that, based on annotation, are unique: the ones that, among the possible 
# compounds, our isolates only have enzymes related to it. 
uniqueCandidates <- unlist(lapply(enzListsPutComp, function(x){
        sums <- sapply(x, function(y) y[length(y)]) > 0
        if(sum(sums) <= 1){
                return(names(x)[sums])
        }
})) 

# Get the KEGG names of the unique candidates
#uniqueMetsNames <- keggList("compound")[names(keggList("compound" )) %in% uniqueCandidates]
uniqueMetsNames <- keggList_cust("compound")[names(keggList_cust("compound" )) %in% uniqueCandidates]
uniqueMetsNames <- sapply(uniqueMetsNames, function(x) unlist(strsplit(x, ";"))[1])
names(uniqueMetsNames) <- names(uniqueCandidates)[match(names(uniqueMetsNames), uniqueCandidates)]
solvedAmbigMets <- uniqueMetsNames

# The unique candidate compound for Ethyl malate 2_? is cpd:C03979, or 2-Dehydro-3-deoxy-L-rhamnonate,
# which already is in the metabolomic dataset, so removed it from the solved putative compounds
# vector
solvedAmbigMets <- solvedAmbigMets[names(solvedAmbigMets) != "Ethyl malate 2_?"]

# Change names of solved putative mets. 
colnames(metFormat)[match(names(solvedAmbigMets), 
                                          colnames(metFormat))] <- solvedAmbigMets

# create a new metabolite dictionary with the names of the solved putative compounds changed, 
# to use it in further analysis. 

dictionary$Consensus[match(names(solvedAmbigMets), dictionary$Consensus)] <- solvedAmbigMets

write.csv(dictionary, file = paste0(outDir, "dictionary_solvedPut.csv"))

####################################################################################################
# Peaks from newest table, which has some added compounds and other removed because were putative. 
# We're keeping them because we use them in the normalization, and remove them afterwards. Just add
# the ones labeled as "Added" in the metabolite dictionary. 
####################################################################################################
mets <-read_xlsx(paste0(dataDir, "metData/metaboliteTable.xlsx"))

strain_names <- unique(mets$mutant)
sampleNames <- mets$sample %>% gsub("Negative|Positive", "", .) %>% gsub("day1", "_", .) %>% unique()


positives <- mets[which(mets$mode == "Positive"), ]
posMat <- NULL
posBatches <- NULL
for(i in 1:length(strain_names)){
        for(j in 1:3){
                if(is.null(posMat) == T){
                        posMat <- positives[positives$sample == paste(strain_names[i], "day1PositiveR", j, sep = ""), 1]
                        posBatches <- positives[positives$sample == paste(strain_names[i], "day1PositiveR", j, sep = ""), 9]
                }else{
                        posMat <- cbind(posMat, positives[positives$sample == paste(strain_names[i], "day1PositiveR", j, sep = ""), 1])
                        posBatches <- cbind(posBatches, positives[positives$sample == paste(strain_names[i], "day1PositiveR", j, sep = ""), 9])
                }
        }
}
posMetNames <- unique(positives$metabolite)
colnames(posMat) <- sampleNames
posMat <- cbind.data.frame(posMetNames, as.data.frame(posMat))
colnames(posMat)[1] <- "Compound Name"
posFormulas <- c()
for(i in seq_along(posMetNames)){
        posFormulas<- c(posFormulas, unique(positives$formula[which(positives$metabolite == posMetNames[i])]))
}

posMat <- cbind(posMat, posFormulas)
colnames(posMat)[86] <- "Formula"
colnames(posBatches) <- sampleNames

negatives <- mets[which(mets$mode == "Negative"), ]
negMat <- NULL
negBatches <- NULL
for(i in 1:length(strain_names)){
        for(j in 1:3){
                if(is.null(negMat) == T){
                        negMat <- negatives[negatives$sample == paste(strain_names[i], "day1NegativeR", j, sep = ""), 1]
                        negBatches <- negatives[negatives$sample == paste(strain_names[i], "day1NegativeR", j, sep = ""), 9]
                }else{
                        negMat <- cbind(negMat, negatives[negatives$sample == paste(strain_names[i], "day1NegativeR", j, sep = ""), 1])
                        negBatches <- cbind(negBatches, negatives[negatives$sample == paste(strain_names[i], "day1NegativeR", j, sep = ""), 9])
                }
        }
}
negMetNames <- unique(negatives$metabolite)
colnames(negMat) <- sampleNames
negMat <- cbind.data.frame(negMetNames, as.data.frame(negMat))
colnames(negMat)[1] <- "Compound Name"
negFormulas <- c()
for(i in seq_along(negMetNames)){
        negFormulas<- c(negFormulas, unique(negatives$formula[which(negatives$metabolite == negMetNames[i])]))
}

negMat <- cbind(negMat, negFormulas)
colnames(negMat)[86] <- "Formula"
colnames(negBatches) <- sampleNames

metabolitePeaks_new <- merge(posMat, negMat, all = T)

# Adapt to metabolomics package format. 
newDatMetNames <- metabolitePeaks_new$`Compound Name`
metsNew <- t(metabolitePeaks_new[, !colnames(metabolitePeaks_new) %in% c("Compound Name", "Formula")])
colnames(metsNew) <- newDatMetNames

metsNewAdded <- as.data.frame(metsNew[, dictionary$State[match(colnames(metsNew), 
                                                               dictionary$New.Data.Names)] == "Added"])


# Hexose diphosphate has a minimim that is very low and is always the same, so substitute by a NA.
metsNewAdded$`hexose diphosphate`[metsNewAdded$`hexose diphosphate` == min(metsNewAdded$`hexose diphosphate`)] <- NA


# Merge the old dataset with the added metabolites in new dataset. 
mets_parsed <- cbind.data.frame(metFormat,
                                metsNewAdded)

write.csv(mets_parsed, file = paste0(outDir, "mets_peakRaw_parsed.csv"))


# Create a dataframe of the batch each strain belongs too (all the replicates belonging to same 
# samples were ran in the same batch, so there's equivalence between batch of sample and batch 
# of strain).
batch1 <- unique(mets$mutant[mets$batch == unique(mets$batch)[1]])
batch2 <- unique(mets$mutant[mets$batch == unique(mets$batch)[2]])

batch1 %in% batch2
batch2 %in% batch1

batchDF <- data.frame(strain = unique(mets$mutant),
                      batch = rep(NA, length(unique(mets$mutant))))

for(i in 1:nrow(batchDF)){
        strain <- as.character(batchDF$strain[i])
        if(strain %in% batch1){
                batchDF$batch[i] <- "batch1"
        }else{
                batchDF$batch[i] <- "batch2"
        }
}

write.csv(batchDF, paste0(outDir, "batchDF.csv"))