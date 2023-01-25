#########################################################################################################
#                                                                                                       #
# Do statistical analyses with the sampled fluxes.                                                      #
#                                                                                                       #
#########################################################################################################

if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(ropls)) install.packages("ropls")
library(ropls)

# Create parsed, which would accept arguments from terminal and pass to the script
parser <- arg_parser("This script does random forest separating animal and human lineages")

parser <- add_argument(parser = parser,
                       arg = c("--resDir",
                               "--whatMod",
                               "--rxnDetsFile",
                               "--mwAlpha",
                               "--medium",
                               "--whatSampSize"),
                       help = c("Directory of results of simulations",
                                "What model set use in the simulations, the one with deletions and all pot deletereous SNPs (delsAllSNPs) or the one with deletions and stop gain SNPs (delsSGSNPs)",
                                "File with reaction details",
                                "Significance threshold for Mann Whitney U test",
                                "Medium where simulations were carried out",
                                "What sampled files to use according to the sample size choosing"),
                       flag = c(F, F, F, F, F, F),
                       default = list("--resDir" = "/storage/PGO/results/mtb/gsmmSims/",
                                      "--whatMod" = "delsAllSNPs",
                                      "--rxnDetsFile" = "/home/guisana/scripts/modsPipe/data/rxnDets/iEK1011_2.0_rxnDets.csv",
                                      "--mwAlpha" = 0.01,
                                      "--medium" = "mi7H9OADCMed",
                                      "--whatSampSize" = "n1000"))

parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################

resDir <- parsed$resDir
whatMod <- parsed$whatMod
medium <- parsed$medium
whatSampSize <- parsed$whatSampSize

sampDir <- sprintf("%s%s/sampling/%s/%s/", resDir, whatMod, medium, whatSampSize)
sampFiles <- list.files(sampDir)
sampFiles <- sampFiles[grepl("_samp.csv", sampFiles)]
rxnDetsFile <- parsed$rxnDetsFile
mwAlpha <- parsed$mwAlpha

oplsdaOutDir <- paste0(sampDir, "oplsda/")
pcaOutDir <- paste0(sampDir, "pca/")
oraOutDir <- paste0(sampDir, "ora/")

if(!dir.exists(oplsdaOutDir)){
        dir.create(oplsdaOutDir)
}
if(!dir.exists(pcaOutDir)){
        dir.create(pcaOutDir)
}
if(!dir.exists(oraOutDir)){
        dir.create(oraOutDir)
}

# Load and parse the data
#########################################################################################################
# Load reaction details
iEK1011_rxnDets <- read.csv(rxnDetsFile, row.names = 1, stringsAsFactors = F)
iEK1011_rxnDets <- as.data.frame(apply(iEK1011_rxnDets, 2, as.character))
iEK1011_rxnDets$subsystem <- as.character(iEK1011_rxnDets$subsystem)
iEK1011_rxnDets$subsystem[grep("Transpor", iEK1011_rxnDets$subsystem)] <- "Transport"
iEK1011_rxnDets$subsystem[grep("Arabinogalactan bioynthesis", iEK1011_rxnDets$subsystem)] <- "Arabinogalactan biosynthesis"

# Load sampled models files 
sampList <- list()
for(f in sampFiles){
        if(grepl("_samp.csv", f)){
                lin <- gsub("_samp.csv", "", f)
        }else{
                lin <- gsub(".csv", "", f)
        }
        flxDist <- read.csv(paste(sampDir, f, sep = ""), row.names = 1)
        print(sprintf("Loading %s from %s...", f, sampDir))
        # flxDist <- flxDist[sample(1:5000, 1000), ]
        sampList[[lin]] <- flxDist
}

lins <- names(sampList)

# Build a single matrix with all the sampled files 
allSampFluxes <- matrix(ncol = length(iEK1011_rxnDets$id), 
                        nrow = length(lins) * 1000, 
                        dimnames = list(make.unique(rep(lins, each = 1000)),
                                        iEK1011_rxnDets$id))

allSampFluxesFile <- sprintf("%sallSampFluxes.txt", sampDir)
        if(!file.exists(allSampFluxesFile)){
        print("Building all lineages sampled fluxes matrix")
        for(i in seq_along(sampList)){
                lin <- names(sampList)[i]
                SL <- sampList[[i]]
                for(rxn in colnames(allSampFluxes)){
                        if(rxn %in% colnames(SL)){
                                allSampFluxes[grep(lin, rownames(allSampFluxes)), colnames(allSampFluxes) == rxn] <- SL[, colnames(SL) == rxn]
                        }else{
                                allSampFluxes[grep(lin, rownames(allSampFluxes)), colnames(allSampFluxes) == rxn] <- rep(0.0, 1000)
                        }
                }
        }
        write.table(allSampFluxes, file = allSampFluxesFile, quote = F, sep = "\t")
        print(sprintf("%s saved in %s", basename(allSampFluxesFile), dirname(allSampFluxesFile)))
}else{
        print(sprintf("%s already exists in %s. Loading %s...", 
                      basename(allSampFluxesFile), 
                      dirname(allSampFluxesFile), 
                      basename(allSampFluxesFile)))
        allSampFluxes <- read.table(allSampFluxesFile, 
                                    header = T, 
                                    sep = "\t", 
                                    row.names = 1)
        allSampFluxes <- as.matrix(allSampFluxes)
}

remRxnsFile <- sprintf("%srxnsRem_%s.txt", resDir, whatMod)
remRxns <- read.csv(remRxnsFile, sep = "\t", header = F)$V1
print(sprintf("rxnsRem_%s.txt loaded from %s", whatMod, resDir))

# Create a samp fluxes dataset where reactions that have been altered directly during model construction
# are removed.

allSampFluxes_noRXNs <- allSampFluxes[, !colnames(allSampFluxes) %in% remRxns]

# Run OPLS-DA
#########################################################################################################

# Obtain lingroup vector
lineages <- rownames(allSampFluxes)
lineages <- gsub("\\..*", "", lineages)
linGroup <- gsub("[[:digit:]]+", "", lineages)

# With all reactions
allSampFluxes_4OPLSDA <- as.data.frame(allSampFluxes)
allSampFluxes_4OPLSDA$linGroup <- factor(linGroup)

oplsda_allRxns_file <- paste0(oplsdaOutDir, "allSampFlxs_OPLSDA.RData")

if(!file.exists(oplsda_allRxns_file)){
        print("Running OPLS-DA of all the sampled fluxes of each model's solution space...")
        pdf(file = paste(oplsdaOutDir, "allSampFluxes_OPLSDA.pdf"))
        allSampFlxs_OPLSDA <- opls(allSampFluxes_4OPLSDA[, 1:(ncol(allSampFluxes_4OPLSDA) - 1)],
                                   allSampFluxes_4OPLSDA$linGroup, 
                                   predI = 1, 
                                   orthoI = 3,
                                   permI = 100
        )
        dev.off()

        save(allSampFlxs_OPLSDA, file = oplsda_allRxns_file)

        print(paste("allSampFlxs_OPLSDA.RData saved at ", oplsdaOutDir, ".", sep = ""))
        print(paste("allSampFlxs_OPLSDA.pdf saved at ", oplsdaOutDir, ".", sep = ""))
}else{
        print(sprintf("%s already exists in %s", 
                      basename(oplsda_allRxns_file), 
                      dirname(oplsda_allRxns_file)))
        load(oplsda_allRxns_file)
}


# With reactions that have not  been directly altered during model building
allSampFluxes_noRXNs_4OPLSDA <- as.data.frame(allSampFluxes_noRXNs)
allSampFluxes_noRXNs_4OPLSDA$linGroup <- factor(linGroup)

oplsda_noRxns_file <- paste0(oplsdaOutDir, "allSampFlxs_noRXNs_OPLSDA.RData")
                                            
if(!file.exists(oplsda_noRxns_file)){
        print("Running OPLS-DA of all the sampled fluxes of each model's solution space...")
        pdf(file = paste(oplsdaOutDir, "allSampFluxes_noRXNs_OPLSDA.pdf"))
        allSampFlxs_noRXNs_OPLSDA <- opls(allSampFluxes_noRXNs_4OPLSDA[, 1:(ncol(allSampFluxes_noRXNs_4OPLSDA) - 1)],
                                          allSampFluxes_noRXNs_4OPLSDA$linGroup, 
                                          predI = 1, 
                                          orthoI = 3,
                                          permI = 100
        )
        dev.off()
        save(allSampFlxs_noRXNs_OPLSDA, file = oplsda_noRxns_file)

        print(paste("allSampFlxs_noRXNs_OPLSDA.RData saved at ", oplsdaOutDir, ".", sep = ""))
        print(paste("allSampFlxs_noRXNs_OPLSDA.pdf saved at ", oplsdaOutDir, ".", sep = ""))
}else{
        print(sprintf("%s already exists in %s", 
                      basename(oplsda_noRxns_file), 
                      dirname(oplsda_noRxns_file)))
        load(oplsda_noRxns_file)
}
if(file.exists(oplsda_noRxns_file) & file.exists(oplsda_allRxns_file)){
        print("OPLS-DA was already done.")
}else{
        print("OPLS-DAs finished!")
}

#
# Do univariate test (Mann Whitney U-test) for determining what reactions have different distributions 
# between animal and human associated models
##########################################################################################################
# Add removed reactions distributions to fluxes dataframes (as reactions are shut down distributions will be all 0).

print("Running univariate analysis between samples belonging to animal associated lineages vs samples of human associated lineages")
addRemovedRXNs <- function(distrLst){
        distrLst_wRemRXNs <- distrLst
        for(i in seq_along(distrLst)){
                toAdd <- as.character(iEK1011_rxnDets$id[!make.names(iEK1011_rxnDets$id) %in% colnames(distrLst[[i]])])
                if(length(toAdd) > 0){
                        for(j in seq_along(toAdd)){
                                rxn2Add <- toAdd[j]
                                distrLst_wRemRXNs[[i]][, rxn2Add] <- rep(0.0, nrow(distrLst[[i]]))
                        }
                }
        }
        return(distrLst_wRemRXNs)
}

distrList <- addRemovedRXNs(sampList)


doMWTest <- function(rxn, alt = "two.sided"){
        flux <- unlist(lapply(distrList, 
                              function(x) x[, grep(paste("^", make.names(rxn), "$", sep = ""), 
                                                   make.names(colnames(x)))]))
        
        A <- flux[grep("A", names(flux))]
        L <- flux[grep("L", names(flux))]
        mwTest <- wilcox.test(A, L, alternative = alt)
        mwTest_pv <- mwTest$p.value
        names(mwTest_pv) <- rxn
        return(mwTest_pv)
}



MWPVals <- sapply(iEK1011_rxnDets$id, doMWTest)

MWPVals_adj <- p.adjust(MWPVals, method = "BH")

MWPVals_adj_sign <- MWPVals_adj[MWPVals_adj < mwAlpha]

# Remove NAs, as they correspond to reactions where all is zero
MWPVals_adj_sign <- MWPVals_adj_sign[!is.na(MWPVals_adj_sign)]

MW_signDF <- data.frame(rxn = names(MWPVals_adj_sign),
                        pValue = MWPVals_adj_sign)
                        
write.csv(MW_signDF, file = paste0(oraOutDir, "mwSignPAdjRxns.csv"))
print(sprintf("mwSignPAdjRxns.csv saved at %s.", oraOutDir))
                        
                        
#
# Do ORA on the altered reactions to see if there are subsytems that are overrepresented among the 
# reactions with stronger loadings.
##########################################################################################################

doSubsystORA <- function(altRXNs){
        altRXNsDF <- data.frame(rxnID = altRXNs,
                                rxnName = iEK1011_rxnDets$name[match(altRXNs, iEK1011_rxnDets$id)],
                                rxnName = iEK1011_rxnDets$subsystem[match(altRXNs, iEK1011_rxnDets$id)])
        subsystTab <- table(iEK1011_rxnDets$subsystem)
        subsystFishPVals <- c()
        subsystChSqPVals <- c()
        for(i in seq_along(subsystTab)){
                subsyst <- names(subsystTab)[i]
                rxnInSubsystSign <- length(altRXNsDF$rxnID[altRXNsDF$subsystem == subsyst])
                rxnSignAll <- length(altRXNsDF$rxnID)
                rxnInSubSyst <- length(iEK1011_rxnDets$id[iEK1011_rxnDets$subsystem == subsyst])
                rxnAll <- length(iEK1011_rxnDets$id)
                contMat <- matrix(c(rxnInSubSyst, 
                                    rxnAll - rxnInSubSyst,
                                    rxnInSubsystSign,
                                    rxnSignAll - rxnInSubsystSign),
                                  ncol = 2,
                                  nrow = 2,
                                  dimnames = list(c("in_subsyst", "not_in_subsyst"),
                                                  c("rxn_not_interest", "rxn_in_interest")))
                subsystFishPVals <- c(subsystFishPVals, fisher.test(contMat)$p.value)
                subsystChSqPVals <- c(subsystChSqPVals, chisq.test(contMat)$p.value)
        }
        oraDF <- data.frame(subsystem = names(subsystTab),
                            fish_pValue = subsystFishPVals,
                            ChSq_pValue = subsystChSqPVals)
        return(oraDF)
}  

# With Mann Whitney significative reactions
mwSignORA <- doSubsystORA(MW_signDF$rxn)

write.csv(mwSignORA, file = paste0(oraOutDir, "mwSignORA.csv"))
print(sprintf("mwSignORA.csv saved at %s.", oraOutDir))

# Removing reactions that were altered directly during model construction.
mwSignORA_noRXNs <- doSubsystORA(MW_signDF$rxn[!MW_signDF$rxn %in% remRxns])

write.csv(mwSignORA_noRXNs, file = paste0(oraOutDir, "mwSignORA_noRXNs.csv"))
print(sprintf("mwSignORA_noRXNs.csv saved at %s.", oraOutDir))

# With VIP significative reactions from OPLS-DA model
vipSign <- data.frame(rxn = names(allSampFlxs_OPLSDA@vipVn[allSampFlxs_OPLSDA@vipVn >= 1]),
                      vip = allSampFlxs_OPLSDA@vipVn[allSampFlxs_OPLSDA@vipVn >= 1])
                      
write.csv(vipSign, file = paste0(oraOutDir, "vipSignRxns.csv"))
print(sprintf("vipSignRxns.csv saved at %s.", oraOutDir))

vipSignORA <- doSubsystORA(vipSign$rxn)

write.csv(vipSignORA, file = paste0(oraOutDir, "vipSignORA.csv"))
print(sprintf("vipSignORA.csv saved at %s.", oraOutDir))

# Removing reactions that were altered directly during model construction.
vipSignORA_noRXNs <- doSubsystORA(vipSign$rxn[!vipSign$rxn %in% remRxns])

write.csv(vipSignORA_noRXNs, file = paste0(oraOutDir, "vipSignORA_noRXNs.csv"))
print(sprintf("vipSignORA_noRXNs.csv saved at %s.", oraOutDir))

print("Over Representation Analysis finished!")

#
# Run PCAs
##########################################################################################################

# With deletions all reactions
print("Running PCA of all the sampled fluxes of each model's solution space...")

allSampFluxesPCA <- prcomp(allSampFluxes)
save(allSampFluxesPCA, file = paste(pcaOutDir, "allSampFluxesPCA.RData"))

print(paste("allSampFlxsPCA.RData saved at ", pcaOutDir, ".", sep = ""))

# Removing reactions that were directly manipulated
print("Running PCA of all the sampled fluxes of each model's solution space, removing reactions that genes that were removed in any model catalyzed...")

allSampFluxes_noRXNs_PCA <- prcomp(allSampFluxes_noRXNs)
save(allSampFluxes_noRXNs_PCA, file = paste(pcaOutDir, "allSampFluxes_noRXNs_PCA.RData"))

print(paste("allSampFlxs_noRXNs_PCA.RData saved at ", pcaOutDir, ".", sep = ""))

# Do PCA plots:
linCols <- c("#c4bf62",
             "#87e0a0",
             "#6db3d6",
             "#f279ce",
             "#ff30c1",
             "#001aff",
             "#8826b5",
             "#ff0000",
             "#871414",
             "#24ad37",
             "#fbff00",
             "#ff9d00",
             "#37ff30")
             
names(linCols) <- lins

pdf(file = paste(pcaOutDir, "allSampFluxesPCA.pdf", sep = ""))
fviz_pca_ind(allSampFluxesPCA, 
             col.ind = lineages,
             habillage = lineages,
             geom = "point") + 
        scale_color_manual(name = "Lineages", 
                           labels = lineages,
                           values = linCols) +
        scale_shape_manual(name = "Lineages", 
                           values = c(rep(2, 4),
                                      rep(19, 9)),
                           labels = sapply(linGroup, function(x) if(x == "A") x <- "Animal" else x <- "Human"))
dev.off()

print(paste("allSampFluxesPCA.pdf saved at ", pcaOutDir, ".", sep = ""))


pdf(file = paste(pcaOutDir, "allSampFluxes_noRXNs_PCA.pdf", sep = ""))
fviz_pca_ind(allSampFluxes_noRXNs_PCA, 
             col.ind = lineages,
             habillage = lineages,
             geom = "point") + 
        scale_color_manual(name = "Lineages", 
                           labels = lineages,
                           values = linCols) +
        scale_shape_manual(name = "Lineages", 
                           values = c(rep(2, 4),
                                      rep(19, 9)),
                           labels = sapply(gsub("\\..*", "", linGroup), function(x) if(x == "A") x <- "Animal" else x <- "Human"))
dev.off()

print(paste("allSampFluxes_noRXNs_PCA.pdf saved at ", pcaOutDir, ".", sep = ""))

print("Principal Component Analysis finished!")