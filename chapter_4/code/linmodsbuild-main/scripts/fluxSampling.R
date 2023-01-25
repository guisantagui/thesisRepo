if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(caret)) install.packages("caret")
library(caret)
if(!require(dplyr)) install.packages("dplyr")
library(dplyr)

whatMod <- "delsAllSNPs"
medium <- "mi7H9OADCChol"
alphaORA <- 0.05

dataDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/data/"
resDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/results/simulations/"
fbaDir <- sprintf("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/results/simulations/%s/FBA/%sMed/", 
                  whatMod,
                  medium)
rxnDetsFile <- paste0(dataDir, "rxnDets/iEK1011_2.0_rxnDets.csv")
remRxnsFile <- sprintf("%s%s/rxnsRem_%s.txt", resDir, whatMod, whatMod)

oraDir <- sprintf("%s%s/sampling/%s/ora/", resDir, whatMod, medium)
sampDir <- sprintf("%s%s/sampling/%s/", resDir, whatMod, medium)
oplsdaDir <- sprintf("%soplsda/", sampDir)

otherModsDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/results/simsOtherMods/"

# Load data
####################################################################################################

# reaction details
rxnDets <- read.csv(rxnDetsFile, row.names = 1, stringsAsFactors = F)
rxnDets$subsystem <- gsub("Transpor", "Transport", rxnDets$subsystem)
rxnDets$subsystem[rxnDets$subsystem == ""] <- "Any subsystem"

# OPLS-DA data 
load(paste0(oplsdaDir, "allSampFlxs_OPLSDA.RData"))
load(paste0(oplsdaDir, "allSampFlxs_noRXNs_OPLSDA.RData"))

# Reactions altered during model building
remRxns <- read.table(remRxnsFile, 
                      header = F, 
                      stringsAsFactors = F)$V1

# Mann Whitney significative reactions
mwSign <- read.csv(paste0(oraDir, "mwSignPAdjRxns.csv"), row.names = 1)
mwSign$rxnName <- rxnDets$name[match(mwSign$rxn, rxnDets$id)]
mwSign$subsystem <- rxnDets$subsystem[match(mwSign$rxn, rxnDets$id)]
mwSign <- mwSign[, c("rxn", "rxnName", "subsystem", "pValue")]

mwSign <- mwSign[!mwSign$rxn %in% remRxns, ]

# VIP significative reactions
vipSign <- read.csv(paste0(oraDir, "vipSignRxns.csv"), row.names = 1)
vipSign$rxnName <- rxnDets$name[match(vipSign$rxn, rxnDets$id)]
vipSign$subsystem <- rxnDets$subsystem[match(vipSign$rxn, rxnDets$id)]
vipSign <- vipSign[, c("rxn", "rxnName", "subsystem", "vip")]

vipSign <- vipSign[!vipSign$rxn %in% remRxns, ]


# ORA results
oraRes <- read.csv(paste0(oraDir, "mwSignORA.csv"), row.names = 1)
oraSign <- oraRes[oraRes$fish_pValue < alphaORA, ]

oraRes_noRXNs <- read.csv(paste0(oraDir, "mwSignORA_noRXNs.csv"), row.names = 1)
oraSign_noRXNs <- oraRes[oraRes_noRXNs$fish_pValue < alphaORA, ]
rownames(oraSign_noRXNs) <- 1:nrow(oraSign_noRXNs)
oraSign_noRXNs[order(oraSign_noRXNs$fish_pValue), ]

write.csv(oraSign_noRXNs, file = sprintf("%soraSign_noRXNs.csv", oraDir))

# ORA vip
oraRes_vip <- read.csv(paste0(oraDir, "vipSignORA.csv"), row.names = 1)
oraSign_vip <- oraRes_vip[oraRes_vip$fish_pValue < alphaORA, ]

oraRes_vip_noRXNs <- read.csv(paste0(oraDir, "vipSignORA_noRXNs.csv"), row.names = 1)
oraSign_vip_noRXNs <- oraRes[oraRes_vip_noRXNs$fish_pValue < alphaORA, ]
rownames(oraSign_vip_noRXNs) <- 1:nrow(oraSign_vip_noRXNs)
oraSign_vip_noRXNs[order(oraSign_vip_noRXNs$fish_pValue), ]

write.csv(oraSign_vip, file = sprintf("%soraSign_vip.csv", oraDir))
write.csv(oraSign_vip_noRXNs, file = sprintf("%soraSign_vip_noRXNs.csv", oraDir))


# Sampled fluxes 
sampFiles <- list.files(sampDir)[grep("_samp.csv", list.files(sampDir))]
sampList <- lapply(sampFiles, function(x) read.csv(sprintf("%s%s", sampDir, x), row.names = 1))
names(sampList) <- gsub("\\_.*", "", sampFiles)

lins <- names(sampList)
# Create a dataframe with all the fluxes in a single file 
allSampFluxes <- matrix(ncol = length(rxnDets$id), 
                        nrow = length(lins) * 1000, 
                        dimnames = list(make.unique(rep(lins, each = 1000)),
                                        rxnDets$id))

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




# Do OPLS-DA score plots
####################################################################################################

getOPLSDAScore <- function(oplsda){
        scores <- data.frame(sample = rownames(oplsda@scoreMN),
                             t1 = oplsda@scoreMN[, 1],
                             o1 = oplsda@orthoScoreMN[, 1],
                             stringsAsFactors = F)
        lineage <- gsub("\\..*", "", scores$sample)
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
        names(linCols) <- unique(lineage)
        #hostAss <- gsub(".", 
        #                "", 
        #                gsub("[[:digit:]]", 
        #                     "", 
        #                     scores$sample), 
        #                fixed = T)
        hostAss <- gsub("[[:digit:]]", 
                        "", 
                        lineage)
        hostAss <- gsub("A", "Animal associated", hostAss)
        hostAss <- gsub("L", "Human associated", hostAss)
        scores$lineage <- lineage
        scores$host <- factor(hostAss)
        # Plot colored according to host association (human or animal)
        plot <-  ggplot(scores, aes(x=t1, y=o1, color=host)) + 
                geom_point(size = 1.5) +
                geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
                geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) + 
                scale_discrete_manual("Host association",
                                      aesthetics = "colour",
                                      values = c("#4287f5", "#f59642")) +
                #labs(col="Rhamnolipid Production") + 
                theme(axis.text.y   = element_text(size=14),
                      axis.text.x   = element_text(size=14),
                      axis.title.y  = element_text(size=14),
                      axis.title.x  = element_text(size=14),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black", fill=NA, size=1)
                ) +
                xlab(sprintf("Predictive component (%s%%)", round(oplsda@modelDF[1, 1]*100))) +
                ylab("First orthogonal component")
        plot
        # Plot colored according to lineage, shape representing host association
        plotLins <-  ggplot(scores, aes(x=t1, y=o1, color=lineage, shape=host)) + 
                geom_point(size = 1.5) +
                geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
                geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) + 
                scale_discrete_manual("Lineage",
                                      aesthetics = "colour",
                                      values = linCols) +
                scale_shape_manual(name = "Host association", 
                                   values = c(2, 19),
                                   labels = c("Animal associated",
                                              "Human associated")) +
                #labs(col="Rhamnolipid Production") + 
                theme(axis.text.y   = element_text(size=14),
                      axis.text.x   = element_text(size=14),
                      axis.title.y  = element_text(size=14),
                      axis.title.x  = element_text(size=14),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black", fill=NA, size=1)
                ) +
                xlab(sprintf("Predictive component (%s%%)", round(oplsda@modelDF[1, 1]*100))) +
                ylab("First orthogonal component")
        out <- list(byHostAss = plot, 
                    byLineage = plotLins)
        return(out)
        #plot
}

oplsda_scorePlots <- getOPLSDAScore(allSampFlxs_OPLSDA)

oplsda_scorePlots$byHostAss
ggsave(paste0(oplsdaDir, "oplsdaScore_plot.pdf"), height = 4, width = 6)

oplsda_scorePlots$byLineage
ggsave(paste0(oplsdaDir, "oplsdaScore_plot_byLin.pdf"), height = 4, width = 6)

oplsda_scorePlots_noRXNs <- getOPLSDAScore(allSampFlxs_noRXNs_OPLSDA)

oplsda_scorePlots_noRXNs$byHostAss
ggsave(paste0(oplsdaDir, "oplsdaScore_noRXNs_plot.pdf"), height = 4, width = 6)

oplsda_scorePlots_noRXNs$byLin
ggsave(paste0(oplsdaDir, "oplsdaScore_noRXNs_plot_byLin.pdf"), height = 4, width = 6)

# Do loadings plot 
####################################################################################################

getLoadPlot <- function(oplsda, filtSubsyst = NULL, filtUniv = NULL, topFilt = NULL, whatTopFilt = "loads"){
        #if(!require(randomcoloR)) install.packages("randomcoloR") 
        #library(randomcoloR)
        #set.seed(200)
        #subsystCols <- distinctColorPalette(length(unique(rxnDets$subsystem)))
        #names(subsystCols) <- unique(rxnDets$subsystem)
        
        #subsystColsDF <- data.frame(subsystem = names(subsystCols),
        #                            colors = subsystCols)
        subsystColsDF <- read.csv(paste0(oplsdaDir, "subsystCols.csv"), 
                                         row.names = 1, 
                                         stringsAsFactors = F)
        subsystCols <- subsystColsDF$colors
        names(subsystCols) <- subsystColsDF$subsystem
        loadsDF <- data.frame(rxnID = rownames(oplsda@loadingMN),
                              rxnName = rxnDets$name[match(rownames(oplsda@loadingMN),
                                                           rxnDets$id)],
                              subsystem = rxnDets$subsystem[match(rownames(oplsda@loadingMN),
                                                                  rxnDets$id)],
                              loading = oplsda@loadingMN[, 1], 
                              stringsAsFactors = F)
        loadsDF <- loadsDF[order(loadsDF$loading, decreasing = T), ]
        if(!is.null(filtSubsyst)){
                loadsDF <- loadsDF[loadsDF$subsystem %in% filtSubsyst, ]
        }
        if(!is.null(filtUniv)){
                loadsDF <- loadsDF[loadsDF$rxnID %in% filtUniv, ]
        }
        if(!is.null(topFilt)){
                if(whatTopFilt == "loads"){
                        topN <- round(topFilt/2)
                        loadsDF <- loadsDF[c(1:topN, (nrow(loadsDF) - topN):nrow(loadsDF)), ]
                }else if(whatTopFilt == "vip"){
                        topVip <- vipSign[order(vipSign$vip, decreasing = T), ]
                        if(!is.null(filtSubsyst)){
                                topVip <- topVip[topVip$subsystem %in% filtSubsyst, ]
                        }
                        topVip <- head(topVip, topFilt)
                        loadsDF <- loadsDF[loadsDF$rxnID %in% topVip$rxn, ]
                }
                
        }
        loadsDFOrd <- data.frame(matrix(ncol = ncol(loadsDF),
                                        nrow = 0,
                                        dimnames = list(NULL, 
                                                        colnames(loadsDF))))
        loadsSubsysts <- unique(loadsDF$subsystem)
        for(i in seq_along(loadsSubsysts)){
                subLoads <- loadsDF[loadsDF$subsystem == loadsSubsysts[i], ]
                loadsDFOrd <- rbind.data.frame(loadsDFOrd, subLoads)
        }
        loadsDF <- loadsDFOrd
        rownames(loadsDF) <- 1:nrow(loadsDF)
        plot <- ggplot(data = loadsDF,
                       aes(x = loading, y = factor(make.unique(as.character(rxnName)), 
                                                   levels = make.unique(as.character(rxnName))), fill = subsystem)) +
                geom_bar(stat = "identity") +
                labs(x = "loading", y = "reaction") +
                scale_fill_manual(values = subsystCols[match(loadsDF$subsystem, names(subsystCols))]) +
                theme_light()
        if(!is.null(topFilt)){
                plot <- plot +
                        theme(axis.text.y = element_text(size=20),
                              axis.text.x = element_text(size=15),
                              panel.background = element_blank(),
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = "black"),
                              axis.line.y = element_line(colour = "black"),
                              panel.border = element_rect(colour = "black", fill=NA, size=1))
        }
        #return(loadsDF)
        plot
        out <- list(plot = plot,
                    DF = loadsDF,
                    subsystCols = subsystColsDF)
        return(out)
}

signLoadsDF <- getLoadPlot(allSampFlxs_OPLSDA, oraSign$subsystem, filtUniv = mwSign$rxn)
signLoadsDF$plot
ggsave(paste0(oplsdaDir, "oplsdaLoads_plot.pdf"), height = 80, width = 40, limitsize = F)


getLoadPlot(allSampFlxs_OPLSDA, oraSign$subsystem, filtUniv = mwSign$rxn, topFilt = 40)
ggsave(paste0(oplsdaDir, "oplsdaLoads_topLoads_plot.pdf"), height =5, width = 8.5, limitsize = F)

signLoads_vipDF <- getLoadPlot(allSampFlxs_OPLSDA, oraSign_vip$subsystem, filtUniv = vipSign$rxn)
signLoads_vipDF$plot
ggsave(paste0(oplsdaDir, "oplsdaLoads_vipSubsysts_plot.pdf"), height = 30, width = 20, limitsize = F)

signLoads_noRXNs_vipDF <- getLoadPlot(allSampFlxs_OPLSDA, oraSign_vip_noRXNs$subsystem, filtUniv = vipSign$rxn)
signLoads_noRXNs_vipDF$plot
ggsave(paste0(oplsdaDir, "oplsdaLoads_vipSubsysts_noRXNs_plot.pdf"), height = 30, width = 20, limitsize = F)


getLoadPlot(allSampFlxs_OPLSDA, oraSign_vip_noRXNs$subsystem, topFilt = 50, whatTopFilt = "vip")
ggsave(paste0(oplsdaDir, "oplsdaLoads_topVips_noRXNs_plot.pdf"), height =12, width = 13.5)




# Create a dataframe of loads unfiltered and set a copy for writing it to a csv 
# for adding more info and using it as a supplementary table
signLoads_vipDF_allRxns <- getLoadPlot(allSampFlxs_OPLSDA, oraSign_vip$subsystem)
signLoads_vipDF_allRxns_noRXNsORA <- getLoadPlot(allSampFlxs_OPLSDA, oraSign_vip_noRXNs$subsystem)

signLoads_vipDF_allRxns_outDF <- getLoadPlot(allSampFlxs_OPLSDA)$DF
signLoads_vipDF_allRxns_outDF <- signLoads_vipDF_allRxns$DF

signLoads_vipDF_allRxns_outDF$VIP <- allSampFlxs_OPLSDA@vipVn[match(signLoads_vipDF_allRxns_outDF$rxnID,
                                                                    names(allSampFlxs_OPLSDA@vipVn))]


signLoads_vipDF_allRxns_outDF$rxn <- rxnDets$reaction[match(signLoads_vipDF_allRxns_outDF$rxnID,
                                                            rxnDets$id)]

signLoads_vipDF_allRxns_outDF <- signLoads_vipDF_allRxns_outDF[, c("rxnID", "rxnName", "rxn", "subsystem", "loading", "VIP")]

write.csv(signLoads_vipDF_allRxns_outDF, file = paste0(oplsdaDir, "oplsdaRxnInfo.csv"))

#subSystCols <- signLoads_vipDF$subsystCols

#subSystCols[subSystCols == "Cholesterol degradation", 2] <- "#4f2046"
#subSystCols[subSystCols == "Glycolysis/Gluconeogenesis", 2] <- "#ffbdfe"
#subSystCols[subSystCols == "Peptidoglycan Metabolism", 2] <- "#001aff"
#subSystCols[subSystCols == "Cysteine and methionine metabolism", 2] <- "#6f7542"
#subSystCols[subSystCols == "Redox Metabolism", 2] <- "#b01c42"
#subSystCols[subSystCols == "Fatty Acid Metabolism", 2] <- "#543082"
#subSystCols[subSystCols == "Alanine, Aspartate, and Glutamate Metabolism", 2] <- "#e78b2e"
#subSystCols[subSystCols == "Arabinogalactan biosynthesis", 2] <- "#44e62b"
#subSystCols[subSystCols == "Pantothenate and CoA Metabolism", 2] <- "#000000"
#subSystCols[subSystCols == "Lipid metabolism", 2] <- "#ff8787"
#subSystCols[subSystCols == "Tyrosine, Tryptophan, and Phenylalanine Metabolism", 2] <- "#665e80"

#write.csv(subSystCols, file = paste0(oplsdaDir, "subsystCols.csv"))

# Do subsystems correlation plot
####################################################################################################

getSubsystProps <- function(loadsDF){
        subsysts <- unique(loadsDF$subsystem)
        subsystsDF <- data.frame(matrix(ncol = 2, 
                                        nrow = 0,
                                        dimnames = list(NULL,
                                                        c("animal_associated",
                                                          "human_associated"))))
        for(i in seq_along(subsysts)){
                subsyst <- subsysts[i]
                trueVec <- loadsDF$loading[loadsDF$subsystem == subsyst] > 0
                hostAssVec <- rep(NA, length(trueVec))
                hostAssVec[trueVec] <- "human_associated"
                hostAssVec[!trueVec] <- "animal_associated"
                animalNum <- sum(hostAssVec == "animal_associated")
                humanNum <- sum(hostAssVec == "human_associated")
                hostCounts <- c(animalNum, humanNum)
                names(hostCounts) <- c("animal_associated", 
                                       "human_associated")
                hostCounts <- hostCounts/sum(hostCounts)
                subSystDF <- data.frame(animal_associated = hostCounts[1], 
                                        human_associated = hostCounts[2])
                rownames(subSystDF) <- subsyst
                subsystsDF <- rbind.data.frame(subsystsDF, subSystDF)
        }
        return(subsystsDF)
}

opls_subsystProps <- getSubsystProps(signLoadsDF$DF)
opls_subsystProps_vip <- getSubsystProps(signLoads_vipDF_allRxns$DF)
opls_subsystProps_vip_noRXNsORA <- getSubsystProps(signLoads_vipDF_allRxns_noRXNsORA$DF)


plotSubsystProps <- function(subsystProps){
        subsystProps$animal_associated <- subsystProps$animal_associated * -1
        propVec <- c()
        subsysVec <- c()
        for(i in 1:nrow(subsystProps)){
                subsystem <- rownames(subsystProps)[i]
                for(j in 1:2){
                        prop <- subsystProps[i, j]
                        propVec <- c(propVec, prop)
                        subsysVec <- c(subsysVec, subsystem)
                }
        }
        plotDF <- data.frame(subsystem = factor(subsysVec, 
                                                levels = unique(subsysVec)),
                             proportion = propVec)
        plot <- ggplot(plotDF, aes(x = subsystem, y = proportion, fill = subsystem)) +
                geom_bar(
                        stat = "identity", position = position_stack(),
                        color = "white") +
                scale_fill_manual(values = signLoadsDF$subsystCols$colors[match(levels(plotDF$subsystem), 
                                                                                signLoadsDF$subsystCols$subsystem)],
                                  name = "Subsystem") +
                coord_flip() + 
                scale_x_discrete(position = "top") +
                theme(axis.text.y = element_text(size=20),
                      axis.text.x = element_text(size=15),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black", fill=NA, size=1))
                
        plot
}

plotSubsystProps(opls_subsystProps)

ggsave(paste0(oplsdaDir, "oplsdaPropSubsyst_plot.pdf"), width = 15, height = 8.5)

plotSubsystProps(opls_subsystProps_vip)

ggsave(paste0(oplsdaDir, "oplsdaPropSubsyst_vip_plot.pdf"), width = 14.5, height = 6.5)

plotSubsystProps(opls_subsystProps_vip_noRXNsORA)

ggsave(paste0(oplsdaDir, "oplsdaPropSubsyst_vip_noRXNsORA_plot.pdf"), width = 14.5, height = 6.5)

# Explore sample distributions
####################################################################################################

# Add removed reactions distributions to fluxes dataframes (as reactions are shut down distributions will be all 0).
addRemovedRXNs <- function(distrLst){
        distrLst_wRemRXNs <- distrLst
        for(i in seq_along(distrLst)){
                toAdd <- as.character(rxnDets$id[!make.names(rxnDets$id) %in% colnames(distrLst[[i]])])
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

getDensPlot <- function(rxn, color = 'byLineage', xlim = NULL, model = "woSNPs", sampSize = 1000){
        lins <- names(distrList)
        if(model == "woSNPs"){
                flux <- unlist(lapply(distrList, 
                                      function(x) x[, grep(paste("^", make.names(rxn), "$", sep = ""), 
                                                           make.names(colnames(x)))]))
        }#else if(model == "wSNPs"){
         #       flux <- unlist(lapply(distrList_wSnps, 
         #                             function(x) x[, grep(paste("^", make.names(rxn), "$", sep = ""), 
         #                                                  make.names(colnames(x)))]))
        #}
        densDF <- data.frame(lin = rep(lins, each = nrow(distrList[[1]])), 
                             flux = flux,
                             linGroup = gsub("[[:digit:]]", "", rep(lins, each = sampSize)))
        if(rxn == "BIOMASS__2"){
                xAxisTitle <- bquote("Specific growth rate"~(h^-1))
        }else{
                xAxisTitle <- sprintf("%s flux (mmol/(gDW·h))", rxn)
        }
        if(color == "byLineage"){
                p <- ggplot(densDF, aes(x=flux, color=lin)) +
                        geom_density() +
                        scale_color_manual(values = c("#c4bf62",
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
                                                      "#37ff30"),
                                           name = "Lineage") +
                        theme_minimal() + 
                        labs(x = xAxisTitle) +
                        theme(axis.text.y = element_text(size=15),
                              axis.text.x = element_text(size=15),
                              axis.title.x = element_text(size = 16),
                              axis.title.y = element_text(size = 16),
                              panel.border = element_rect(colour = "black", 
                                                          fill=NA, 
                                                          size=1))
        }else if(color == "byLinGroup"){
                p <- ggplot(densDF, aes(x=flux, color=linGroup)) +
                        geom_density() +
                        scale_color_manual(values = c("red",
                                                      "blue"),
                                           name = "Host association",
                                           labels = c("Animal associated",
                                                      "Human associated")) +
                        theme_minimal() +
                        labs(x = xAxisTitle) +
                        theme(axis.text.y = element_text(size=15),
                              axis.text.x = element_text(size=15),
                              axis.title.x = element_text(size = 16),
                              axis.title.y = element_text(size = 16),
                              panel.border = element_rect(colour = "black", 
                                                          fill=NA, 
                                                          size=1))
        }
        histPlot <- ggplot(densDF, aes(x = flux)) +
                geom_histogram(data = subset(densDF, linGroup == "A"), fill = "blue", alpha = 0.3, bins = 100) +
                geom_histogram(data = subset(densDF, linGroup == "L"), fill = "red", alpha = 0.3, bins = 100) +
                theme_minimal() + 
                labs(x = xAxisTitle) +
                theme(axis.text.y = element_text(size=15),
                      axis.text.x = element_text(size=15),
                      axis.title.x = element_text(size = 16),
                      axis.title.y = element_text(size = 16),
                      panel.border = element_rect(colour = "black", 
                                                  fill=NA, 
                                                  size=1))
        if(!is.null(xlim)){
                p <- p + xlim(xlim[1], xlim[2])
                histPlot <- histPlot + xlim(xlim[1], xlim[2])
        }
        out <- list(DF = densDF,
                    densPlot = p,
                    histPlot = histPlot)
        return(out)
}

getRXNsInSubsyst <- function(subsyst, output = "rxnID"){
        if(output == "rxnID"){
                out <- rxnDets$id[rxnDets$subsystem == subsyst]
        }else if(output == "rxnID_names"){
                out <- rxnDets[rxnDets$subsystem == subsyst, c("id", "name", "subsystem")]
        }else if(output == "all"){
                out <- rxnDets[rxnDets$subsystem == subsyst, ]
        }
        return(out)
}

cholDegRXNs <- getRXNsInSubsyst("Cholesterol degradation")

pptGMetRXNs <- getRXNsInSubsyst("Peptidoglycan Metabolism")

biomassDistrs <-  getDensPlot("BIOMASS__2", "byLineage")

biomassDistrs$histPlot
biomassDistrs$densPlot
ggsave(paste0(oplsdaDir, "BM_dens_byLin.pdf"), width = 6, height = 4)

sum(distrList$A1$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$A2$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$A3$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$A4$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$L1$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$L2$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$L3$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$L4$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$L5$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$L6$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$L7$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$L8$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)
sum(distrList$L9$BIOMASS__2 < 0.05)/length(distrList$A1$BIOMASS__2)


# Do a PCA with 10 samples of each lineage with the biomass lower than 
# 0.04 to see if they aggregate according to host association

altSubsystRxns <- unlist(sapply(as.character(oraSign_vip$subsystem), getRXNsInSubsyst))

#altSubsystRxns <- altSubsystRxns[altSubsystRxns %in% vipSign$rxn]

lowGrowthSamps <- data.frame(matrix(nrow = 0, 
                                    #ncol = ncol(distrList$A1),
                                    #ncol = length(as.character(vipSign$rxn)),
                                    ncol = length(altSubsystRxns),
                                    #dimnames = list(NULL, 
                                    #                colnames(distrList$A1)),
                                    #dimnames = list(NULL,
                                    #                as.character(vipSign$rxn)),
                                    dimnames = list(NULL, altSubsystRxns)))
for(i in seq_along(distrList)){
        lin <- names(distrList)[i]
        linSamps <- distrList[[lin]]
        #linSampsFilt <- linSamps[linSamps$BIOMASS__2 < 0.04, ]
        #linSampsFilt <- linSamps[linSamps$BIOMASS__2 < 0.04, as.character(vipSign$rxn)]
        linSampsFilt <- linSamps[linSamps$BIOMASS__2 < 0.025, make.names(altSubsystRxns)]
        set.seed(234)
        linSampsFilt <- linSampsFilt[sample(1:nrow(linSampsFilt), 200), ]
        #linSampsFilt <- linSampsFilt[, colnames(distrList$A1)]
        rownames(linSampsFilt) <- paste(lin, rownames(linSampsFilt), sep = ".")
        lowGrowthSamps <- rbind.data.frame(lowGrowthSamps, linSampsFilt)
}

lowGrowthSampsPCA <- prcomp(lowGrowthSamps)




factoextra::fviz_pca_ind(lowGrowthSampsPCA,
                         col.ind = gsub("\\..*", "", rownames(lowGrowthSamps)),
                         habillage = gsub("\\..*", "", rownames(lowGrowthSamps)),
                         geom = "point") +
        scale_color_manual(name = "Lineages", 
                           labels = gsub("\\..*", "", rownames(lowGrowthSamps)),
                           values = c("#c4bf62",
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
                                      "#37ff30")) +
        scale_shape_manual(name = "Lineages", 
                           values = c(rep(2, 4),
                                      rep(19, 9)),
                           #labels = sapply(linGroup, function(x) if(x == "A") x <- "Animal" else x <- "Human")
                           )




getLowGrowthPCA <- function(bmThrshld = 0.04, nsamps = 10){
        if(!require(factoextra)) install.packages("factoextra")
        library(factoextra)
        lowGrowthSamps <- data.frame(matrix(nrow = 0, 
                                            ncol = ncol(distrList$A1),
                                            dimnames = list(NULL, 
                                                            colnames(distrList$A1))))
        for(i in seq_along(distrList)){
                lin <- names(distrList)[i]
                linSamps <- distrList[[lin]]
                linSampsFilt <- linSamps[linSamps$BIOMASS__2 < bmThrshld, ]
                set.seed(123)
                linSampsFilt <- linSampsFilt[sample(1:nrow(linSampsFilt), nsamps), ]
                linSampsFilt <- linSampsFilt[, colnames(distrList$A1)]
                rownames(linSampsFilt) <- paste(lin, rownames(linSampsFilt), sep = ".")
                lowGrowthSamps <- rbind.data.frame(lowGrowthSamps, linSampsFilt)
        }
        lineages <- gsub("\\..*", "", rownames(lowGrowthSamps))
        linGroup <- gsub("[[:digit:]]", "", lineages)
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
        names(linCols) <- unique(lineages)
        lowGrowthSampsPCA <- prcomp(lowGrowthSamps)
        pcaPlot <- fviz_pca_ind(lowGrowthSampsPCA,
                                col.ind = lineages,
                                habillage = lineages,
                                geom = "point") +
                scale_color_manual(name = "Lineages", 
                                   labels = lineages,
                                   values = linCols) +
                scale_shape_manual(name = "Lineages", 
                                   values = c(rep(2, 4),
                                              rep(19, 9)),
                                   labels = sapply(linGroup, 
                                                   function(x) if(x == "A") x <- "Animal" else x <- "Human"))
        pcaPlot
        out <- list(pcaPlot = pcaPlot,
                    pca = lowGrowthSampsPCA,
                    input = lowGrowthSamps)
        return(out)
}


# Get the minimal sample size that we can perform opls-da
lowGrPCA <- getLowGrowthPCA(bmThrshld = 0.015, nsamps = 5)
#lowGrPCA <- getLowGrowthPCA(bmThrshld = 0.05, nsamps = 5)

lowGrPCA$pcaPlot

lowGr4OPLSDA <- lowGrPCA$input

lowGr4OPLSDA$linGroup <- gsub("[[:digit:]]", "", gsub(".", "", rownames(lowGr4OPLSDA), fixed = T))
        
ropls::opls(lowGr4OPLSDA[, 1:(ncol(lowGr4OPLSDA)-1)],
            lowGr4OPLSDA$linGroup,
            orthoI = 3,
            predI = 1)


mi7H9OADC_chol_mets <- c('EX_glu__L_e',
                         'EX_cu2_e',
                         'EX_btn_e',
                         'EX_pydxn_e',
                         'EX_ca2_e',
                         'EX_mg2_e',
                         'EX_h_e',
                         'EX_k_e',
                         'EX_nh4_e',
                         'EX_h2o_e',
                         'EX_pi_e',
                         'EX_cl_e',
                         'EX_o2_e',
                         'EX_na1_e',
                         'EX_so4_e',
                         'EX_cit_e',
                         'EX_fe3_e',
                         'EX_glyc_e',
                         'EX_glc__D_e',
                         'EX_ocdca_e',
                         'EX_chsterol_e')

View(lowGrPCA$input[, colnames(lowGrPCA$input) %in% mi7H9OADC_chol_mets])


# Compute mean of the sampled fluxes
lowGrMeans <- data.frame(matrix(nrow = 0,
                                ncol = ncol(lowGrPCA$input),
                                dimnames = list(NULL,
                                                colnames(lowGrPCA$input))))
for(l in unique(gsub("\\..*", "", rownames(lowGrPCA$input)))){
        subMat <- lowGrPCA$input[grep(l, rownames(lowGrPCA$input)), ]
        subMatMean <- data.frame(matrix(apply(subMat, 2, mean), 
                                        nrow = 1, ncol = ncol(subMat),
                                        dimnames = list(l,
                                                        colnames(subMat))))
        lowGrMeans <- rbind.data.frame(lowGrMeans,
                                       subMatMean)
}

mi7H9OADC_chol_med <- data.frame(metabolites = c('EX_glu__L_e',
                                                 'EX_cu2_e',
                                                 'EX_btn_e',
                                                 'EX_pydxn_e',
                                                 'EX_ca2_e',
                                                 'EX_mg2_e',
                                                 'EX_h_e',
                                                 'EX_k_e',
                                                 'EX_nh4_e',
                                                 'EX_h2o_e',
                                                 'EX_pi_e',
                                                 'EX_cl_e',
                                                 'EX_o2_e',
                                                 'EX_na1_e',
                                                 'EX_so4_e',
                                                 'EX_cit_e',
                                                 'EX_fe3_e',
                                                 'EX_glyc_e',
                                                 'EX_glc__D_e',
                                                 'EX_ocdca_e',
                                                 'EX_chsterol_e'),
                                 bound = c(1,
                                           1000,
                                           1,
                                           1,
                                           1000,
                                           1000,
                                           1000,
                                           1000,
                                           10,
                                           1000,
                                           1,
                                           1000,
                                           20,
                                           1000,
                                           1000,
                                           1,
                                           5,
                                           1,
                                           1,
                                           1,
                                           1))

# Load FBA results:
fba7h9OADCchol <- read.csv(paste0(fbaDir, 
                                  sprintf("fbaFlxs_%s.csv",
                                          gsub("mi", "",medium))),
                           row.names = 1)

# Compare fluxes in fba and in samples for medium metabolites for each lineage
compFba <- function(mets){
        diffDFLst <- list()
        for(i in seq_along(distrList)){
                lin <- names(distrList)[i]
                distrLin <- distrList[[lin]]
                linMetDF <- data.frame(matrix(nrow = 0, 
                                              ncol = 2,
                                              dimnames = list(NULL,
                                                              c("met", "fluxDiff"))))
                for(j in 1:length(mets)){
                        met <- mets[j]
                        sampVec <- distrLin[, as.character(met)]
                        fbaLinMet <- fba7h9OADCchol[lin, as.character(met)]
                        #minimum <- min(c(sampVec, fbaLinMet))
                        #if(minimum < 0){
                        #        sampVec <- sampVec + abs(minimum)
                        #        fbaLinMet <- fbaLinMet + abs(minimum)
                        #}
                        #metVecDiff <- log2((sampVec + 0.000001)/(fbaLinMet + 0.000001))
                        metVecDiff <- sampVec - fbaLinMet
                        metDiffDF <- data.frame(met = rep(met, length(metVecDiff)),
                                                fluxDiff = metVecDiff)
                        linMetDF <- rbind.data.frame(linMetDF, metDiffDF)
                }
                linMetDF$met <- as.factor(linMetDF$met)
                diffDFLst[[lin]] <- linMetDF
        }
        return(diffDFLst)
}

sampFbaCompMedFlxs <- compFba(as.character(mi7H9OADC_chol_med$metabolites))

boxPlotDiffs <- function(diffFlxs, title = NULL){
        plot <- ggplot(diffFlxs, aes(x=met, y=fluxDiff)) + 
                geom_boxplot() +
                coord_flip() +
                theme(plot.title = element_text(size = 30, face = "bold"),
                      axis.text.y = element_text(size=20),
                      axis.text.x = element_text(size=15),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                labs(x = "Exchange reaction",
                     y = "Flux difference (sampled - FBA)") +
                if(!is.null(title)){
                        print(title)
                        plot <- plot 
                                ggtitle(title)
                }
        plot 
}

boxPlotDiffs(sampFbaCompMedFlxs$A1, "A1")
ggsave(paste0(oplsdaDir, "A1_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$A2, "A2")
ggsave(paste0(oplsdaDir, "A2_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$A3, "A3")
ggsave(paste0(oplsdaDir, "A3_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$A4, "A4")
ggsave(paste0(oplsdaDir, "A4_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$L1, "L1")
ggsave(paste0(oplsdaDir, "L1_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$L2, "L2")
ggsave(paste0(oplsdaDir, "L2_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$L3, "L3")
ggsave(paste0(oplsdaDir, "L3_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$L4, "L4")
ggsave(paste0(oplsdaDir, "L4_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$L5, "L5")
ggsave(paste0(oplsdaDir, "L5_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$L6, "L6")
ggsave(paste0(oplsdaDir, "L6_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$L7, "L7")
ggsave(paste0(oplsdaDir, "L7_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$L8, "L8")
ggsave(paste0(oplsdaDir, "L8_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)
boxPlotDiffs(sampFbaCompMedFlxs$L9, "L9")
ggsave(paste0(oplsdaDir, "L9_sampFBADiff7h9OADCChol.pdf"), 
       width = 4.5, height = 7)

# Try to do a PCA of sampled and fba fluxes to see if they separate.

fbaSampsDF <- fba7h9OADCchol
for(i in seq_along(distrList)){
        lin <- names(distrList)[i]
        distrLin <- distrList[[lin]]
        distrLin <- distrLin[, colnames(fbaSampsDF)]
        rownames(distrLin) <- paste(lin, 
                                    rownames(distrLin), 
                                    sep = ".")
        fbaSampsDF <- rbind.data.frame(fbaSampsDF, distrLin)
}

fbaSampsDF_pca <- prcomp(fbaSampsDF)

fbaOrSamp <- factor(c(rep("fba", 13),
                      rep("samp", (nrow(fbaSampsDF)-13))))

fviz_pca_ind(fbaSampsDF_pca,
             col.ind = fbaOrSamp,
             habillage = fbaOrSamp,
             geom = "point")


fba7h9OADCchol[, colnames(fba7h9OADCchol) %in% mi7H9OADC_chol_med$metabolites]
View(fba7h9OADCchol[4, ])


# Do OPLS-DA separating fba and sampled to identify the fluxes that are different
fbaSampsOPLSDA <- ropls::opls(fbaSampsDF, fbaOrSamp, predI = 1,
                              orthoI = 3)

vipSignFbaSamps <- sort(fbaSampsOPLSDA@vipVn[fbaSampsOPLSDA@vipVn > 1], decreasing = T)

fbaSamps_signRxnsOPLSDA <- data.frame(rxn = names(vipSignFbaSamps),
                                      rxnName = rxnDets$name[match(names(vipSignFbaSamps), 
                                                                   rxnDets$id)],
                                      subsystem = rxnDets$subsystem[match(names(vipSignFbaSamps), 
                                                                          rxnDets$id)],
                                      VIP = vipSignFbaSamps,
                                      loads = fbaSampsOPLSDA@loadingMN[match(names(vipSignFbaSamps), 
                                                                             rownames(fbaSampsOPLSDA@loadingMN))])

View(fbaSamps_signRxnsOPLSDA)


fbaSamps_OPLSDA_loads <- data.frame(rxn = rownames(fbaSampsOPLSDA@loadingMN),
                                    loads = fbaSampsOPLSDA@loadingMN[, 1])

fbaSamps_OPLSDA_loads <- fbaSamps_OPLSDA_loads[order(fbaSamps_OPLSDA_loads$loads, decreasing = T), ]

fbaSamps_OPLSDA_loads$rxnName <- rxnDets$name[match(fbaSamps_OPLSDA_loads$rxn, rxnDets$id)]
fbaSamps_OPLSDA_loads$subsystem <- rxnDets$subsystem[match(fbaSamps_OPLSDA_loads$rxn, rxnDets$id)]
fbaSamps_OPLSDA_loads$VIP <- fbaSampsOPLSDA@vipVn[match(fbaSamps_OPLSDA_loads$rxn, names(fbaSampsOPLSDA@vipVn))]



# Do ORA 
iEK1011_rxnDets <- rxnDets
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

fbaSampORA <- doSubsystORA(fbaSamps_signRxnsOPLSDA$rxn)

fbaSampORA_sign <- fbaSampORA[fbaSampORA$fish_pValue < 0.05, ]

fbaSamps_signRxnsOPLSDA_signORA <- fbaSamps_signRxnsOPLSDA[fbaSamps_signRxnsOPLSDA$subsystem %in% fbaSampORA_sign$subsystem, ]

fbaSamps_signRxnsOPLSDA_signORA[order(fbaSamps_signRxnsOPLSDA_signORA$loads, decreasing = T), ]


View(fbaSamps_signRxnsOPLSDA[order(fbaSamps_signRxnsOPLSDA$loads, decreasing = T), ])

# Get fold change from FBA result to the means to see what's different
fcFBASampDF <- data.frame(matrix(nrow = 0,
                                 ncol = ncol(lowGrMeans),
                                 dimnames = list(NULL,
                                                 colnames(lowGrMeans))))
for(i in 1:nrow(lowGrMeans)){
        lin <- rownames(lowGrMeans)[i]
        linVec <- c()
        for(j in 1:ncol(lowGrMeans)){
                rxn <- colnames(lowGrMeans)[j]
                #fc <- (lowGrMeans[i, j] - fba7h9OADCchol[i, j])/fba7h9OADCchol[i, j]
                fc <- (lowGrMeans[i, j] - fba7h9OADCchol[i, j])#/fba7h9OADCchol[i, j]
                linVec <- c(linVec, fc)
        }
        linDF <- data.frame(matrix(linVec,
                                   nrow = 1,
                                   ncol = ncol(fcFBASampDF),
                                   dimnames = list(lin,
                                                   colnames(lowGrMeans))))
        fcFBASampDF <- rbind.data.frame(fcFBASampDF, linDF)
        
        
}

# Get the ones that are more different
fcFBASampDF_greatest <- fcFBASampDF[, apply(fcFBASampDF, 2, function(x) sum(x > 10 | x < -10) > 10)]

rxnDets[match(colnames(fcFBASampDF_greatest), rxnDets$id), ]

rxnDets[rxnDets$id == 'GLUDxi', ]
rxnDets[rxnDets$id == 'GTHOr', ]







rxnID <- cholDegRXNs[1]
rxnDF <- getDensPlot(rxnID, "byLinGroup")
rxnVec <- unlist(lapply(distrList, 
                        function(x) x[, grep(paste("^", make.names(rxnID), "$", sep = ""), 
                                             make.names(colnames(x)))]))

hist(rxnVec[grep("A", names(rxnVec))], breaks = 100, xlim = c(-0.001, 0.6))
hist(rxnVec[grep("L", names(rxnVec))], breaks = 100, xlim = c(-0.001, 0.6))







max <- which.max(density(rxnVec[grep("L", names(rxnVec))])$y)
xInt <- density(rxnVec[grep("L", names(rxnVec))])$x[max]


cholDistrs



ggplot(rxnDF$DF, aes(x=flux, color=linGroup)) +
        geom_density() +
        scale_color_manual(values = c("red",
                                      "blue")) + 
        geom_vline(xintercept = xInt)




















density(rxnVec[grep("L", names(rxnVec))])$x[which.max(density(rxnVec[grep("L", names(rxnVec))])$y)]

density(rxnVec[grep("A", names(rxnVec))])















getDensPlot(cholDegRXNs[1], "byLinGroup")
getDensPlot(cholDegRXNs[2], "byLinGroup")
getDensPlot(cholDegRXNs[3], "byLinGroup")
getDensPlot(cholDegRXNs[4], "byLinGroup")
getDensPlot(cholDegRXNs[5], "byLinGroup")
getDensPlot(cholDegRXNs[6], "byLinGroup")
getDensPlot(cholDegRXNs[7], "byLinGroup")
getDensPlot(cholDegRXNs[8], "byLinGroup")
getDensPlot(cholDegRXNs[9], "byLinGroup")
getDensPlot(cholDegRXNs[10], "byLinGroup")
getDensPlot(cholDegRXNs[11], "byLinGroup")
getDensPlot(cholDegRXNs[12], "byLinGroup")
getDensPlot(cholDegRXNs[13], "byLinGroup")
getDensPlot(cholDegRXNs[14], "byLinGroup")
getDensPlot(cholDegRXNs[15], "byLinGroup")
getDensPlot(cholDegRXNs[16], "byLinGroup")
getDensPlot(cholDegRXNs[17], "byLinGroup")
getDensPlot(cholDegRXNs[18], "byLinGroup")
getDensPlot(cholDegRXNs[19], "byLinGroup")
getDensPlot(cholDegRXNs[20], "byLinGroup")
getDensPlot(cholDegRXNs[21], "byLinGroup")
getDensPlot(cholDegRXNs[22], "byLinGroup")
getDensPlot(cholDegRXNs[23], "byLinGroup")
getDensPlot(cholDegRXNs[24], "byLinGroup")
getDensPlot(cholDegRXNs[25], "byLinGroup")
getDensPlot(cholDegRXNs[26], "byLinGroup")
getDensPlot(cholDegRXNs[27], "byLinGroup")
getDensPlot(cholDegRXNs[28], "byLinGroup")
getDensPlot(cholDegRXNs[29], "byLinGroup")
getDensPlot(cholDegRXNs[30], "byLinGroup")
getDensPlot(cholDegRXNs[31], "byLinGroup")
getDensPlot(cholDegRXNs[32], "byLinGroup")
getDensPlot(cholDegRXNs[33], "byLinGroup")
getDensPlot(cholDegRXNs[34], "byLinGroup")
getDensPlot(cholDegRXNs[35], "byLinGroup")
getDensPlot(cholDegRXNs[36], "byLinGroup", xlim = c(0, 0.4))
getDensPlot(cholDegRXNs[37], "byLinGroup")
getDensPlot(cholDegRXNs[38], "byLinGroup")
getDensPlot(cholDegRXNs[39], "byLinGroup")
getDensPlot(cholDegRXNs[40], "byLinGroup")

getDensPlot(pptGMetRXNs[1], "byLinGroup", xlim = c(0, 0.000000000001))
getDensPlot(pptGMetRXNs[2], "byLinGroup", xlim = 0.001)
getDensPlot(pptGMetRXNs[3], "byLinGroup", xlim = 0.001)
getDensPlot(pptGMetRXNs[4], "byLinGroup", xlim = 0.001)
getDensPlot(pptGMetRXNs[5], "byLinGroup", xlim = 0.001)
getDensPlot(pptGMetRXNs[6], "byLinGroup", xlim = 5E-5)
getDensPlot(pptGMetRXNs[7], "byLinGroup", xlim = 5E-5)
getDensPlot(pptGMetRXNs[7], "byLinGroup", xlim = 5E-5)
getDensPlot(pptGMetRXNs[7], "byLinGroup", xlim = 5E-5)
getDensPlot(pptGMetRXNs[7], "byLinGroup", xlim = 5E-5)
getDensPlot(pptGMetRXNs[7], "byLinGroup", xlim = 5E-5)


# Explore flux distributions of the human/animal samples that are more separated one from another
####################################################################################################

# This function, given a oplsda object, gets the N samples with more extreme score values in p1
# component (the animal and human samples more different to each other) and returns the sampled fluxes

getTopSamps <- function(oplsda, top){
        lins <- names(distrList)
        sampNames <- rownames(oplsda@scoreMN)
        scores <- oplsda@scoreMN
        animal <- oplsda@scoreMN[grep("A", rownames(oplsda@scoreMN))]
        names(animal) <- sampNames[grep("A", sampNames)]
        human <- oplsda@scoreMN[grep("L", rownames(oplsda@scoreMN))]
        names(human) <- sampNames[grep("L", sampNames)]
        animal <- sort(animal)
        human <- sort(human, decreasing = T)
        animal <- animal[1:top]
        animTopLins <- unique(gsub("\\..*", "", names(animal)))
        human <- human[1:top]
        humnTopLins <- unique(gsub("\\..*", "", names(human)))
        getTopSampDF <- function(topLins, topSamps){
                topSampFluxDF <- data.frame(matrix(ncol = ncol(distrList[[1]]),
                                                   nrow = 0, 
                                                   dimnames = list(NULL,
                                                                   colnames(distrList[[1]]))))
                for(i in seq_along(topLins)){
                        lin <- topLins[i]
                        linSamps <- names(topSamps)[grep(lin, names(topSamps))]
                        # Convert to names within sampled fluxes DF (just the numbers)
                        sampNameDistr <- gsub(".*\\.", "", linSamps)
                        sampNameDistr[sampNameDistr == lin] <- "0"
                        linSampDF <- distrList[[lin]][sampNameDistr, ]
                        rownames(linSampDF) <- linSamps
                        linSampDF <- linSampDF[, colnames(topSampFluxDF)]
                        topSampFluxDF <- rbind.data.frame(topSampFluxDF, linSampDF)
                }
                return(topSampFluxDF)
        }
        topSF_animal <- getTopSampDF(animTopLins, animal)
        topSF_human <- getTopSampDF(humnTopLins, human)
        out <- list(topSF_animal = topSF_animal,
                    topSF_human = topSF_human)
        return(out)
}


oplsdaExtremSamps <- getTopSamps(allSampFlxs_OPLSDA, 1000)

oplsdaExtremSamps$topSF_animal$BIOMASS__2
oplsdaExtremSamps$topSF_human$BIOMASS__2


boxplot(oplsdaExtremSamps$topSF_animal$EX_hia_e,
        oplsdaExtremSamps$topSF_human$EX_hia_e)

boxplot(oplsdaExtremSamps$topSF_animal$MMSAD2,
        oplsdaExtremSamps$topSF_human$MMSAD2)

boxplot(oplsdaExtremSamps$topSF_animal$BIOMASS__2,
        oplsdaExtremSamps$topSF_human$BIOMASS__2)


hist(oplsdaExtremSamps$topSF_animal$BIOMASS__2, breaks = 10, xlim = c(0, 0.2))

hist(oplsdaExtremSamps$topSF_human$BIOMASS__2, breaks = 10, xlim = c(0, 0.2))

max(oplsdaExtremSamps$topSF_animal$BIOMASS__2)

max(oplsdaExtremSamps$topSF_human$BIOMASS__2)


plotTst <- getDensPlot("BIOMASS__2", "byLinGroup")
plotTst$histPlot

max(plotTst$DF[plotTst$DF$linGroup == "L", ]$flux)
max(plotTst$DF[plotTst$DF$lin == "A4", ]$flux)


# Load iSB619 sampled fluxes (model of S. aureus) to see how is biomass distribution

iSB <- read.csv(paste0(otherModsDir, "iSB619_samp.csv"))

hist(iSB$BIOMASS_SA_8a, breaks = 100)

iJN <- read.csv(paste0(otherModsDir, "iJN1463_samp.csv"))

hist(iJN$BIOMASS_KT2440_WT3, breaks = 10)


# Do confusion matrix of OPLS-DSA
predOPLSDA <- ropls::predict(allSampFlxs_OPLSDA, allSampFluxes)

linGroup <- factor(gsub("[[:digit:]]" ,"" ,gsub("\\..*", "", rownames(allSampFluxes))), 
                   levels = c("A", "L"))

confMatDels <- confusionMatrix(predOPLSDA, linGroup)

tableDels <- data.frame(confMatDels$table)

plotTableDels <- tableDels %>%
        mutate(goodbad = ifelse(tableDels$Prediction == tableDels$Reference, "good", "bad")) %>%
        group_by(Reference) %>%
        mutate(prop = Freq/sum(Freq))

p <- ggplot(data = plotTableDels, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = prop)) +
        geom_tile() +
        geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1, size = 20) +
        scale_fill_manual(values = c(good = "green", bad = "red")) +
        theme_bw() +
        xlim(rev(levels(tableDels$Reference)))
p
