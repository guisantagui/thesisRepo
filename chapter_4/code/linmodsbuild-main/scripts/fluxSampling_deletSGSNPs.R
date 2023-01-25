if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

whatMod <- "delsSGSNPs"
medium <- "mi7H9OADCChol"

dataDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/data/"
resDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/results/simulations/"
rxnDetsFile <- paste0(dataDir, "rxnDets/iEK1011_2.0_rxnDets.csv")
remRxnsFile <- sprintf("%s%s/rxnsRem_delsSGSNPs.txt", resDir, whatMod)

oraDir <- sprintf("%s%s/sampling/%s/ora/", resDir, whatMod, medium)
oplsdaDir <- sprintf("%s%s/sampling/%s/oplsda/", resDir, whatMod, medium)


# Load data
####################################################################################################

# ORA data 
rxnDets <- read.csv(rxnDetsFile, row.names = 1, stringsAsFactors = F)
rxnDets$subsystem <- gsub("Transpor", "Transport", rxnDets$subsystem)
rxnDets$subsystem[rxnDets$subsystem == ""] <- "Any subsystem"

remRxns <- read.table(remRxnsFile, 
                      header = F, 
                      stringsAsFactors = F)$V1

mwSign <- read.csv(paste0(oraDir, "mwSignPAdjRxns.csv"), row.names = 1)
mwSign$rxnName <- rxnDets$name[match(mwSign$rxn, rxnDets$id)]
mwSign$subsystem <- rxnDets$subsystem[match(mwSign$rxn, rxnDets$id)]
mwSign <- mwSign[, c("rxn", "rxnName", "subsystem", "pValue")]

mwSign <- mwSign[!mwSign$rxn %in% remRxns, ]

dim(mwSign[mwSign$pValue < 0.001, ])

dim(mwSign)

dim(rxnDets)

oraRes <- read.csv(paste0(oraDir, "mwSignORA.csv"), row.names = 1)
oraSign <- oraRes[oraRes$fish_pValue < 0.01, ]

oraRes_noRXNs <- read.csv(paste0(oraDir, "mwSignORA_noRXNs.csv"), row.names = 1)
oraSign_noRXNs <- oraRes[oraRes_noRXNs$fish_pValue < 0.01, ]
rownames(oraSign_noRXNs) <- 1:nrow(oraSign_noRXNs)
oraSign_noRXNs[order(oraSign_noRXNs$fish_pValue), ]

write.csv(oraSign_noRXNs, file = sprintf("%soraSign_noRXNs.csv", oraDir))

dim(oraSign_noRXNs)

# OPLS-DA data 
load(paste0(oplsdaDir, "allSampFlxs_OPLSDA.RData"))
load(paste0(oplsdaDir, "allSampFlxs_noRXNs_OPLSDA.RData"))

getOPLSDAScore <- function(oplsda){
        scores <- data.frame(sample = rownames(oplsda@scoreMN),
                             t1 = oplsda@scoreMN[, 1],
                             o1 = oplsda@orthoScoreMN[, 1],
                             stringsAsFactors = F)
        hostAss <- gsub(".", 
                        "", 
                        gsub("[[:digit:]]", 
                             "", 
                             scores$sample), 
                        fixed = T)
        hostAss <- gsub("A", "Animal associated", hostAss)
        hostAss <- gsub("L", "Human associated", hostAss)
        scores$host <- factor(hostAss)
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
}

getOPLSDAScore(allSampFlxs_OPLSDA)
ggsave(paste0(oplsdaDir, "oplsdaScore_plot.pdf"), height = 4, width = 6)

getOPLSDAScore(allSampFlxs_noRXNs_OPLSDA)
ggsave(paste0(oplsdaDir, "oplsdaScore_noRXNs_plot.pdf"), height = 4, width = 6)


# Do loadings plot 

getLoadPlot <- function(oplsda, filtSubsyst, filtUniv = NULL, topFilt = NULL){
        if(!require(randomcoloR)) install.packages("randomcoloR") 
        library(randomcoloR)
        set.seed(200)
        subsystCols <- distinctColorPalette(length(unique(rxnDets$subsystem)))
        names(subsystCols) <- unique(rxnDets$subsystem)
        loadsDF <- data.frame(rxnID = rownames(oplsda@loadingMN),
                              rxnName = rxnDets$name[match(rownames(oplsda@loadingMN),
                                                           rxnDets$id)],
                              subsystem = rxnDets$subsystem[match(rownames(oplsda@loadingMN),
                                                                  rxnDets$id)],
                              loading = oplsda@loadingMN[, 1], 
                              stringsAsFactors = F)
        loadsDF <- loadsDF[order(loadsDF$loading, decreasing = T), ]
        loadsDF <- loadsDF[loadsDF$subsystem %in% filtSubsyst, ]
        if(!is.null(filtUniv)){
                loadsDF <- loadsDF[loadsDF$rxnID %in% filtUniv, ]
        }
        if(!is.null(topFilt)){
                loadsDF <- loadsDF[c(1:topFilt, (nrow(loadsDF) - topFilt):nrow(loadsDF)), ]
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
        #return(loadsDF)
        plot
        out <- list(plot = plot,
                    DF = loadsDF,
                    subsystCols = subsystCols)
        return(out)
}

signLoadsDF <- getLoadPlot(allSampFlxs_OPLSDA, oraSign$subsystem, filtUniv = mwSign$rxn)
signLoadsDF$plot
ggsave(paste0(oplsdaDir, "oplsdaLoads_plot.pdf"), height = 80, width = 40, limitsize = F)

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
                scale_fill_manual(values = signLoadsDF$subsystCols[match(levels(plotDF$subsystem), 
                                                                         names(signLoadsDF$subsystCols))],
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

signLoadsDF$subsystCols

plotSubsystProps(opls_subsystProps)

ggsave(paste0(oplsdaDir, "oplsdaPropSubsyst_plot.pdf"), width = 15, height = 8.5)


signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[1], ] # Glycolisis, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[2], ] # Transport, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[3], ] # Redox metabolism, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[4], ] # Alternate carbon, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[5], ] # Oxidative phosphorylation, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[6], ] # Cholesterol degradation, higher in humans
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[7], ] # Membrane metabolism, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[8], ] # Extracellular exchange, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[9], ] # Glyoxylate metabolism, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[10], ] # Glycerophospholipid metabolism, some positive, some negative, but towards human
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[11], ] # Purine and Pyrimidine Biosynthesis, some positive, some negative, but towards human
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[12], ] # TCA cycle, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[13], ] # Pyruvate metabolism, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[14], ] # PPP, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[15], ] # Fatty acid metabolism, some positive, some negative, but towards animals
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[16], ] # Polyprenyl metabolism, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[17], ] # Cofactor and Prosthetic Group Biosynthesis, mostly positive, few negative but stronger
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[18], ] # Folate metabolism, some positive some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[19], ] # Mycolic acid pathway, mostly positive, but 3 negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[20], ] # Alanine, Aspartate, and Glutamate Metabolism, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[21], ] # Mycolic acid biosynthesis, all positive
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[22], ] # Pantothenate and CoA Metabolism, all positive
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[23], ] # Riboflavin Metabolism, Mostly positive
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[24], ] # Arabinogalactan biosynthesis, positive all except one
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[25], ] # Peptidoglycan Metabolism, all positive
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[26], ] # Sulfur metabolism, one negative two positives
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[27], ] # Lipid metabolism, mostly positive
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[28], ] # Nucleotide Salvage Pathway, some positive some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[29], ] # Intracellular demand, some positive some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[30], ] # Valine, Leucine, and Isoleucine Metabolism, some positive some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[31], ] # Cysteine and methionine metabolism, some positive some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[32], ] # Other Amino Acid Metabolism, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[33], ] # Arginine and Proline Metabolism, some positive, some negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[34], ] # Tyrosine, Tryptophan, and Phenylalanine Metabolism, all negative
signLoadsDF$DF[signLoadsDF$DF$subsystem == signSubsysts[35], ] # Histidine Metabolism, all negative




mwSign$rxn


getDensPlot <- function(rxn, color = 'byLineage', xlim = NULL, model = "woSNPs"){
        if(model == "woSNPs"){
                flux <- unlist(lapply(distrList, 
                                      function(x) x[, grep(paste("^", make.names(rxn), "$", sep = ""), 
                                                           make.names(colnames(x)))]))
        }else if(model == "wSNPs"){
                flux <- unlist(lapply(distrList_wSnps, 
                                      function(x) x[, grep(paste("^", make.names(rxn), "$", sep = ""), 
                                                           make.names(colnames(x)))]))
        }
        densDF <- data.frame(lin = rep(lins, each = nrow(A1_fluxDist)), 
                             flux = flux,
                             linGroup = gsub("[[:digit:]]", "", rep(lins, each = 5000)))
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
                                                      "#37ff30"))
        }else if(color == "byLinGroup"){
                p <- ggplot(densDF, aes(x=flux, color=linGroup)) +
                        geom_density() +
                        scale_color_manual(values = c("red",
                                                      "blue"))
        }
        if(!is.null(xlim)){
                p <- p + xlim(xlim[1], xlim[2])
        }
        p
}
