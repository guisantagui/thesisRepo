if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(ggpubr)) install.packages("ggpubr")
library(ggpubr)
if(!require(ggpattern))  install_github("coolbutuseless/ggpattern")
library(ggpattern)
if(!require(ropls))  install.packages("ropls")
library(ropls)
if(!require(nnet))  install.packages("nnet")
library(nnet)

pcrResDir <- "C:/Users/Guillem/Documents/PhD/comput/data/qPCR/results/"
dataDir <- "C:/Users/Guillem/Documents/PhD/comput/data/"
citoDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcPCRs/data/"
plotDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcPCRs/plots/"
if(!dir.exists(plotDir)){
        dir.create(plotDir)
}
resDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcPCRs/results/"
if(!dir.exists(resDir)){
        dir.create(resDir)
}

infectTableFile <- paste0(citoDir,
                          "tabla_infecciones_completa_updated2.xlsx")

# Load functions
source("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcPCRs/scripts/pcrAnalFnctns.R")


# Load file with flow citometry data.
##################################################################################################################
citoMap <- data.frame(readxl::read_xlsx(infectTableFile, sheet = 2))

citoMap <- citoMap[citoMap$ID_cito != "NO_CITO", ]

naIdxs <- which(!is.na(citoMap$ID_RNA))
for(i in 1:length(naIdxs)){
        num <- naIdxs[i]
        if(i < length(naIdxs)){
                vecRNA <- citoMap$ID_RNA[num:(((naIdxs)[i + 1])-1)]
                vecMicro <- citoMap$ID_micro[num:(((naIdxs)[i + 1])-1)]
                vecInfDate <- citoMap$Fecha_INF[num:(((naIdxs)[i + 1])-1)]
                vecExp <- citoMap$Experimento[num:(((naIdxs)[i + 1])-1)]
        }else{
                vecRNA <- citoMap$ID_RNA[num:nrow(citoMap)]
                vecMicro <- citoMap$ID_micro[num:nrow(citoMap)]
                vecInfDate <- citoMap$Fecha_INF[num:nrow(citoMap)]
                vecExp <- citoMap$Experimento[num:nrow(citoMap)]
        }
        vecRNA <- rep(vecRNA[1], length(vecRNA))
        vecMicro <- rep(vecMicro[1], length(vecMicro))
        vecInfDate <- rep(vecInfDate[1], length(vecInfDate))
        vecExp <- rep(vecExp[1], length(vecExp))
        citoMap$ID_RNA[num:((num + length(vecRNA)) - 1)] <- vecRNA
        citoMap$ID_micro[num:((num + length(vecMicro)) - 1)] <- vecMicro
        citoMap$Fecha_INF[num:((num + length(vecInfDate)) - 1)] <- vecInfDate
        citoMap$Experimento[num:((num + length(vecExp)) - 1)] <- vecExp
}

#citoMap$ID_RNA <- sapply(citoMap$ID_RNA, getFullRNums)

# Filter to keep just day 6
citoMap_d6 <- citoMap[citoMap$Día_cito == 6, ]

# Filter to keep just day 2
citoMap <- citoMap[citoMap$Día_cito == 2, ]

rnaIDsStrings <- unique(citoMap$ID_RNA)
for(i in seq_along(unique(citoMap$Experimento))){
        exp <- unique(citoMap$Experimento)[i]
        print(unique(citoMap$ID_RNA))
        rnaIDString <- rnaIDsStrings[i]
        
        cases <- paste0(citoMap$Cepas[citoMap$Experimento == exp],
                        citoMap$Célula[citoMap$Experimento == exp])
        uniqCases <- unique(cases)
        if(rnaIDString != "NO_RNA"){
                rnaIDs <- strsplit(rnaIDString, " - ")[[1]]
                idRNASub <- rnaIDs[match(cases, uniqCases)]
        }else{
                idRNASub <- rep(unique(citoMap$ID_RNA)[i], length(cases))
        }
        
        citoMap$ID_RNA[citoMap$Experimento == exp] <- idRNASub
}

# Isolate only rows that have values for the four types of death.
citoMap <- keepCompSamps(citoMap)


citoMap_d6 <- keepCompSamps(citoMap_d6)

# Load oligo info dataframe
oligoInfo <- as.data.frame(readxl::read_xlsx(paste0(dataDir, "RNA_info.xlsx")))
colnames(oligoInfo) <- make.names(colnames(oligoInfo))
oligoInfo <- oligoInfo[oligoInfo$RNA.FRACTION == "Prokaryote" | oligoInfo$RNA.FRACTION == "Total", ]
oligoInfo$INFECTION.DATE <- as.Date(oligoInfo$INFECTION.DATE, origin = "1899-12-30")
oligoInfo$SAMPLING.DATE <- as.Date(oligoInfo$SAMPLING.DATE, origin = "1899-12-30")

# Load and parse results data...

# Results of day 1 RNA extractions

# The first one contains standard curves for pcrs of each gene
pcrRes <- parsePcrRes(paste0(pcrResDir, "2021-04-26.xlsx"))
pcrRes20210506 <- parsePcrRes(paste0(pcrResDir, "2021-05-06.xlsx"))
pcrRes20210510 <- parsePcrRes(paste0(pcrResDir, "2021-05-10.xlsx"))
pcrRes20210511 <- parsePcrRes(paste0(pcrResDir, "2021-05-11.xlsx"))
pcrRes20210512 <- parsePcrRes(paste0(pcrResDir, "2021-05-12.xlsx"))
pcrRes20211029 <- parsePcrRes(paste0(pcrResDir, "2021-10-29.xlsx"))
pcrRes20211209 <- parsePcrRes(paste0(pcrResDir, "2021-12-09.xlsx"))
pcrRes20211213 <- parsePcrRes(paste0(pcrResDir, "2021-12-13.xlsx"))
pcrRes20211221 <- parsePcrRes(paste0(pcrResDir, "2021-12-21.xlsx"))

# Results of day 3 RNA extractions
pcrRes20220216 <- parsePcrRes(paste0(pcrResDir, "2022-02-16.xlsx"))
# pcrRes20220221 <- parsePcrRes(paste0(pcrResDir, "2022-02-21.xlsx")) # 20220221 and 20220302 are the same, but the first one got sprinkles when removing
pcrRes20220302 <- parsePcrRes(paste0(pcrResDir, "2022-03-02.xlsx"))   # a plastic and therefore we repeated it
pcrRes20220222 <- parsePcrRes(paste0(pcrResDir, "2022-02-22.xlsx"))
# Only bovis plate 
pcrRes20220311 <- parsePcrRes(paste0(pcrResDir, "2022-03-11.xlsx"))

# Do scales dataframe, fit linear model to log10 concentrations and obtain efficiency and duplication
# factor
resScales <- data.frame(copiesInPcr = as.numeric(sapply(pcrRes20211221$Sample.Name[c(1:84, 90:95)], 
                                                        function(x) strsplit(x, 
                                                                             " ")[[1]][2])),
                        replicate = paste0("R", sapply(pcrRes20211221$Sample.Name[c(1:84, 90:95)], 
                                                       function(x) strsplit(x, 
                                                                            " ")[[1]][3])),
                        log2CopiesInPcr = log2(as.numeric(sapply(pcrRes20211221$Sample.Name[c(1:84, 90:95)], 
                                                                 function(x) strsplit(x, 
                                                                                      " ")[[1]][2]))),
                        log10CopiesInPcr = log10(as.numeric(sapply(pcrRes20211221$Sample.Name[c(1:84, 90:95)], 
                                                                   function(x) strsplit(x, 
                                                                                        " ")[[1]][2]))),
                        Ct = as.numeric(pcrRes20211221$CT[c(1:84, 90:95)]),
                        primerPair = as.character(sapply(pcrRes20211221$Sample.Name[c(1:84, 90:95)], 
                                                         function(x) strsplit(x, 
                                                                              " ")[[1]][1])))

resScales$replicate <- gsub("NA", "1", resScales$replicate)
resScales$copiesInPcr <- as.numeric(gsub(60, 0, resScales$copiesInPcr))
resScales$log10CopiesInPcr <- log10(resScales$copiesInPcr)
resScales$log2CopiesInPcr <- log2(resScales$copiesInPcr)
resScales$primerPair <- as.character(resScales$primerPair)
resScales$primerPair <- gsub("Mb1684c", "Mb1784c", resScales$primerPair)
resScales$primerPair[resScales$primerPair == "Rv2349"] <- "Rv2349c"



my.formula <- y ~ x


regPlot <- ggplot(resScales, aes(x = log10CopiesInPcr, y = Ct, color = primerPair)) +
        geom_point() +
        geom_smooth(method = "lm", se = F, formula = my.formula) +
        scale_color_manual(values = c("green", "blue", "red", "purple", "orange")) + 
        stat_regline_equation(label.x = c(2, 2, 2, 2, 2), 
                              label.y = c(15, 14.2, 13.4, 12.6, 11.8),
                              aes(label =  paste(..eq.label.., ..rr.label.., ..adj.rr.label.., sep = "~~~~"))) +
        theme_bw(base_size = 14)

regPlot

ggsave(paste0(plotDir, "standCurves.pdf"), height = 6, width = 7)
ggsave(paste0(plotDir, "standCurves.png"), height = 6, width = 7)

# Fit linear models in order to extract slopes
sigALM_4Fact <- lm(data = resScales[resScales$primerPair == "sigA" & !is.na(resScales$Ct), ], 
                       Ct~log10CopiesInPcr)
Rv2349cLM_4Fact <- lm(data = resScales[resScales$primerPair == "Rv2349c" & !is.na(resScales$Ct), ], 
                          Ct~log10CopiesInPcr)
Rv2350cLM_4Fact <- lm(data = resScales[resScales$primerPair == "Rv2350c" & !is.na(resScales$Ct), ], 
                          Ct~log10CopiesInPcr)
Rv2351cLM_4Fact <- lm(data = resScales[resScales$primerPair == "Rv2351c" & !is.na(resScales$Ct), ], 
                          Ct~log10CopiesInPcr)
Mb1784cLM_4Fact <- lm(data = resScales[resScales$primerPair == "Mb1784c" & !is.na(resScales$Ct), ], 
                          Ct~log10CopiesInPcr)

# Get efficiencies and amplification factors from slopes
sigA_eff <- getEff(sigALM_4Fact)
Rv2349c_eff <- getEff(Rv2349cLM_4Fact)
Rv2350c_eff <- getEff(Rv2350cLM_4Fact)
Rv2351c_eff <- getEff(Rv2351cLM_4Fact)
Mb1784c_eff <- getEff(Mb1784cLM_4Fact)

effs <- list(sigA = sigA_eff, 
             Rv2349c = Rv2349c_eff, 
             Rv2350c = Rv2350c_eff, 
             Rv2351c = Rv2351c_eff,
             Mb1784c = Mb1784c_eff)


# Day 1
results20210426 <- parseSamps(pcrRes)
results20210506 <- parseSamps(pcrRes20210506)
results20210510 <- parseSamps(pcrRes20210510)
results20210511 <- parseSamps(pcrRes20210511)
results20210512 <- parseSamps(pcrRes20210512)
results20211029 <- parseSamps(pcrRes20211029)
results20211209 <- parseSamps(pcrRes20211209)
results20211213 <- parseSamps(pcrRes20211213)
results20211221 <- parseSamps(pcrRes20211221)


# Only bovis plate 
results20220311_D1 <- parseSamps(pcrRes20220311, filtD1 = T, filtD3 = F)

# Day 1
results20210426_clean <- remBovis(cleanSamps(results20210426)) # Added rem bovis, to add bovis samples later (from bovis only plate, 20220311)
results20210506_clean <- remBovis(cleanSamps(results20210506)) # Added rem bovis, to add bovis samples later (from bovis only plate, 20220311)
results20210510_clean <- remBovis(cleanSamps(results20210510)) # Added rem bovis, to add bovis samples later (from bovis only plate, 20220311)
results20210511_clean <- remBovis(cleanSamps(results20210511)) # Added rem bovis, to add bovis samples later (from bovis only plate, 20220311)
# Remove R0025 from this experiment, as the other plc genes were measured in the 
# experiment of 04/26 and this sample was removed because of high Ct in sigA
results20210511_clean <- results20210511_clean[results20210511_clean$oligo != "R0025", ]

# Remove bovis in these experiment, as it is from before we located the contamination
results20210512_clean <- remBovis(cleanSamps(results20210512))
results20211029_clean <- remBovis(cleanSamps(results20211029)) # Added rem bovis, to add bovis samples later (from bovis only plate, 20220311)
results20211209_clean <- remBovis(cleanSamps(results20211209)) # Added rem bovis, to add bovis samples later (from bovis only plate, 20220311)
results20211213_clean <- remBovis(cleanSamps(results20211213)) # Added rem bovis, to add bovis samples later (from bovis only plate, 20220311)
results20211221_clean <- remBovis(cleanSamps(results20211221)) # Added rem bovis, to add bovis samples later (from bovis only plate, 20220311)


results20210426_clean$CtNorm <- rep(getSigACPosCts(pcrRes), nrow(results20210426_clean))
results20210506_clean$CtNorm <- rep(getSigACPosCts(pcrRes20210506), nrow(results20210506_clean))
results20210510_clean$CtNorm <- rep(getSigACPosCts(pcrRes20210510), nrow(results20210510_clean))
results20210511_clean$CtNorm <- rep(getSigACPosCts(pcrRes20210511), nrow(results20210511_clean))
results20210512_clean$CtNorm <- rep(getSigACPosCts(pcrRes20210512), nrow(results20210512_clean))
results20211029_clean$CtNorm <- rep(getSigACPosCts(pcrRes20211029), nrow(results20211029_clean))
results20211209_clean$CtNorm <- rep(getSigACPosCts(pcrRes20211209), nrow(results20211209_clean))
results20211213_clean$CtNorm <- rep(getSigACPosCts(pcrRes20211213), nrow(results20211213_clean))
results20211221_clean$CtNorm <- rep(getSigACPosCts(pcrRes20211221), nrow(results20211221_clean))



# Only bovis plate 
results20220311_D1_clean <- cleanSamps(results20220311_D1)
results20220311_D1_clean$CtNorm <- rep(getSigACPosCts(pcrRes20220311), nrow(results20220311_D1_clean))


# Merge all experiments 

results <- rbind.data.frame(#getRelExp(results20210426_clean),  # Remove this experiment because there is no PlcD expression
        getRelExp(results20210506_clean),
        getRelExp(results20210510_clean),
        getRelExp(results20210511_clean),
        getRelExp(results20210512_clean),
        getRelExp(results20211029_clean),
        getRelExp(results20211209_clean),
        getRelExp(results20211213_clean), 
        getRelExp(results20211221_clean),
        getRelExp(results20220311_D1_clean)) # Added bovis from separated plate


# Remove R0043, as has a NA in Rv2349 Ct while Ct is low in the other samples of same case
results <- results[results$oligo != "R0043", ]


table(paste(results$lin, results$host, sep = "_"))/c(1, 1, 4, 4, 4, 4, 4, 4)
# There are 4 replicates in all the cases with the exception of bovis/bomac. 

allRefCts <- rbind.data.frame(#results20210426_clean[results20210426_clean$primerPair == "sigA", ],
        results20210506_clean[results20210506_clean$primerPair == "sigA", ],
        results20210510_clean[results20210510_clean$primerPair == "sigA", ],
        #results20210511[results20210511$primerPair == "sigA", ],
        results20210512_clean[results20210512_clean$primerPair == "sigA", ],
        results20211029_clean[results20211029_clean$primerPair == "sigA", ],
        results20211209_clean[results20211209_clean$primerPair == "sigA", ],
        results20211213_clean[results20211213_clean$primerPair == "sigA", ],
        results20211221_clean[results20211221_clean$primerPair == "sigA", ],
        results20220311_D1_clean[results20220311_D1_clean$primerPair == "sigA", ])


# Remove R0043 from here too.
allRefCts <- allRefCts[allRefCts$oligo != "R0043", ]

allRefCts$lin <- gsub("M. bovis", "Bovis", allRefCts$lin)


bestSamples <- getBestSamps(allRefCts, n = 3)


table(paste(allRefCts[allRefCts$oligo %in% bestSamples, ]$lin, allRefCts[allRefCts$oligo %in% bestSamples, ]$host))


print(sprintf("Percentage of selected samples with Ct above 30: %s%%", 
              as.character(sum(allRefCts[allRefCts$oligo %in% bestSamples, ]$Ct > 30)/length(allRefCts[allRefCts$oligo %in% bestSamples, ]$Ct)*100)))


hist(allRefCts[allRefCts$oligo %in% bestSamples, ]$Ct)

results_inCito <- results[results$oligo %in% unlist(sapply(citoMap$ID_RNA, function(x) strsplit(x, ", "))), ]

results <- results[results$oligo %in% bestSamples, ]


# Substitute NAs in relative expression by zeros
results$relExp[is.na(results$relExp)] <- 0

# Calculate relative expression normalized to the value of sigA in positive control, to remove any variability 
# between plates
results$relExp_plateNorm <- results$relExp/effs$sigA$ampFact^(-results$CtNorm)


write.csv(results, file = paste0(resDir, "bestSamplesAll.csv"))

######################################################################################################################
#                                                                                                                    #
# Expression boxplots and univariate tests                                                                           #
#                                                                                                                    #
######################################################################################################################

#
# Do expression boxplots of each strain in each type of macrphage
######################################################################################################################

plotAll <- doBoxPlots(results, plateNorm = T)

plotAll$plot

ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots_plateNorm.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots_plateNorm.pdf"), height = 7.5, width = 6.75)

#
# Univariate analyses
######################################################################################################################

doShap("plcA", "Chimp", "THP-1")
doShap("plcA", "Chimp", "BOMAC")

doShap("plcB", "Chimp", "THP-1")
doShap("plcB", "Chimp", "BOMAC")

doShap("plcC", "Chimp", "THP-1")
doShap("plcC", "Chimp", "BOMAC")

doShap("plcD", "Chimp", "THP-1")
doShap("plcD", "Chimp", "BOMAC")
# All Normal


doShap("plcA", "L6", "THP-1")
doShap("plcA", "L6", "BOMAC")

doShap("plcB", "L6", "THP-1")
doShap("plcB", "L6", "BOMAC")

doShap("plcC", "L6", "THP-1")
doShap("plcC", "L6", "BOMAC")

doShap("plcD", "L6", "THP-1")
doShap("plcD", "L6", "BOMAC")
# All Normal except plcA in L6


doShap("plcA", "L5", "THP-1")
doShap("plcA", "L5", "BOMAC")

doShap("plcB", "L5", "THP-1")
doShap("plcB", "L5", "BOMAC")

doShap("plcC", "L5", "THP-1")
doShap("plcC", "L5", "BOMAC")

doShap("plcD", "L5", "THP-1")
doShap("plcD", "L5", "BOMAC")
# All Normal

doShap("plcD", "Bovis", "THP-1")
doShap("plcD", "Bovis", "BOMAC")
# All Normal


doMW("plcA", "Chimp")
doMW("plcB", "Chimp")
doMW("plcC", "Chimp")
doMW("plcD", "Chimp")
# Any significatives

doMW("plcA", "L6")
doMW("plcB", "L6")
doMW("plcC", "L6")
doMW("plcD", "L6")
# Any significatives

doMW("plcA", "L5")
doMW("plcB", "L5")
doMW("plcC", "L5")
doMW("plcD", "L5")
# Any significatives

doMW("plcD", "Bovis")


doTTst("plcA", "Chimp")
doTTst("plcB", "Chimp")
doTTst("plcC", "Chimp")
doTTst("plcD", "Chimp")
# Any significatives

doTTst("plcA", "L6")
doTTst("plcB", "L6")
doTTst("plcC", "L6")
doTTst("plcD", "L6")
# Any significatives

doTTst("plcA", "L5")
doTTst("plcB", "L5")
doTTst("plcC", "L5")
doTTst("plcD", "L5")
# Any significatives

doTTst("plcA", "Bovis")
doTTst("plcB", "Bovis")
doTTst("plcC", "Bovis")
doTTst("plcD", "Bovis")

# Statistics with plate normalized values
doShap("plcA", "Chimp", "THP-1", plateNorm = T)
doShap("plcA", "Chimp", "BOMAC", plateNorm = T)

doShap("plcB", "Chimp", "THP-1", plateNorm = T)
doShap("plcB", "Chimp", "BOMAC", plateNorm = T)

doShap("plcC", "Chimp", "THP-1", plateNorm = T)
doShap("plcC", "Chimp", "BOMAC", plateNorm = T)

doShap("plcD", "Chimp", "THP-1", plateNorm = T)
doShap("plcD", "Chimp", "BOMAC", plateNorm = T)
# All Normal except plcA in chimp in Bomac


doShap("plcA", "L6", "THP-1", plateNorm = T)
doShap("plcA", "L6", "BOMAC", plateNorm = T)

doShap("plcB", "L6", "THP-1", plateNorm = T)
doShap("plcB", "L6", "BOMAC", plateNorm = T)

doShap("plcC", "L6", "THP-1", plateNorm = T)
doShap("plcC", "L6", "BOMAC", plateNorm = T)

doShap("plcD", "L6", "THP-1", plateNorm = T)
doShap("plcD", "L6", "BOMAC", plateNorm = T)
# All Normal except plcA in L6 in BOMAC


doShap("plcA", "L5", "THP-1", plateNorm = T)
doShap("plcA", "L5", "BOMAC", plateNorm = T)

doShap("plcB", "L5", "THP-1", plateNorm = T)
doShap("plcB", "L5", "BOMAC", plateNorm = T)

doShap("plcC", "L5", "THP-1", plateNorm = T)
doShap("plcC", "L5", "BOMAC", plateNorm = T)

doShap("plcD", "L5", "THP-1", plateNorm = T)
doShap("plcD", "L5", "BOMAC", plateNorm = T)
# All Normal

doShap("plcD", "Bovis", "THP-1", plateNorm = T)
doShap("plcD", "Bovis", "BOMAC", plateNorm = T)
# All Normal


doMW("plcA", "Chimp", plateNorm = T)
doMW("plcB", "Chimp", plateNorm = T)
doMW("plcC", "Chimp", plateNorm = T)
doMW("plcD", "Chimp", plateNorm = T)
# Any significatives

doMW("plcA", "L6", plateNorm = T)
doMW("plcB", "L6", plateNorm = T)
doMW("plcC", "L6", plateNorm = T)
doMW("plcD", "L6", plateNorm = T)
# Any significatives

doMW("plcA", "L5", plateNorm = T)
doMW("plcB", "L5", plateNorm = T)
doMW("plcC", "L5", plateNorm = T)
doMW("plcD", "L5", plateNorm = T)
# Any significatives

doMW("plcD", "Bovis", plateNorm = T)


doTTst("plcA", "Chimp", plateNorm = T)
doTTst("plcB", "Chimp", plateNorm = T)
doTTst("plcC", "Chimp", plateNorm = T)
doTTst("plcD", "Chimp", plateNorm = T)
# Any significatives

doTTst("plcA", "L6", plateNorm = T)
doTTst("plcB", "L6", plateNorm = T)
doTTst("plcC", "L6", plateNorm = T)
doTTst("plcD", "L6", plateNorm = T)
# Any significatives

doTTst("plcA", "L5", plateNorm = T)
doTTst("plcB", "L5", plateNorm = T)
doTTst("plcC", "L5", plateNorm = T)
doTTst("plcD", "L5", plateNorm = T)
# plcA and plcD are significative in L5

doTTst("plcD", "Bovis", plateNorm = T)

######################################################################################################################
#                                                                                                                    #
# Multivariate tests                                                                                                 #
#                                                                                                                    #
######################################################################################################################

#
# Data parsing for obtaining object appropriate for multivariate analysis
######################################################################################################################

# Create numeric matrixes of expression data, where rows are samples and columns are expression of a plc gene. 
expMat <- getExpMat(results, plateNorm = T)


# Sample the citometry data and add to the matrix

# RNA vs citometry day 1
expCitMat <- getExpCitMat(citoMat = citoMap,
                          expMat = expMat,
                          seed = 555)

# RNA vs citometry day 6
expCitMat_d6 <- getExpCitMat(citoMat = citoMap_d6,
                             expMat = expMat,
                             seed = 777) # In the previous analysis it was 777


# Get random sets of expression/citometry data at day1 vs day2 and day 1 vs day 6
# with the aim of checking if the results of the regressions are prevalent with 
# different random seeds
randExpCitMatLst <- list()
for(i in 1:1000){
        randExpCitMatLst[[i]] <- getExpCitMat(citoMat = citoMap,
                                              expMat = expMat,
                                              seed = i)
}

randExpCitMatLst_d1d6 <- list()
for(i in 1:1000){
        randExpCitMatLst_d1d6[[i]] <- getExpCitMat(citoMat = citoMap_d6,
                                                   expMat = expMat,
                                                   seed = i)
}

# Generate numeric matrix of the day 1 and 6 citometry data for performing PCAs
citoNum_d1 <- createNumMat(citoMap)
citoNum_d6 <- createNumMat(citoMap_d6)

######################################################################################################################
#                                                                                                                    #
# PCAs                                                                                                               #
#                                                                                                                    #
######################################################################################################################

citD1_pca <- prcomp(citoNum_d1, scale. = T)

citD6_pca <- prcomp(citoNum_d6, scale. = T)


expCitMat_d6_pca <- prcomp(expCitMat_d6, scale. = T)

citD1Biplot <- pcBiplot(citD1_pca)
citD1Biplot

ggsave(paste0(plotDir, "citD1PCA_fancy.pdf"), height = 4, width = 5)
ggsave(paste0(plotDir, "citD1PCA_fancy.png"), height = 4, width = 5)

citD6Biplot <- pcBiplot(citD6_pca)
citD6Biplot

pcBiplot(citD6_pca, x="PC2", y = "PC3")

ggsave(paste0(plotDir, "citD6PCA_fancy.pdf"), height = 4, width = 5)
ggsave(paste0(plotDir, "citD6PCA_fancy.png"), height = 4, width = 5)

ggarrange(citD1Biplot, 
          citD6Biplot,
          legend = "bottom",
          common.legend = T)

ggsave(paste0(plotDir, "citD1_D6Biplots.pdf"), height = 3.5, width = 6.3)
ggsave(paste0(plotDir, "citD1_D6Biplots.png"), height = 3.5, width = 6.3)

pcBiplot(expCitMat_d6_pca)

ggsave(paste0(plotDir, "expCitD6PCA_fancy.pdf"), height = 4, width = 5)
ggsave(paste0(plotDir, "expCitD6PCA_fancy.png"), height = 4, width = 5)

######################################################################################################################
#                                                                                                                    #
# Linear models of day 1 post infection expression data Vs day 6 post infection citometry data                       #
#                                                                                                                    #
######################################################################################################################

# Do PCA of citometry at day 6 using only necrosis and viability to obtain the virulence measure
expCitMat_d6_scaled <- apply(expCitMat_d6, 2, function(x) (x - mean(x))/sd(x))
citD6_necrViablPCA <- prcomp(expCitMat_d6_scaled[, colnames(expCitMat_d6_scaled) %in% c("necrotic", "viable")], 
                             scale. = F, 
                             center = F)

citD6_necrViablPCA$rotation

apply(citoNum_d6, 2, function(x) (x - mean(x))/sd(x))[, 3:4] %*% citD6_necrViablPCA$rotation

pcBiplot(citD6_necrViablPCA, labels = T)


citD6_necrViablPCA <- prcomp(expCitMat_d6_scaled[, colnames(expCitMat_d6_scaled) %in% c("necrotic", "viable")], 
                             scale. = F, 
                             center = F)

citD6_necrViablPCA$rotation

apply(citoNum_d6, 2, function(x) (x - mean(x))/sd(x))[, 3:4] %*% citD6_necrViablPCA$rotation


# Get PC1, which will be our measure for virulence
expCitMat_d6_virulence <- (expCitMat_d6_scaled[, colnames(expCitMat_d6_scaled) %in% c("necrotic", "viable")] %*% citD6_necrViablPCA$rotation)[, 1]

expCitMat_d6_scaled <- data.frame(expCitMat_d6_scaled)

expCitMat_d6_scaled$virulence <- expCitMat_d6_virulence

expCitMat_d6_scaled$host_preference <- apply(expCitMat_d6_scaled[, 5:8], 1, function(x) sum(x * (citD6_pca$rotation[c(1:4), 2])))


# Fit linear models to each type of death
lateApFit <- lm(formula("late_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(lateApFit)
viableFit <- lm(formula("viable ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(viableFit)
earlApFit <- lm(formula("early_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(earlApFit)
necrotFit <- lm(formula("necrotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(necrotFit)
hostPrFit_d6 <- lm(formula("host_preference ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(hostPrFit_d6)
virulFit <- lm(formula("virulence ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(virulFit)

car::vif(virulFit)


doMergedAddVarPlots(virulFit)


getResVsFitPlot <- function(lMod, standResids = F, CI = 0.95){
        fitVals <- lMod$fitted.values
        resids <- lMod$residuals
        if(standResids){{
                resids <- resids/sd(resids)
        }}
        bactCols <- data.frame(Bacteria = c("L5", "L6", "Bovis", "Chimp"),
                               legendName = c("L5", "L6", "*M. bovis*", "Chimpanzee Bacillus"),
                               color = c("#871414", "#24ad37", "#f279ce", "#c4bf62"),
                               stringsAsFactors = F)
        plotDF <- data.frame(sample = sapply(names(resids), function(x) strsplit(x, split = "_")[[1]][3]),
                             Bacteria = sapply(names(resids), function(x) strsplit(x, split = "_")[[1]][1]),
                             Host = sapply(names(resids), function(x) strsplit(x, split = "_")[[1]][2]),
                             fitted = fitVals,
                             residuals = resids,
                             stringsAsFactors = F)
        plotDF$Host <- gsub("BOMAC", "BoMac", plotDF$Host)
        plotDF$Host <- gsub("THP-1", "THP1", plotDF$Host)
        plotDF$Host <- factor(plotDF$Host, levels = c("THP1", "BoMac"))
        bactLevls <- c("L5", "L6", "Bovis", "Chimp")[c("L5", "L6", "Bovis", "Chimp") %in% plotDF$Bacteria]
        plotDF$Bacteria <- factor(plotDF$Bacteria, levels = bactLevls)
        plotCols <- bactCols$color[match(levels(plotDF$Bacteria), bactCols$Bacteria)]
        cInt <- c(qnorm((1-CI)/2),
                  qnorm(CI + (1-CI)/2))
        plot <- ggplot(data = plotDF, mapping = aes(x = fitted, 
                                                    y = residuals,
                                                    col = Bacteria,
                                                    shape = Host)) +
                scale_discrete_manual("Bacteria",
                                      aesthetics = "colour",
                                      values = plotCols,
                                      labels = bactCols$legendName[match(levels(plotDF$Bacteria),
                                                                         bactCols$Bacteria)]) +
                geom_point() +
                geom_hline(yintercept = 0, alpha = 0.3) + 
                geom_text_repel(aes(x = fitted, 
                              y = residuals,
                              col = Bacteria,
                              label = sample),
                          data = plotDF[plotDF$residuals > cInt[2] | plotDF$residuals < cInt[1], ]) +
                theme(axis.text.y = element_text(size=10),
                      axis.text.x = element_text(size=10),
                      axis.title.x = element_text(size=10),
                      axis.title.y = element_text(size=10),
                      legend.text = ggtext::element_markdown(size = 8),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black", fill=NA, size=1))
        outliers <- plotDF$sample[plotDF$residuals > cInt[2] | plotDF$residuals < cInt[1]]
        if(length(outliers) > 0){
                out <- list(plot = plot, 
                            outliers = outliers)
                return(out)
        }else{
                return(plot)
        }
}

resFitPlot <- getResVsFitPlot(virulFit, standResids = T)$plot
outliers <- getResVsFitPlot(virulFit, standResids = T)$outliers

resFitPlot

virulFit_noOLs <- lm(formula("virulence ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"),
                     data = expCitMat_d6_scaled[!grepl(paste(outliers, collapse = "|"),
                                                       rownames(expCitMat_d6_scaled)), ])
summary(virulFit_noOLs)

doMergedAddVarPlots(virulFit_noOLs)

getResVsFitPlot(virulFit_noOLs)


summary(lm(formula("virulence ~ Rv2349c + Rv2350c"),
           data = expCitMat_d6_scaled[!grepl(paste(outliers, collapse = "|"),
                                             rownames(expCitMat_d6_scaled)), ]))


doMergedAddVarPlots(lm(formula("virulence ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"),
                       data = expCitMat_d6_scaled[!grepl(paste(outliers, collapse = "|"),
                                                         rownames(expCitMat_d6_scaled)) & !grepl("Bovis", rownames(expCitMat_d6_scaled)), ]))


summary(lm(formula("host_preference ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled[!grepl("Bovis", rownames(expCitMat_d6_scaled)), ]))

doMergedAddVarPlots(lm(formula("host_preference ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled[!grepl("Bovis", rownames(expCitMat_d6_scaled)), ]))

summary(lm(formula("virulence ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled[!grepl("Bovis", rownames(expCitMat_d6_scaled)), ]))


doMergedAddVarPlots(lm(formula("virulence ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), 
                       data = expCitMat_d6_scaled[!grepl("Bovis", rownames(expCitMat_d6_scaled)) & !grepl("R0182", rownames(expCitMat_d6_scaled)), ]))

doMergedAddVarPlots(lm(formula("virulence ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), 
                       data = expCitMat_d6_scaled[!grepl("R0182", rownames(expCitMat_d6_scaled)), ]))


summary(lm(formula("virulence ~ Rv2349c"), data = expCitMat_d6_scaled))

summary(lm(formula("virulence ~ Rv2349c"), data = expCitMat_d6_scaled[!grepl("Bovis", rownames(expCitMat_d6_scaled)), ]))

doMergedAddVarPlots(virulFit)



doRegPlot(lm(formula("virulence ~ Rv2349c"),
             data = expCitMat_d6_scaled[!grepl("Bovis", rownames(expCitMat_d6_scaled)), ]),
          varY = "virulence",
          varX = "Rv2349c")

doRegPlot(lm(formula("virulence ~ Rv2349c"),
             data = expCitMat_d6_scaled),
          varY = "virulence",
          varX = "Rv2349c")

getResVsFitPlot(lm(formula("virulence ~ Rv2349c"),
                   data = expCitMat_d6_scaled))

getResVsFitPlot(lm(formula("virulence ~ Rv2349c"),
                   data = expCitMat_d6_scaled[!grepl("Bovis", rownames(expCitMat_d6_scaled)), ]))


summary(lm(formula("virulence ~ Rv2349c"),
           data = expCitMat_d6_scaled))





getRandRegsFreqs <- function(randMatLst, getVirulence = T, getHostPr = T, alpha = 0.05, indVars = c("Rv2349c", 
                                                                                                    "Rv2350c",
                                                                                                    "Rv2351c",
                                                                                                    "Mb1784c"),
                             noBov = F){
        indFormula <- paste(indVars, collapse = " + ")
        freqMat <- data.frame(matrix(rep(0, 4*(1 + length(indVars))),
                                     nrow = 4,
                                     ncol = 1 + length(indVars),
                                     dimnames = list(c("lateAp",
                                                       "viable",
                                                       "earlAp",
                                                       "necrot"),
                                                     c("modelCount",
                                                       indVars))))
        if(getHostPr){
                freqMat <- rbind.data.frame(freqMat,
                                            data.frame(matrix(rep(0, 1 + length(indVars)),
                                                              ncol = 1 + length(indVars),
                                                              nrow = 1,
                                                              dimnames = list("host_preference",
                                                                              colnames(freqMat)))))
        }
        if(getVirulence){
                freqMat <- rbind.data.frame(freqMat,
                                            data.frame(matrix(rep(0, 1 + length(indVars)),
                                                              ncol = 1 + length(indVars),
                                                              nrow = 1,
                                                              dimnames = list("virule",
                                                                              colnames(freqMat)))))
        }
        for(i in seq_along(randMatLst)){
                randMat <- randMatLst[[i]]
                if(noBov){
                        randMat <- randMat[!grepl("Bovis", rownames(randMat)), ]
                }
                randMat <- as.data.frame(scale(randMat))
                if(getVirulence == T){
                        #randMat$virulence <- getVirulence(randMat$necrotic, randMat$viable, rotAngle = 40)
                        # Obtain virulence as the combination of necrosis and viability that most spreads the data
                        randMat$virulence <- (as.matrix(randMat[, colnames(randMat) %in% c("necrotic", "viable")]) %*% citD6_necrViablPCA$rotation)[, 1]
                        
                        virule <- lm(formula(paste0("virulence ~", indFormula)), data = randMat)
                        viruleSum <- summary(virule)
                        virulePVal <- pf(viruleSum$fstatistic[1],
                                         viruleSum$fstatistic[2],
                                         viruleSum$fstatistic[3],
                                         lower.tail = F)
                        viruleGenePVals <- viruleSum[[4]][, 4][2:(1 + length(indVars))]
                        if(virulePVal <= alpha){
                                freqMat["virule", 1] <- freqMat["virule", 1] + 1
                        }
                        freqMat["virule", 2:(1 + length(indVars))][viruleGenePVals <= alpha] <- freqMat["virule", 2:(1 + length(indVars))][viruleGenePVals <= alpha] + 1
                }
                if(getHostPr == T){
                        # Obtain the projection on the 2nd component of the PCA of the 4 types of death, as it separates 
                        # the infections according to host specificity
                        randCit <- randMat[, colnames(randMat) %in% c("early_apoptotic", "late_apoptotic", "necrotic", "viable")]
                        
                        rotMat <- citD6_pca$rotation
                        rotMat <- as.matrix(rotMat[match(colnames(randCit), rownames(rotMat)), "PC2"])
                        randMat$host_preference <- (as.matrix(randCit) %*% rotMat)[, 1]
                        
                        hostPr <- lm(formula("host_preference ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = randMat)
                        hostPrSum <- summary(hostPr)
                        hostPrPVal <- pf(hostPrSum$fstatistic[1],
                                         hostPrSum$fstatistic[2],
                                         hostPrSum$fstatistic[3],
                                         lower.tail = F)
                        hostPrGenePVals <- hostPrSum[[4]][, 4][2:(1 + length(indVars))]
                        if(hostPrPVal <= alpha){
                                freqMat["host_preference", 1] <- freqMat["host_preference", 1] + 1
                        }
                        freqMat["host_preference", 2:(1 + length(indVars))][hostPrGenePVals <= alpha] <- freqMat["host_preference", 2:(1 + length(indVars))][hostPrGenePVals <= alpha] + 1
                }
                
                lateAp <- lm(formula("late_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = randMat)
                lateApSum <- summary(lateAp)
                lateApPVal <- pf(lateApSum$fstatistic[1],
                                 lateApSum$fstatistic[2],
                                 lateApSum$fstatistic[3],
                                 lower.tail = F)
                lateApGenePVals <- lateApSum[[4]][, 4][2:(1 + length(indVars))]
                if(lateApPVal <= alpha){
                        freqMat["lateAp", 1] <- freqMat["lateAp", 1] + 1
                }
                freqMat["lateAp", 2:(1 + length(indVars))][lateApGenePVals <= alpha] <- freqMat["lateAp", 2:(1 + length(indVars))][lateApGenePVals <= alpha] + 1
                
                viable <- lm(formula("viable ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = randMat)
                viableSum <- summary(viable)
                viablePVal <- pf(viableSum$fstatistic[1],
                                 viableSum$fstatistic[2],
                                 viableSum$fstatistic[3],
                                 lower.tail = F)
                viableGenePVals <- viableSum[[4]][, 4][2:(1 + length(indVars))]
                if(viablePVal <= alpha){
                        freqMat["viable", 1] <- freqMat["viable", 1] + 1
                }
                freqMat["viable", 2:(1 + length(indVars))][viableGenePVals <= alpha] <- freqMat["viable", 2:(1 + length(indVars))][viableGenePVals <= alpha] + 1
                
                earlAp <- lm(formula("early_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = randMat)
                earlApSum <- summary(earlAp)
                earlApPVal <- pf(earlApSum$fstatistic[1],
                                 earlApSum$fstatistic[2],
                                 earlApSum$fstatistic[3],
                                 lower.tail = F)
                earlApGenePVals <- earlApSum[[4]][, 4][2:(1 + length(indVars))]
                if(earlApPVal <= alpha){
                        freqMat["earlAp", 1] <- freqMat["earlAp", 1] + 1
                }
                freqMat["earlAp", 2:(1 + length(indVars))][earlApGenePVals <= alpha] <- freqMat["earlAp", 2:(1 + length(indVars))][earlApGenePVals <= alpha] + 1
                
                necrot <- lm(formula("necrotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = randMat)
                necrotSum <- summary(necrot)
                necrotPVal <- pf(necrotSum$fstatistic[1],
                                 necrotSum$fstatistic[2],
                                 necrotSum$fstatistic[3],
                                 lower.tail = F)
                necrotGenePVals <- necrotSum[[4]][, 4][2:(1 + length(indVars))]
                if(necrotPVal <= alpha){
                        freqMat["necrot", 1] <- freqMat["necrot", 1] + 1
                }
                freqMat["necrot", 2:(1 + length(indVars))][necrotGenePVals <= alpha] <- freqMat["necrot", 2:(1 + length(indVars))][necrotGenePVals <= alpha] + 1
                
                #virule <- lm(formula("virulence ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = randMat)
                #viruleSum <- summary(virule)
                #virulePVal <- pf(viruleSum$fstatistic[1],
                #                 viruleSum$fstatistic[2],
                #                 viruleSum$fstatistic[3],
                #                 lower.tail = F)
                #viruleGenePVals <- viruleSum[[4]][, 4][2:5]
                #if(virulePVal <= alpha){
                #        freqMat[5, 1] <- freqMat[5, 1] + 1
                #}
                #freqMat[5, 2:5][viruleGenePVals <= alpha] <- freqMat[5, 2:5][viruleGenePVals <= alpha] + 1
        }
        return(freqMat)
}



getRandRegsFreqs(randExpCitMatLst_d1d6, alpha = .05)

getRandRegsFreqs(randExpCitMatLst_d1d6, alpha = .05, indVars = "Rv2349c")
getRandRegsFreqs(randExpCitMatLst_d1d6, alpha = .05, indVars = "Rv2350c")
getRandRegsFreqs(randExpCitMatLst_d1d6, alpha = .05, indVars = "Rv2351c")
getRandRegsFreqs(randExpCitMatLst_d1d6, alpha = .05, indVars = "Mb1784c")

getRandRegsFreqs(randExpCitMatLst_d1d6, alpha = .05, noBov = T)

getRandRegsFreqs(randExpCitMatLst_d1d6, alpha = .06, indVars = "Rv2349c", noBov = T)
getRandRegsFreqs(randExpCitMatLst_d1d6, alpha = .05, indVars = "Rv2350c", noBov = T)
getRandRegsFreqs(randExpCitMatLst_d1d6, alpha = .05, indVars = "Rv2351c", noBov = T)
getRandRegsFreqs(randExpCitMatLst_d1d6, alpha = .05, indVars = "Mb1784c", noBov = T)



######################################################################################################################
#                                                                                                                    #
# OPLS models of day 1 post infection expression data Vs day 6 post infection citometry data                         #
#                                                                                                                    #
######################################################################################################################


oplsda_virulExp <- opls(x = expCitMat_d6_scaled[, colnames(expCitMat_d6_scaled) %in% c("Rv2349c",
                                                                                       "Rv2350c",
                                                                                       "Rv2351c",
                                                                                       "Mb1784c")],
                        y = expCitMat_d6_scaled$virulence,
                        orthoI = 2,
                        predI = 1,
                        permI  = 20)

oplsda_virulExp@loadingMN

oplsda_virulExp@vipVn


oplsda_virulExp_noBov <- opls(x = expCitMat_d6_scaled[!grepl("Bovis", rownames(expCitMat_d6_scaled)), 
                                                      colnames(expCitMat_d6_scaled) %in% c("Rv2349c",
                                                                                           "Rv2350c",
                                                                                           "Rv2351c",
                                                                                           "Mb1784c")],
                              y = expCitMat_d6_scaled$virulence[!grepl("Bovis", rownames(expCitMat_d6_scaled))],
                              orthoI = 2,
                              predI = 1,
                              permI  = 20)

oplsda_virulExp@loadingMN

oplsda_virulExp@vipVn

oplsda_virulExp@modelDF
oplsda_virulExp@summaryDF[1, 8]
countsDF[, colnames(countsDF) %in% names(oplsda_virulExp@vipVn[oplsda_virulExp@vipVn > 1])] + 1

countsDF <- data.frame(matrix(rep(0, 5),
                              ncol = 5, 
                              nrow = 1,
                              dimnames = list("counts",
                                              c("model",
                                                "Rv2349c",
                                                "Rv2350c",
                                                "Rv2351c",
                                                "Mb1784c"))))

countsDF$model
for(i in seq_along(randExpCitMatLst_d1d6)){
        print(i)
        randDF <- randExpCitMatLst_d1d6[[i]]
        randDF <- as.data.frame(scale(randDF))
        print(prcomp(randDF[, colnames(randDF) %in% c("necrotic",
                                                      "viable")],
                     scale. = F,
                     center = F)$rotation)
        randDF$virulence <- (as.matrix(randDF[, colnames(randDF) %in% c("necrotic", "viable")]) %*% citD6_necrViablPCA$rotation)[, 1]
        oplsRand <- opls(x = randDF[!grepl("Bovis", rownames(randDF)), 
                                    colnames(randDF) %in% c("Rv2349c",
                                                            "Rv2350c",
                                                            "Rv2351c",
                                                            "Mb1784c")],
                         y = randDF$virulence[!grepl("Bovis", rownames(randDF))],
                         orthoI = 2,
                         predI = 1,
                         permI  = 20)
        pR2Y <- oplsRand@summaryDF[1, 7]
        pQ2Y <- oplsRand@summaryDF[1, 8]
        vips <- oplsRand@vipVn
        if(pR2Y <= 0.05 & pQ2Y <= 0.05){
                countsDF$model <- countsDF$model + 1
        }
        countsDF[, colnames(countsDF) %in% names(vips[vips > 1])] <- countsDF[, colnames(countsDF) %in% names(vips[vips > 1])] + 1
}

