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

# Curves of cholesterol/mbt/mcr7 genes
pcrRes20220728 <- parsePcrRes(paste0(pcrResDir, "2022-07-28_stands.xlsx"))

# Results of cholesterol/mbt/mcr7 genes
pcrRes20220729_1 <- parsePcrRes(paste0(pcrResDir, "2022-07-29_1.xlsx"))
pcrRes20220729_2 <- parsePcrRes(paste0(pcrResDir, "2022-07-29_2.xlsx"))

# Do scales dataframe, fit linear model to log10 concentrations and obtain efficiency and duplication
# factor
resScales <- data.frame(copiesInPcr = as.numeric(sapply(pcrRes$Sample.Name[1:40], 
                                                        function(x) strsplit(x, 
                                                                             " ")[[1]][2])),
                        log2CopiesInPcr = log2(as.numeric(sapply(pcrRes$Sample.Name[1:40], 
                                                                 function(x) strsplit(x, 
                                                                                      " ")[[1]][2]))),
                        log10CopiesInPcr = log10(as.numeric(sapply(pcrRes$Sample.Name[1:40], 
                                                                   function(x) strsplit(x, 
                                                                                        " ")[[1]][2]))),
                        Ct = as.numeric(pcrRes$CT[1:40]),
                        primerPair = as.character(sapply(pcrRes$Sample.Name[1:40], 
                                                         function(x) strsplit(x, 
                                                                              " ")[[1]][1])))

resScales_new <- data.frame(copiesInPcr = as.numeric(sapply(pcrRes20211221$Sample.Name[c(1:84, 90:95)], 
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

resScalesNuGenes <- data.frame(copiesInPcr = as.numeric(sapply(pcrRes20220728$Sample.Name[c(1:70, 73:82, 85:94)], 
                                                               function(x) strsplit(x, 
                                                                                    " ")[[1]][2])),
                               replicate = paste0("R", sapply(pcrRes20220728$Sample.Name[c(1:70, 73:82, 85:94)], 
                                                              function(x) strsplit(x, 
                                                                                   " ")[[1]][3])),
                               log2CopiesInPcr = log2(as.numeric(sapply(pcrRes20220728$Sample.Name[c(1:70, 73:82, 85:94)], 
                                                                        function(x) strsplit(x, 
                                                                                             " ")[[1]][2]))),
                               log10CopiesInPcr = log10(as.numeric(sapply(pcrRes20220728$Sample.Name[c(1:70, 73:82, 85:94)], 
                                                                          function(x) strsplit(x, 
                                                                                               " ")[[1]][2]))),
                               Ct = as.numeric(pcrRes20220728$CT[c(1:70, 73:82, 85:94)]),
                               primerPair = as.character(sapply(pcrRes20220728$Sample.Name[c(1:70, 73:82, 85:94)], 
                                                                function(x) strsplit(x, 
                                                                                     " ")[[1]][1])))

resScales_new$replicate <- gsub("NA", "1", resScales_new$replicate)
resScales_new$copiesInPcr <- as.numeric(gsub(60, 0, resScales_new$copiesInPcr))
resScales_new$log10CopiesInPcr <- log10(resScales_new$copiesInPcr)
resScales_new$log2CopiesInPcr <- log2(resScales_new$copiesInPcr)
resScales_new$primerPair <- as.character(resScales_new$primerPair)
resScales_new$primerPair <- gsub("Mb1684c", "Mb1784c", resScales_new$primerPair)
resScales_new$primerPair[resScales_new$primerPair == "Rv2349"] <- "Rv2349c"


resScales_new[is.na(resScales_new$Ct), ]
rownames(resScales) <- pcrRes$Sample.Name[1:40]

resScales$primerPair[grep("Rv2349", resScales$primerPair)] <- "Rv2349c"


resScalesNuGenes$replicate <- gsub("NA", "1", resScalesNuGenes$replicate)


my.formula <- y ~ x

regPlot <- ggplot(resScales[resScales$copiesInPcr %in% c(3e+02:3e+06), ], aes(x = log10CopiesInPcr, y = Ct, color = primerPair)) +
        geom_point() +
        geom_smooth(method = "lm", se = F, formula = my.formula) +
        scale_color_manual(values = c("green", "blue", "red", "purple")) + 
        stat_regline_equation(label.x = c(2, 2, 2, 2), 
                              label.y = c(15, 14.2, 13.4, 12.6),
                              aes(label =  paste(..eq.label.., ..rr.label.., ..adj.rr.label.., sep = "~~~~")))

regPlot

regPlot_new <- ggplot(resScales_new, aes(x = log10CopiesInPcr, y = Ct, color = primerPair)) +
        geom_point() +
        geom_smooth(method = "lm", se = F, formula = my.formula) +
        scale_color_manual(values = c("green", "blue", "red", "purple", "orange")) + 
        stat_regline_equation(label.x = c(2, 2, 2, 2, 2), 
                              label.y = c(15, 14.2, 13.4, 12.6, 11.8),
                              aes(label =  paste(..eq.label.., ..rr.label.., ..adj.rr.label.., sep = "~~~~"))) +
        theme_bw(base_size = 14)

regPlot_new

ggsave(paste0(plotDir, "standCurves.pdf"), height = 6, width = 7)
ggsave(paste0(plotDir, "standCurves.png"), height = 6, width = 7)

regPlot_nuGenes <- ggplot(resScalesNuGenes, aes(x = log10CopiesInPcr, y = Ct, color = primerPair)) +
        geom_point() +
        geom_smooth(method = "lm", se = F, formula = my.formula) +
        scale_color_manual(values = c("black", "brown", "pink", "magenta", "grey", "orange")) + 
        stat_regline_equation(label.x = c(2, 2, 2, 2, 2, 2), 
                              label.y = c(15, 14.2, 13.4, 12.6, 11.8, 11),
                              aes(label =  paste(..eq.label.., ..rr.label.., ..adj.rr.label.., sep = "~~~~"))) +
        theme_bw(base_size = 14)

regPlot_nuGenes

ggsave(paste0(plotDir, "standCurves_nuGenes.pdf"), height = 6, width = 7)
ggsave(paste0(plotDir, "standCurves_nuGenes.png"), height = 6, width = 7)


ggplot(resScalesNuGenes[resScalesNuGenes$primerPair == "sigA", ], aes(x = log10CopiesInPcr, y = Ct)) +
        geom_point() +
        geom_smooth(method = "lm", se = F, formula = my.formula) +
        stat_regline_equation(label.x = 2, 
                              label.y = 15,
                              aes(label =  paste(..eq.label.., ..rr.label.., ..adj.rr.label.., sep = "~~~~"))) +
        theme_bw(base_size = 14)


# Fit linear models in order to extract slopes
sigALM_4Fact <- lm(data = resScales[1:5, ], Ct~log10CopiesInPcr)
Rv2349cLM_4Fact <- lm(data = resScales[11:15, ], Ct~log10CopiesInPcr)
Rv2350cLM_4Fact <- lm(data = resScales[21:25, ], Ct~log10CopiesInPcr)
Rv2351cLM_4Fact <- lm(data = resScales[31:35, ], Ct~log10CopiesInPcr)

sigALM_4Fact_new <- lm(data = resScales_new[resScales_new$primerPair == "sigA" & !is.na(resScales_new$Ct), ], 
                       Ct~log10CopiesInPcr)
Rv2349cLM_4Fact_new <- lm(data = resScales_new[resScales_new$primerPair == "Rv2349c" & !is.na(resScales_new$Ct), ], 
                          Ct~log10CopiesInPcr)
Rv2350cLM_4Fact_new <- lm(data = resScales_new[resScales_new$primerPair == "Rv2350c" & !is.na(resScales_new$Ct), ], 
                          Ct~log10CopiesInPcr)
Rv2351cLM_4Fact_new <- lm(data = resScales_new[resScales_new$primerPair == "Rv2351c" & !is.na(resScales_new$Ct), ], 
                          Ct~log10CopiesInPcr)
Mb1784cLM_4Fact_new <- lm(data = resScales_new[resScales_new$primerPair == "Mb1784c" & !is.na(resScales_new$Ct), ], 
                          Ct~log10CopiesInPcr)

sigALM_4Fact_new2 <- lm(data = resScalesNuGenes[resScalesNuGenes$primerPair == "sigA" & !is.na(resScalesNuGenes$Ct), ], 
                        Ct~log10CopiesInPcr)
Rv1106LM_4Fact <- lm(data = resScalesNuGenes[resScalesNuGenes$primerPair == "Rv1106" & !is.na(resScalesNuGenes$Ct), ], 
                     Ct~log10CopiesInPcr)
Rv2379cLM_4Fact <- lm(data = resScalesNuGenes[resScalesNuGenes$primerPair == "Rv2379c" & !is.na(resScalesNuGenes$Ct), ], 
                      Ct~log10CopiesInPcr)
Rv2383cLM_4Fact <- lm(data = resScalesNuGenes[resScalesNuGenes$primerPair == "Rv2383c" & !is.na(resScalesNuGenes$Ct), ], 
                      Ct~log10CopiesInPcr)
Rv3545cLM_4Fact <- lm(data = resScalesNuGenes[resScalesNuGenes$primerPair == "Rv3545c" & !is.na(resScalesNuGenes$Ct), ], 
                      Ct~log10CopiesInPcr)
mcr7LM_4Fact <- lm(data = resScalesNuGenes[resScalesNuGenes$primerPair == "mcr7" & !is.na(resScalesNuGenes$Ct), ], 
                   Ct~log10CopiesInPcr)

# Get efficiencies and amplification factors from slopes
sigA_eff <- getEff(sigALM_4Fact)
Rv2349c_eff <- getEff(Rv2349cLM_4Fact)
Rv2350c_eff <- getEff(Rv2350cLM_4Fact)
Rv2351c_eff <- getEff(Rv2351cLM_4Fact)

sigA_eff_new <- getEff(sigALM_4Fact_new)
Rv2349c_eff_new <- getEff(Rv2349cLM_4Fact_new)
Rv2350c_eff_new <- getEff(Rv2350cLM_4Fact_new)
Rv2351c_eff_new <- getEff(Rv2351cLM_4Fact_new)
Mb1784c_eff_new <- getEff(Mb1784cLM_4Fact_new)

sigA_eff_new2 <- getEff(sigALM_4Fact_new2)
Rv1106_eff <- getEff(Rv1106LM_4Fact)
Rv2379c_eff <- getEff(Rv2379cLM_4Fact)
Rv2383c_eff <- getEff(Rv2383cLM_4Fact)
Rv3545c_eff <- getEff(Rv3545cLM_4Fact)
mcr7_eff <- getEff(mcr7LM_4Fact)


effs <- list(sigA = sigA_eff, 
             Rv2349c = Rv2349c_eff, 
             Rv2350c = Rv2350c_eff, 
             Rv2351c = Rv2351c_eff)

# Do another one with the new efficiencies to not need to change the other functions
effs <- list(sigA = sigA_eff_new, 
             Rv2349c = Rv2349c_eff_new, 
             Rv2350c = Rv2350c_eff_new, 
             Rv2351c = Rv2351c_eff_new,
             Mb1784c = Mb1784c_eff_new)

lapply(list("ampFact", "efficiency"), function(x) sigA_eff_new2[[x]]/sigA_eff_new2[[x]] * sigA_eff_new[[x]])



# Adjust efficiencies between the two sets of parameters, as the experiments were performed
# on different dates, to have the same efficiency and amplification factor in sigA gene
effsNuGenes <- list(sigA = as.list(sapply(1:2, function(x) sigA_eff_new2[[x]]/sigA_eff_new2[[x]] * sigA_eff_new[[x]])), 
                    Rv1106 = as.list(sapply(1:2, function(x) Rv1106_eff[[x]]/sigA_eff_new2[[x]] * sigA_eff_new[[x]])), 
                    Rv2379c = as.list(sapply(1:2, function(x) Rv2379c_eff[[x]]/sigA_eff_new2[[x]] * sigA_eff_new[[x]])), 
                    Rv2383c = as.list(sapply(1:2, function(x) Rv2383c_eff[[x]]/sigA_eff_new2[[x]] * sigA_eff_new[[x]])),
                    Rv3545c = as.list(sapply(1:2, function(x) Rv3545c_eff[[x]]/sigA_eff_new2[[x]] * sigA_eff_new[[x]])),
                    mcr7 = as.list(sapply(1:2, function(x) mcr7_eff[[x]]/sigA_eff_new2[[x]] * sigA_eff_new[[x]])))

for(n in names(effsNuGenes)){
        names(effsNuGenes[[n]]) <- c("ampFact", "efficiency")
}

effs$Rv1106 <- effsNuGenes$Rv1106
effs$Rv2379c <- effsNuGenes$Rv2379c
effs$Rv2383c <- effsNuGenes$Rv2383c
effs$Rv3545c <- effsNuGenes$Rv3545c
effs$mcr7 <- effsNuGenes$mcr7

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

# Day 1 new genes
results20220729_1 <- parseSamps(pcrRes20220729_1)
results20220729_2 <- parseSamps(pcrRes20220729_2)
# Day 3
results20220216 <- parseSamps(pcrRes20220216, filtD1 = F, filtD3 = T)
# results20220221 <- parseSamps(pcrRes20220221, filtD1 = F, filtD3 = T)
results20220302 <- parseSamps(pcrRes20220302, filtD1 = F, filtD3 = T)  # Repeated 20220221 because of possible contamination
results20220222 <- parseSamps(pcrRes20220222, filtD1 = F, filtD3 = T)

# Only bovis plate 
results20220311_D3 <- parseSamps(pcrRes20220311, filtD1 = F, filtD3 = T)
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

# Day 1, new genes
results20220729_1_clean <- cleanSamps(results20220729_1)
results20220729_2_clean <- cleanSamps(results20220729_2)


results20210426_clean$CtNorm <- rep(getSigACPosCts(pcrRes), nrow(results20210426_clean))
results20210506_clean$CtNorm <- rep(getSigACPosCts(pcrRes20210506), nrow(results20210506_clean))
results20210510_clean$CtNorm <- rep(getSigACPosCts(pcrRes20210510), nrow(results20210510_clean))
results20210511_clean$CtNorm <- rep(getSigACPosCts(pcrRes20210511), nrow(results20210511_clean))
results20210512_clean$CtNorm <- rep(getSigACPosCts(pcrRes20210512), nrow(results20210512_clean))
results20211029_clean$CtNorm <- rep(getSigACPosCts(pcrRes20211029), nrow(results20211029_clean))
results20211209_clean$CtNorm <- rep(getSigACPosCts(pcrRes20211209), nrow(results20211209_clean))
results20211213_clean$CtNorm <- rep(getSigACPosCts(pcrRes20211213), nrow(results20211213_clean))
results20211221_clean$CtNorm <- rep(getSigACPosCts(pcrRes20211221), nrow(results20211221_clean))

results20220729_1_clean$CtNorm <- rep(getSigACPosCts(pcrRes20220729_1), nrow(results20220729_1_clean))
results20220729_2_clean$CtNorm <- rep(getSigACPosCts(pcrRes20220729_2), nrow(results20220729_2_clean))


# Day 3
results20220216_clean <- remBovis(cleanSamps(results20220216)) # Added rem bovis, to add bovis samples later (from bovis only plate, 20220311)
#results20220221_clean <- remBovis(cleanSamps(results20220221)) # Added rem bovis, to add bovis samples later (from bovis only plate, 20220311)
results20220302_clean <- remBovis(cleanSamps(results20220302)) # Repeated 20220221 because of possible contamination
results20220222_clean <- remBovis(cleanSamps(results20220222)) # Added rem bovis, to add bovis samples later (from bovis only plate, 20220311)

# Day 3
results20220216_clean$CtNorm <- rep(getSigACPosCts(pcrRes20220216), nrow(results20220216_clean))
# results20220221_clean$CtNorm <- rep(getSigACPosCts(pcrRes20220221), nrow(results20220221_clean))
results20220302_clean$CtNorm <- rep(getSigACPosCts(pcrRes20220302), nrow(results20220302_clean)) # Repeated 20220221 because of possible contamination
results20220222_clean$CtNorm <- rep(getSigACPosCts(pcrRes20220222), nrow(results20220222_clean))

# Only bovis plate 
results20220311_D1_clean <- cleanSamps(results20220311_D1)
results20220311_D3_clean <- cleanSamps(results20220311_D3)

results20220311_D1_clean$CtNorm <- rep(getSigACPosCts(pcrRes20220311), nrow(results20220311_D1_clean))
results20220311_D3_clean$CtNorm <- rep(getSigACPosCts(pcrRes20220311), nrow(results20220311_D3_clean))

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

resultsD3 <- rbind.data.frame(getRelExp(results20220216_clean),
                              #getRelExp(results20220221_clean),
                              getRelExp(results20220302_clean), # Repeated 20220221 because of possible contamination
                              getRelExp(results20220222_clean),
                              getRelExp(results20220311_D3_clean)) # Added bovis from separated plate

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

allRefCtsD3 <- rbind.data.frame(results20220216_clean[results20220216_clean$primerPair == "sigA", ],
                                #results20220221_clean[results20220221_clean$primerPair == "sigA", ],
                                results20220302_clean[results20220302_clean$primerPair == "sigA", ],
                                results20220222_clean[results20220222_clean$primerPair == "sigA", ],
                                results20220311_D3_clean[results20220311_D3_clean$primerPair == "sigA", ])

# Remove R0043 from here too.
allRefCts <- allRefCts[allRefCts$oligo != "R0043", ]

allRefCts$lin <- gsub("M. bovis", "Bovis", allRefCts$lin)
allRefCtsD3$lin <- gsub("M. bovis", "Bovis", allRefCtsD3$lin)


bestSamples <- getBestSamps(allRefCts, n = 3)

bestSamplesD3 <- getBestSamps(allRefCtsD3, n = 3)

table(paste(allRefCts[allRefCts$oligo %in% bestSamples, ]$lin, allRefCts[allRefCts$oligo %in% bestSamples, ]$host))

table(paste(allRefCtsD3[allRefCtsD3$oligo %in% bestSamplesD3, ]$lin, allRefCtsD3[allRefCtsD3$oligo %in% bestSamplesD3, ]$host))

print(sprintf("Percentage of selected samples with Ct above 30: %s%%", 
              as.character(sum(allRefCts[allRefCts$oligo %in% bestSamples, ]$Ct > 30)/length(allRefCts[allRefCts$oligo %in% bestSamples, ]$Ct)*100)))

print(sprintf("Percentage of selected samples with Ct above 30 in day 3: %s%%", 
              as.character(sum(allRefCtsD3[allRefCtsD3$oligo %in% bestSamplesD3, ]$Ct > 30)/length(allRefCtsD3[allRefCtsD3$oligo %in% bestSamplesD3, ]$Ct)*100)))

hist(allRefCts[allRefCts$oligo %in% bestSamples, ]$Ct)

results_inCito <- results[results$oligo %in% unlist(sapply(citoMap$ID_RNA, function(x) strsplit(x, ", "))), ]

results <- results[results$oligo %in% bestSamples, ]

resultsD3 <- resultsD3[resultsD3$oligo %in% bestSamplesD3, ]

# Add new genes to global results

rbind.data.frame(getRelExp(results20220729_1_clean),
                 getRelExp(results20220729_2_clean))


results <- rbind.data.frame(results,
                            getRelExp(results20220729_1_clean),  # Add new genes experiments
                            getRelExp(results20220729_2_clean))

# Substitute NAs in relative expression by zeros
results$relExp[is.na(results$relExp)] <- 0

resultsD3$relExp[is.na(resultsD3$relExp)] <- 0

results$relExp_plateNorm <- results$relExp/effs$sigA$ampFact^(-results$CtNorm)

resultsD3$relExp_plateNorm <- resultsD3$relExp/effs$sigA$ampFact^(-resultsD3$CtNorm)




write.csv(results, file = paste0(resDir, "bestSamplesAll.csv"))

write.csv(resultsD3, file = paste0(resDir, "bestSamplesAllD3.csv"))

######################################################################################################################
#                                                                                                                    #
# Expression boxplots and univariate tests                                                                           #
#                                                                                                                    #
######################################################################################################################

#
# Do expression boxplots of each strain in each type of macrphage
######################################################################################################################

plotAll <- doBoxPlots(results)

plotAll_D3 <- doBoxPlots(resultsD3)

plotAll$plot

ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots.pdf"), height = 7.5, width = 6.75)

doBoxPlots(results, plateNorm = T)$plot

ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots_plateNorm.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots_plateNorm.pdf"), height = 7.5, width = 6.75)

plotAll_D3$plot

ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots_D3.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots_D3.pdf"), height = 7.5, width = 6.75)

doBoxPlots(resultsD3, plateNorm = T)$plot

ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots_D3_plateNorm.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots_D3_plateNorm.pdf"), height = 7.5, width = 6.75)


# New genes boxplots
doExpBoxPlot(results[results$primerPair == "Rv1106", ])
doExpBoxPlot(results[results$primerPair == "Rv2379c", ])
doExpBoxPlot(results[results$primerPair == "Rv2383c", ])
doExpBoxPlot(results[results$primerPair == "Rv3545c", ])
doExpBoxPlot(results[results$primerPair == "mcr7", ])


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

# Try with just new
results_new <- rbind.data.frame(getRelExp(results20211029),
                                getRelExp(results20211209),
                                getRelExp(results20211213),
                                getRelExp(results20211221))

allRefCts_new <- rbind.data.frame(results20211029[results20211029$primerPair == "sigA", ],
                                  results20211209[results20211209$primerPair == "sigA", ],
                                  results20211213[results20211213$primerPair == "sigA", ],
                                  results20211221[results20211221$primerPair == "sigA", ])


table(paste(results_new$lin, results_new$host, sep = "_"))

# Take 3 samples of each case
bestSamples_new <- getBestSamps(allRefCts_new, n = 3)

results_new <- results_new[results_new$oligo %in% bestSamples_new, ]

results_new$relExp[is.na(results_new$relExp)] <- 0

write.csv(results, file = paste0(resDir, "bestSamplesNew.csv"))

plotNew <- doBoxPlots(results_new)

plotNew$plot

ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots_new.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "relExpPlcABCD_sepBoxplots_new.pdf"), height = 7.5, width = 6.75)

# Statistics 
doShap("plcA", "Chimp", "THP-1", which = "new")
doShap("plcA", "Chimp", "BOMAC", which = "new")

doShap("plcB", "Chimp", "THP-1", which = "new")
doShap("plcB", "Chimp", "BOMAC", which = "new")

doShap("plcC", "Chimp", "THP-1", which = "new")
doShap("plcC", "Chimp", "BOMAC", which = "new")

doShap("plcD", "Chimp", "THP-1", which = "new")
doShap("plcD", "Chimp", "BOMAC", which = "new")
# All Normal


doShap("plcA", "L6", "THP-1", which = "new")
doShap("plcA", "L6", "BOMAC", which = "new")

doShap("plcB", "L6", "THP-1", which = "new")
doShap("plcB", "L6", "BOMAC", which = "new")

doShap("plcC", "L6", "THP-1", which = "new")
doShap("plcC", "L6", "BOMAC", which = "new")

doShap("plcD", "L6", "THP-1", which = "new")
doShap("plcD", "L6", "BOMAC", which = "new")
# All Normal except plcA in L6


doShap("plcA", "L5", "THP-1", which = "new")
doShap("plcA", "L5", "BOMAC", which = "new")

doShap("plcB", "L5", "THP-1", which = "new")
doShap("plcB", "L5", "BOMAC", which = "new")

doShap("plcC", "L5", "THP-1", which = "new")
doShap("plcC", "L5", "BOMAC", which = "new")

doShap("plcD", "L5", "THP-1", which = "new")
doShap("plcD", "L5", "BOMAC", which = "new")
# All Normal

doShap("plcD", "Bovis", "THP-1", which = "new")
doShap("plcD", "Bovis", "BOMAC", which = "new")
# All Normal


doMW("plcA", "Chimp", which = "new")
doMW("plcB", "Chimp", which = "new")
doMW("plcC", "Chimp", which = "new")
doMW("plcD", "Chimp", which = "new")
# Any significatives

doMW("plcA", "L6", which = "new")
doMW("plcB", "L6", which = "new")
doMW("plcC", "L6", which = "new")
doMW("plcD", "L6", which = "new")
# Any significatives

doMW("plcA", "L5", which = "new")
doMW("plcB", "L5", which = "new")
doMW("plcC", "L5", which = "new")
doMW("plcD", "L5", which = "new")
# Any significatives

doMW("plcD", "Bovis", which = "new")


doTTst("plcA", "Chimp", which = "new")
doTTst("plcB", "Chimp", which = "new")
doTTst("plcC", "Chimp", which = "new")
doTTst("plcD", "Chimp", which = "new")
# Any significatives

doTTst("plcA", "L6", which = "new")
doTTst("plcB", "L6", which = "new")
doTTst("plcC", "L6", which = "new")
doTTst("plcD", "L6", which = "new")
# Any significatives

doTTst("plcA", "L5", which = "new")
doTTst("plcB", "L5", which = "new")
doTTst("plcC", "L5", which = "new")
doTTst("plcD", "L5", which = "new")
# Any significatives

doTTst("plcD", "Bovis", which = "new")

# Create numeric matrixes of expression data, where rows are samples and columns are expression of a plc gene. 
expMat <- getExpMat(results, plateNorm = T)
expMatD3 <- getExpMat(resultsD3, plateNorm = T)

# Sample the citometry data and add to the matrix

expCitMat <- getExpCitMat(citoMat = citoMap,
                          expMat = expMat,
                          seed = 555)

expCitMat_d6 <- getExpCitMat(citoMat = citoMap_d6,
                             expMat = expMat)

expCitMat_d3_d6 <- getExpCitMat(citoMat = citoMap_d6,
                                expMat = expMatD3)

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


# Try to see if the sum of the expression of all the plc genes has differences.
totalPlcD1 <- apply(expMat, 1, sum)

totalPlcD1 <- data.frame(sample = sapply(names(totalPlcD1), function(x) strsplit(x, "_")[[1]][3]),
                         lin = sapply(names(totalPlcD1), function(x) strsplit(x, "_")[[1]][1]),
                         host = sapply(names(totalPlcD1), function(x) strsplit(x, "_")[[1]][2]),
                         plcTotExp = totalPlcD1)

ggplot(totalPlcD1, aes(x = lin, fill = lin, y = plcTotExp)) +
        geom_boxplot_pattern(pattern_color = "black",
                             pattern_fill = "black",
                             pattern_spacing = 0.015,
                             aes(pattern = host)) + 
        scale_pattern_manual(values = c("none", "stripe"),
                             labels = c("BOMAC", "THP-1")) +
        labs(title = "total phospholipase C",
             x = "Strain",
             y = "Expression (relative to *sigA*)") +
        #scale_fill_manual(values = c("#c4bf62", "#24ad37", "#871414"), 
        #                  labels = c("Chimp", "L6", "L5")) +
        scale_fill_manual(values = c("#f279ce", "#c4bf62", "#24ad37", "#871414"), 
                          labels = c("Bovis", "Chimp", "L6", "L5")) +
        theme(title = ggtext::element_markdown(),
              axis.title.y = ggtext::element_markdown(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.grid.major = element_line(colour = "#d4d4d4"))

shapiro.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "Bovis" & totalPlcD1$host == "THP-1"])
shapiro.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "Bovis" & totalPlcD1$host == "BOMAC"])

shapiro.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "Chimp" & totalPlcD1$host == "THP-1"])
shapiro.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "Chimp" & totalPlcD1$host == "BOMAC"])

shapiro.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "L6" & totalPlcD1$host == "THP-1"])
shapiro.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "L6" & totalPlcD1$host == "BOMAC"])

shapiro.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "L5" & totalPlcD1$host == "THP-1"])
shapiro.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "L5" & totalPlcD1$host == "BOMAC"])


t.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "Bovis" & totalPlcD1$host == "THP-1"],
       totalPlcD1$plcTotExp[totalPlcD1$lin == "Bovis" & totalPlcD1$host == "BOMAC"])

t.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "Chimp" & totalPlcD1$host == "THP-1"],
       totalPlcD1$plcTotExp[totalPlcD1$lin == "Chimp" & totalPlcD1$host == "BOMAC"])

t.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "L6" & totalPlcD1$host == "THP-1"],
       totalPlcD1$plcTotExp[totalPlcD1$lin == "L6" & totalPlcD1$host == "BOMAC"])

t.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "L5" & totalPlcD1$host == "THP-1"],
       totalPlcD1$plcTotExp[totalPlcD1$lin == "L5" & totalPlcD1$host == "BOMAC"])


wilcox.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "Bovis" & totalPlcD1$host == "THP-1"],
            totalPlcD1$plcTotExp[totalPlcD1$lin == "Bovis" & totalPlcD1$host == "BOMAC"])

wilcox.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "Chimp" & totalPlcD1$host == "THP-1"],
            totalPlcD1$plcTotExp[totalPlcD1$lin == "Chimp" & totalPlcD1$host == "BOMAC"])

wilcox.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "L6" & totalPlcD1$host == "THP-1"],
            totalPlcD1$plcTotExp[totalPlcD1$lin == "L6" & totalPlcD1$host == "BOMAC"])

wilcox.test(totalPlcD1$plcTotExp[totalPlcD1$lin == "L5" & totalPlcD1$host == "THP-1"],
            totalPlcD1$plcTotExp[totalPlcD1$lin == "L5" & totalPlcD1$host == "BOMAC"])

# Nothing significative

######################################################################################################################
#                                                                                                                    #
# PCAs                                                                                                               #
#                                                                                                                    #
######################################################################################################################

citMat_pca <- prcomp(expCitMat[, !colnames(expCitMat) %in% c("Rv2349c", "Rv2350c", "Rv2351c", "Mb1784c")],
                     scale. = T)

expCitMat_pca <- prcomp(expCitMat, scale. = T)

expCitMat_d6_pca <- prcomp(expCitMat_d6, scale. = T)

expCitMat_d3_d6_pca <- prcomp(expCitMat_d3_d6, scale. = T)

citD1_pca <- prcomp(citoNum_d1, scale. = T)

citD6_pca <- prcomp(citoNum_d6, scale. = T)

screePlot(citMat_pca)
screePlot(expCitMat_pca, nComps = 8)
screePlot(expCitMat_d6_pca, nComps = 8)
screePlot(expCitMat_d3_d6_pca, nComps = 8)
screePlot(citD1_pca)
screePlot(citD6_pca)

pcBiplot(expCitMat_pca)

ggsave(paste0(plotDir, "expCitPCA_fancy.pdf"), height = 4, width = 5)
ggsave(paste0(plotDir, "expCitPCA_fancy.png"), height = 4, width = 5)

pcBiplot(expCitMat_d6_pca)

ggsave(paste0(plotDir, "expCitD6PCA_fancy.pdf"), height = 4, width = 5)
ggsave(paste0(plotDir, "expCitD6PCA_fancy.png"), height = 4, width = 5)

pcBiplot(expCitMat_d3_d6_pca)

ggsave(paste0(plotDir, "expCitD3D6PCA_fancy.pdf"), height = 4, width = 5)
ggsave(paste0(plotDir, "expCitD3D6PCA_fancy.png"), height = 4, width = 5)

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

ggarrange(citD1Biplot, 
          citD6Biplot,
          legend = "right",
          common.legend = T)

ggsave(paste0(plotDir, "citD1_D6Biplots_rightLeg.pdf"), height = 3.5, width = 8.5)
ggsave(paste0(plotDir, "citD1_D6Biplots_rightLeg.png"), height = 3.5, width = 8.5)


pcBiplot(citMat_pca)

ggsave(paste0(plotDir, "citPCA_fancy.pdf"), height = 4, width = 5)
ggsave(paste0(plotDir, "citPCA_fancy.png"), height = 4, width = 5)

expMat_pca <- prcomp(expCitMat[, colnames(expCitMat) %in% c("Rv2349c", "Rv2350c", "Rv2351c", "Mb1784c")],
                     scale. = T)

pcBiplot(expMat_pca)

ggsave(paste0(plotDir, "expPCA_fancy.pdf"), height = 4, width = 5)
ggsave(paste0(plotDir, "expPCA_fancy.png"), height = 4, width = 5)

cor.test(expCitMat$Rv2349c,  expCitMat$late_apoptotic)
cor.test(expCitMat$Rv2350c,  expCitMat$late_apoptotic)

cor.test(expCitMat$Rv2349c,  expCitMat$necrotic)
cor.test(expCitMat$Rv2350c,  expCitMat$necrotic)

# Do PCA of expression data without bovis
expMatNoBov_pca <- prcomp(expCitMat[sapply(rownames(expCitMat), function(x) strsplit(x, "_")[[1]][1] != "Bovis"), 
                                    colnames(expCitMat) %in% c("Rv2349c", "Rv2350c", "Rv2351c", "Mb1784c")],
                          scale. = T)

pcBiplot(expMatNoBov_pca)

ggsave(paste0(plotDir, "expPCA_noBov_fancy.pdf"), height = 4, width = 5)
ggsave(paste0(plotDir, "expPCA_noBov_fancy.png"), height = 4, width = 5)

# Do PCA only of necrosis and viability to see if we can obtain a one dimensional measure for virulence

pcBiplot(prcomp(expCitMat[, colnames(expCitMat) %in% c("necrotic", "viable")], scale. = T, center = T))

pcBiplot(prcomp(citoNum_d6[, colnames(citoNum_d6) %in% c("necrotic", "viable")], scale. = T, center = T))



# Function for rotating the data a certain angle (in degrees).
rotateMat <- function(inDat, angle, varX, varY, sampMax = "R0003", sampMin = "R0211"){
        inDat <- inDat[, c(varX, varY)]
        angle <- angle*pi/180
        inDat <- as.matrix(inDat)
        rotMat <- matrix(c(cos(angle),
                           sin(angle),
                           -sin(angle),
                           cos(angle)),
                         nrow = 2, 
                         ncol = 2)
        rotData <- t(rotMat %*% t(inDat) )
        newMin <- rotData[grep(sampMin, rownames(rotData)), 1]
        newMax <- rotData[grep(sampMax, rownames(rotData)), 1]
        print(newMax + newMin)
        return(rotData)
        
}

# We estimated (roughly and visually through rotation) that a counter clockwise rotation of 40 degree of the 
# viability/necrosis plot produced a maximal spread of the data accross the x axis (not in terms of maximim difference)
# between max and min, but in terms of having them separated but the intermediate phenotypes in the middle. With it 
# We obtained the following function for getting a virulence measure given necrosis and apoptosis

getVirulence <- function(necrosis, viability, rotAngle = 40){
        virulence <- necrosis * cos(rotAngle*pi/180) - viability * sin(rotAngle*pi/180)
        return(virulence)
}



plot(rotateMat(expCitMat, 60, varX = "necrotic", varY = "viable", ), ylim = c(-6, 6), xlim = c(-6, 6))

plot(rotateMat(citoNum_d6, 40,
               varX = "necrotic",
               varY = "viable",
               sampMax = "L5_THP-1_R3",
               sampMin = "Bovis_THP-1_R4"), ylim = c(-6, 6), xlim = c(-6, 6))

plot(expCitMat$necrotic, expCitMat$viable, ylim = c(-6, 6), xlim = c(-6, 6))

as.matrix(expCitMat)

plotViabl <- function(inDat){
        data$Bacteria <- sapply(data$obsnames, function(x) strsplit(as.character(x), "_")[[1]][1])
        data$Host <- sapply(data$obsnames, function(x) strsplit(as.character(x), "_")[[1]][2])
        data$Host <- factor(data$Host, levels = c("THP-1", "BOMAC"))
        bactLevls <- c("L5", "L6", "Bovis", "Chimp")[c("L5", "L6", "Bovis", "Chimp") %in% data$Bacteria]
        data$Bacteria <- factor(data$Bacteria, levels = bactLevls)
}

######################################################################################################################
#                                                                                                                    #
# Linear models of day 1 post infection expression data Vs day 2 post infection citometry data                       #
#                                                                                                                    #
######################################################################################################################

# Multiple regression of BOMAC and THP-1 at the same time 
######################################################################################################################

# Scale and center the data
expCitMat$virulence <- getVirulence(expCitMat$necrotic, expCitMat$viable, rotAngle = 40)
expCitMat_scaled <- as.data.frame(scale(expCitMat))
 
# Fit linear models to each type of death
lateApFit <- lm(formula("late_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_scaled)
summary(lateApFit)
viableFit <- lm(formula("viable ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_scaled)
summary(viableFit)
earlApFit <- lm(formula("early_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_scaled)
summary(earlApFit)
necrotFit <- lm(formula("necrotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_scaled)
summary(necrotFit)
virulFit <- lm(formula("virulence ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_scaled)
summary(virulFit)

# Try to do the fit excluding bovis, as it doesnt have plcABC genes and is a set of strains with low
# virulence. 
virulFit_noBov <- lm(formula("virulence ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), 
                     data = expCitMat_scaled[!grepl("bovis", tolower(rownames(expCitMat_scaled))), ])

summary(virulFit_noBov)

virulFitPlcC <- lm(formula("virulence ~ Rv2349c"), data = expCitMat_scaled)

virulOPLSDA <- opls(expCitMat_scaled[, colnames(expCitMat_scaled) %in% c("Rv2349c", "Rv2350c", "Rv2351c", "Mb1784c")], 
                    expCitMat_scaled$virulence,
                    orthoI = 2,
                    predI = 1)
getVipVn(virulOPLSDA)
getLoadingMN(virulOPLSDA)


expMat_wNuGenes <- getExpMat(results, 
                             plateNorm = T, 
                             genes = c("Rv1106", "Rv2349c", "Rv2350c", "Rv2351c", "Mb1784c", "Rv2379c", "Rv2383c", "Rv3545c", "mcr7"))

expMat_wNuGenes$virulence <- expCitMat$virulence[match(rownames(expMat_wNuGenes),
                                                       rownames(expCitMat))]

expMat_wNuGenes <- expMat_wNuGenes[expMat_wNuGenes$mcr7 != 0, ]

expMat_wNuGenes <- expMat_wNuGenes[, colnames(expMat_wNuGenes) %in% c("Rv1106", "Rv2349c", "Rv2350c", "Rv2351c", "Mb1784c", "Rv2379c", "Rv2383c", "Rv3545c", "mcr7", "virulence")]

expMat_wNuGenes_scaled <- as.data.frame(scale(expMat_wNuGenes))

virulFit_wNuGenes <- lm(formula = "virulence ~ .", data = expMat_wNuGenes_scaled)
summary(virulFit_wNuGenes)

doAddVarPlots1(virulFit_wNuGenes, variable = "Rv2379c", formatGeneNames = F)
doAddVarPlots1(virulFit_wNuGenes, variable = "Rv3545c", formatGeneNames = F)
doAddVarPlots1(virulFit_wNuGenes, variable = "mcr7", formatGeneNames = F)


doMergedAddVarPlots(virulFit_wNuGenes)


virulOPLSDA_wNuGenes <- opls(expMat_wNuGenes_scaled[, colnames(expMat_wNuGenes_scaled) %in% c("Rv1106", "Rv2349c", "Rv2350c", "Rv2351c", "Mb1784c", "Rv2379c", "Rv2383c", "Rv3545c", "mcr7")], 
                             expMat_wNuGenes_scaled$virulence,
                             orthoI = 2,
                             predI = 1)



data.frame(apply(expMat_wNuGenes, 2, function(x) (x - mean(x))/sd(x)))


# Fit models to the randomly generated expression-citometry datasets, to see what models and genes are 
# statistically significant at higher proportions
getRandRegsFreqs <- function(randMatLst, getVirulence = T){
        freqMat <- data.frame(matrix(rep(0, 5*5),
                                     nrow = 5,
                                     ncol = 5,
                                     dimnames = list(c("lateAp",
                                                       "viable",
                                                       "earlAp",
                                                       "necrot",
                                                       "virule"),
                                                     c("modelCount",
                                                       "Rv2349c",
                                                       "Rv2350c",
                                                       "Rv2351c",
                                                       "Mb1784c"))))
        for(i in seq_along(randMatLst)){
                randMat <- randMatLst[[i]]
                if(getVirulence == T){
                        randMat$virulence <- getVirulence(randMat$necrotic, randMat$viable, rotAngle = 40)
                }
                randMat <- as.data.frame(scale(randMat))
                
                lateAp <- lm(formula("late_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = randMat)
                lateApSum <- summary(lateAp)
                lateApPVal <- pf(lateApSum$fstatistic[1],
                                 lateApSum$fstatistic[2],
                                 lateApSum$fstatistic[3],
                                 lower.tail = F)
                lateApGenePVals <- lateApSum[[4]][, 4][2:5]
                if(lateApPVal <= 0.05){
                        freqMat[1, 1] <- freqMat[1, 1] + 1
                }
                freqMat[1, 2:5][lateApGenePVals <= 0.05] <- freqMat[1, 2:5][lateApGenePVals <= 0.05] + 1
                
                viable <- lm(formula("viable ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = randMat)
                viableSum <- summary(viable)
                viablePVal <- pf(viableSum$fstatistic[1],
                                 viableSum$fstatistic[2],
                                 viableSum$fstatistic[3],
                                 lower.tail = F)
                viableGenePVals <- viableSum[[4]][, 4][2:5]
                if(viablePVal <= 0.05){
                        freqMat[2, 1] <- freqMat[2, 1] + 1
                }
                freqMat[2, 2:5][viableGenePVals <= 0.05] <- freqMat[2, 2:5][viableGenePVals <= 0.05] + 1
                
                earlAp <- lm(formula("early_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = randMat)
                earlApSum <- summary(earlAp)
                earlApPVal <- pf(earlApSum$fstatistic[1],
                                 earlApSum$fstatistic[2],
                                 earlApSum$fstatistic[3],
                                 lower.tail = F)
                earlApGenePVals <- earlApSum[[4]][, 4][2:5]
                if(earlApPVal <= 0.05){
                        freqMat[3, 1] <- freqMat[3, 1] + 1
                }
                freqMat[3, 2:5][earlApGenePVals <= 0.05] <- freqMat[3, 2:5][earlApGenePVals <= 0.05] + 1
                
                necrot <- lm(formula("necrotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = randMat)
                necrotSum <- summary(necrot)
                necrotPVal <- pf(necrotSum$fstatistic[1],
                                 necrotSum$fstatistic[2],
                                 necrotSum$fstatistic[3],
                                 lower.tail = F)
                necrotGenePVals <- necrotSum[[4]][, 4][2:5]
                if(necrotPVal <= 0.05){
                        freqMat[4, 1] <- freqMat[4, 1] + 1
                }
                freqMat[4, 2:5][necrotGenePVals <= 0.05] <- freqMat[4, 2:5][necrotGenePVals <= 0.05] + 1
                
                virule <- lm(formula("virulence ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = randMat)
                viruleSum <- summary(virule)
                virulePVal <- pf(viruleSum$fstatistic[1],
                                 viruleSum$fstatistic[2],
                                 viruleSum$fstatistic[3],
                                 lower.tail = F)
                viruleGenePVals <- viruleSum[[4]][, 4][2:5]
                if(virulePVal <= 0.05){
                        freqMat[5, 1] <- freqMat[5, 1] + 1
                }
                freqMat[5, 2:5][viruleGenePVals <= 0.05] <- freqMat[5, 2:5][viruleGenePVals <= 0.05] + 1
        }
        return(freqMat)
}

getRandRegsFreqs(randExpCitMatLst)

getRandRegsFreqs(randExpCitMatLst_d1d6)



# The result of necrosis and virulence fit is the same in all of the iterations. Late apoptosis 
# is significant in 680 out of 1000 iterations. The rest don't happen a lot.

# Do added variable plots for the models.
lateApAll <- doMergedAddVarPlots(lateApFit)

lateApAll

ggsave(paste0(plotDir, "lateAp_addVar.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "lateAp_addVar.pdf"), height = 7.5, width = 6.75)


earlApAll <- doMergedAddVarPlots(earlApFit)

earlApAll

ggsave(paste0(plotDir, "earlAp_addVar.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "earlAp_addVar.pdf"), height = 7.5, width = 6.75)


viableAll <- doMergedAddVarPlots(viableFit)

viableAll

ggsave(paste0(plotDir, "viable_addVar.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "viable_addVar.pdf"), height = 7.5, width = 6.75)


necrotAll <- doMergedAddVarPlots(necrotFit)

necrotAll

ggsave(paste0(plotDir, "necrot_addVar.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "necrot_addVar.pdf"), height = 7.5, width = 6.75)


virulAll <- doMergedAddVarPlots(virulFit)

virulAll

ggsave(paste0(plotDir, "virul_addVar.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "virul_addVar.pdf"), height = 7.5, width = 6.75)


# Multiple regression of BOMAC infections
######################################################################################################################

expCitMat_bomac <- expCitMat[grep("BOMAC", rownames(expCitMat)), ]

# Scale and center the data
expCitMat_bomac_scaled <- as.data.frame(scale(expCitMat_bomac))

# Fit linear models to each type of death
lateApFit_bomac <- lm(formula("late_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_bomac_scaled)
summary(lateApFit_bomac)
viableFit_bomac <- lm(formula("viable ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_bomac_scaled)
summary(viableFit_bomac)
earlApFit_bomac <- lm(formula("early_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_bomac_scaled)
summary(earlApFit_bomac)
necrotFit_bomac <- lm(formula("necrotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_bomac_scaled)
summary(necrotFit_bomac)

# Do added variable plots for the models.

lateApAll_bomac <- doMergedAddVarPlots(lateApFit_bomac)

lateApAll_bomac

ggsave(paste0(plotDir, "lateAp_addVar_bomac.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "lateAp_addVar_bomac.pdf"), height = 7.5, width = 6.75)


earlApAll_bomac <- doMergedAddVarPlots(earlApFit_bomac)

earlApAll_bomac

ggsave(paste0(plotDir, "earlAp_addVar_bomac.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "earlAp_addVar_bomac.pdf"), height = 7.5, width = 6.75)


viableAll_bomac <- doMergedAddVarPlots(viableFit_bomac)

viableAll_bomac

ggsave(paste0(plotDir, "viable_addVar_bomac.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "viable_addVar_bomac.pdf"), height = 7.5, width = 6.75)


necrotAll_bomac <- doMergedAddVarPlots(necrotFit_bomac)

necrotAll_bomac

ggsave(paste0(plotDir, "necrot_addVar_bomac.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "necrot_addVar_bomac.pdf"), height = 7.5, width = 6.75)

# Multiple regression of THP-1 infections
######################################################################################################################

expCitMat_thp1 <- expCitMat[grep("THP-1", rownames(expCitMat)), ]

# Scale and center the data
expCitMat_thp1_scaled <- as.data.frame(scale(expCitMat_thp1))

# Fit linear models to each type of death
lateApFit_thp1 <- lm(formula("late_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_thp1_scaled)
summary(lateApFit_thp1)
viableFit_thp1 <- lm(formula("viable ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_thp1_scaled)
summary(viableFit_thp1)
earlApFit_thp1 <- lm(formula("early_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_thp1_scaled)
summary(earlApFit_thp1)
necrotFit_thp1 <- lm(formula("necrotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_thp1_scaled)
summary(necrotFit_thp1)

# Do added variable plots for the models.

lateApAll_thp1 <- doMergedAddVarPlots(lateApFit_thp1)

lateApAll_thp1

ggsave(paste0(plotDir, "lateAp_addVar_thp1.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "lateAp_addVar_thp1.pdf"), height = 7.5, width = 6.75)


earlApAll_thp1 <- doMergedAddVarPlots(earlApFit_thp1)

earlApAll_thp1

ggsave(paste0(plotDir, "earlAp_addVar_thp1.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "earlAp_addVar_thp1.pdf"), height = 7.5, width = 6.75)


viableAll_thp1 <- doMergedAddVarPlots(viableFit_thp1)

viableAll_thp1

ggsave(paste0(plotDir, "viable_addVar_thp1.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "viable_addVar_thp1.pdf"), height = 7.5, width = 6.75)


necrotAll_thp1 <- doMergedAddVarPlots(necrotFit_thp1)

necrotAll_thp1

ggsave(paste0(plotDir, "necrot_addVar_thp1.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "necrot_addVar_thp1.pdf"), height = 7.5, width = 6.75)

#
# Fit linear models to each type of cell death, separating by strain
######################################################################################################################

expCitMat_L5 <- expCitMat[grep("L5", rownames(expCitMat)), ]
expCitMat_L6 <- expCitMat[grep("L6", rownames(expCitMat)), ]
expCitMat_bov <- expCitMat[grep("Bovis", rownames(expCitMat)), ]
expCitMat_cmp <- expCitMat[grep("Chimp", rownames(expCitMat)), ]

# Scale and center the data
expCitMat_L5_scaled <- as.data.frame(scale(expCitMat_L5))

# Fit linear models to each type of death
lateApFit_L5 <- lm(formula("late_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_L5_scaled)
summary(lateApFit_L5)
viableFit_L5 <- lm(formula("viable ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_L5_scaled)
summary(viableFit_L5)
earlApFit_L5 <- lm(formula("early_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_L5_scaled)
summary(earlApFit_L5)
necrotFit_L5 <- lm(formula("necrotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_L5_scaled)
summary(necrotFit_L5)

# Scale and center the data
expCitMat_L6_scaled <- as.data.frame(scale(expCitMat_L6))

# Fit linear models to each type of death
lateApFit_L6 <- lm(formula("late_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_L6_scaled)
summary(lateApFit_L6)
viableFit_L6 <- lm(formula("viable ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_L6_scaled)
summary(viableFit_L6)
earlApFit_L6 <- lm(formula("early_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_L6_scaled)
summary(earlApFit_L6)
necrotFit_L6 <- lm(formula("necrotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_L6_scaled)
summary(necrotFit_L6)

# Scale and center the data
expCitMat_bov_scaled <- as.data.frame(scale(expCitMat_bov))

expCitMat_bov_scaled <- expCitMat_bov_scaled[, colnames(expCitMat_bov_scaled) != "Rv2349c"]

# Fit linear models to each type of death
lateApFit_bov <- lm(formula("late_apoptotic ~ Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_bov_scaled)
summary(lateApFit_bov)
viableFit_bov <- lm(formula("viable ~ Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_bov_scaled)
summary(viableFit_bov)
earlApFit_bov <- lm(formula("early_apoptotic ~ Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_bov_scaled)
summary(earlApFit_bov)
necrotFit_bov <- lm(formula("necrotic ~ Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_bov_scaled)
summary(necrotFit_bov)

# Scale and center the data
expCitMat_cmp_scaled <- as.data.frame(scale(expCitMat_cmp))

# Fit linear models to each type of death
lateApFit_cmp <- lm(formula("late_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_cmp_scaled)
summary(lateApFit_cmp)
viableFit_cmp <- lm(formula("viable ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_cmp_scaled)
summary(viableFit_cmp)
earlApFit_cmp <- lm(formula("early_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_cmp_scaled)
summary(earlApFit_cmp)
necrotFit_cmp <- lm(formula("necrotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_cmp_scaled)
summary(necrotFit_cmp)

# Do AddVar plots
lateApAll_L5 <- doMergedAddVarPlots(lateApFit_L5)

lateApAll_L5

ggsave(paste0(plotDir, "lateAp_addVar_L5.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "lateAp_addVar_L5.pdf"), height = 7.5, width = 6.75)

lateApAddvar2_L5 <- doAddVarPlots2(lateApFit_L5, expCitMat_L5_scaled)
lateApAddvar2_L5$globPlot

ggsave(paste0(plotDir, "lateAp_addVar2_L5.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "lateAp_addVar2_L5.pdf"), height = 7.5, width = 6.75)


earlApAll_L5 <- doMergedAddVarPlots(earlApFit_L5)

earlApAll_L5

ggsave(paste0(plotDir, "earlAp_addVar_L5.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "earlAp_addVar_L5.pdf"), height = 7.5, width = 6.75)




earlApAddvar2_L5 <- doAddVarPlots2(earlApFit_L5, expCitMat_L5_scaled)
earlApAddvar2_L5$globPlot

ggsave(paste0(plotDir, "earlAp_addVar2_L5.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "earlAp_addVar2_L5.pdf"), height = 7.5, width = 6.75)


necrotAll_L5 <- doMergedAddVarPlots(necrotFit_L5)

necrotAll_L5

ggsave(paste0(plotDir, "necrot_addVar_L5.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "necrot_addVar_L5.pdf"), height = 7.5, width = 6.75)


viableAll_L5 <- doMergedAddVarPlots(viableFit_L5)

viableAll_L5

ggsave(paste0(plotDir, "viable_addVar_L5.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "viable_addVar_L5.pdf"), height = 7.5, width = 6.75)


lateApAll_L6 <- doMergedAddVarPlots(lateApFit_L6)

lateApAll_L6

ggsave(paste0(plotDir, "lateAp_addVar_L6.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "lateAp_addVar_L6.pdf"), height = 7.5, width = 6.75)


earlApAll_L6 <- doMergedAddVarPlots(earlApFit_L6)

earlApAll_L6

ggsave(paste0(plotDir, "earlAp_addVar_L6.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "earlAp_addVar_L6.pdf"), height = 7.5, width = 6.75)


necrotAll_L6 <- doMergedAddVarPlots(necrotFit_L6)

necrotAll_L6

ggsave(paste0(plotDir, "necrot_addVar_L6.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "necrot_addVar_L6.pdf"), height = 7.5, width = 6.75)


viableAll_L6 <- doMergedAddVarPlots(viableFit_L6)

viableAll_L6

ggsave(paste0(plotDir, "viable_addVar_L6.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "viable_addVar_L6.pdf"), height = 7.5, width = 6.75)


lateApAll_bov <- doMergedAddVarPlots(lateApFit_bov)

lateApAll_bov

ggsave(paste0(plotDir, "lateAp_addVar_bov.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "lateAp_addVar_bov.pdf"), height = 7.5, width = 6.75)


earlApAll_bov <- doMergedAddVarPlots(earlApFit_bov)

earlApAll_bov

ggsave(paste0(plotDir, "earlAp_addVar_bov.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "earlAp_addVar_bov.pdf"), height = 7.5, width = 6.75)

necrotAll_bov <- doMergedAddVarPlots(necrotFit_bov)

necrotAll_bov

ggsave(paste0(plotDir, "necrot_addVar_bov.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "necrot_addVar_bov.pdf"), height = 7.5, width = 6.75)


viableAll_bov <- doMergedAddVarPlots(viableFit_bov)

viableAll_bov

ggsave(paste0(plotDir, "viable_addVar_bov.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "viable_addVar_bov.pdf"), height = 7.5, width = 6.75)


lateApAll_cmp <- doMergedAddVarPlots(lateApFit_cmp)

lateApAll_cmp

ggsave(paste0(plotDir, "lateAp_addVar_cmp.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "lateAp_addVar_cmp.pdf"), height = 7.5, width = 6.75)


earlApAll_cmp <- doMergedAddVarPlots(earlApFit_cmp)

earlApAll_cmp

ggsave(paste0(plotDir, "earlAp_addVar_cmp.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "earlAp_addVar_cmp.pdf"), height = 7.5, width = 6.75)


necrotAll_cmp <- doMergedAddVarPlots(necrotFit_cmp)

necrotAll_cmp

ggsave(paste0(plotDir, "necrot_addVar_cmp.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "necrot_addVar_cmp.pdf"), height = 7.5, width = 6.75)


viableAll_cmp <- doMergedAddVarPlots(viableFit_cmp)

viableAll_cmp

ggsave(paste0(plotDir, "viable_addVar_cmp.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "viable_addVar_cmp.pdf"), height = 7.5, width = 6.75)

######################################################################################################################
#                                                                                                                    #
# OPLS-DA                                                                                                            #
#                                                                                                                    #
######################################################################################################################

# Use OPLS-DA models to determine association to host. 
oplsDA_expVsNecr <- opls(expCitMat_scaled[1:4], y = expCitMat_scaled$necrotic, orthoI = 3, predI = 1, permI = 100)

getVipVn(oplsDA_expVsNecr)

expCitMatHost <- expCitMat
expCitMatHost <- as.data.frame(scale(expCitMatHost))

expCitMatHost$host <- sapply(rownames(expCitMatHost), function(x) strsplit(x, "_")[[1]][2])

oplsDA_expVsHost <- opls(expCitMatHost[1:4], y = expCitMatHost$host, orthoI = 3, predI = 1, permI = 100)

getVipVn(oplsDA_expVsHost)

######################################################################################################################
#                                                                                                                    #
# Logistic Regressions                                                                                               #
#                                                                                                                    #
######################################################################################################################

# Multiple logistic regression with all the strains at the same time 
######################################################################################################################

expMat_4Logit <- expCitMat_scaled[, 1:4]
expMat_4Logit$host <- factor(sapply(rownames(expMat_4Logit), function(x) strsplit(x, "_")[[1]][2]))

expMat_4Logit_logit <- doLogit(".", expMat_4Logit)

summary(expMat_4Logit_logit)

# Univariate logistic regression with all the strains at the same time 
######################################################################################################################

expMat_4Logit_Mb1784c_logit <- doLogit("Mb1784c", expMat_4Logit)
expMat_4Logit_Rv2349c_logit <- doLogit("Rv2349c", expMat_4Logit)
expMat_4Logit_Rv2350c_logit <- doLogit("Rv2350c", expMat_4Logit)
expMat_4Logit_Rv2351c_logit <- doLogit("Rv2351c", expMat_4Logit)

summary(expMat_4Logit_Mb1784c_logit)

# Plot the result of the regressions

# Mb1784c
expMat_4Logit_Mb1784c_logitPlot <- plotLogit(expMat_4Logit_Mb1784c_logit)

expMat_4Logit_Mb1784c_logitPlot

ggsave(paste0(plotDir, "Mb1784_Host_logit.png"), height = 4.5, width = 5.5)
ggsave(paste0(plotDir, "Mb1784_Host_logit.pdf"), height = 4.5, width = 5.5)

# Rv2349c
expMat_4Logit_Rv2349c_logitPlot <- plotLogit(expMat_4Logit_Rv2349c_logit)

expMat_4Logit_Rv2349c_logitPlot

ggsave(paste0(plotDir, "Rv2349c_Host_logit.png"), height = 4.5, width = 5.5)
ggsave(paste0(plotDir, "Rv2349c_Host_logit.pdf"), height = 4.5, width = 5.5)

# Rv2350c
expMat_4Logit_Rv2350c_logitPlot <- plotLogit(expMat_4Logit_Rv2350c_logit)

expMat_4Logit_Rv2350c_logitPlot

ggsave(paste0(plotDir, "Rv2350c_Host_logit.png"), height = 4.5, width = 5.5)
ggsave(paste0(plotDir, "Rv2350c_Host_logit.pdf"), height = 4.5, width = 5.5)

# Rv2351c
expMat_4Logit_Rv2351c_logitPlot <- plotLogit(expMat_4Logit_Rv2350c_logit)

expMat_4Logit_Rv2351c_logitPlot

ggsave(paste0(plotDir, "Rv2351c_Host_logit.png"), height = 4.5, width = 5.5)
ggsave(paste0(plotDir, "Rv2351c_Host_logit.pdf"), height = 4.5, width = 5.5)

ggarrange(expMat_4Logit_Rv2351c_logitPlot,
          expMat_4Logit_Rv2350c_logitPlot + 
                  rremove("ylab"), 
          expMat_4Logit_Rv2349c_logitPlot, 
          expMat_4Logit_Mb1784c_logitPlot + 
                  rremove("ylab"),
          legend = "bottom",
          common.legend = T)


ggsave(paste0(plotDir, "plcABCD_Host_logit.png"), height = 9, width = 10)
ggsave(paste0(plotDir, "plcABCD_Host_logit.pdf"), height = 9, width = 10)

# Multivariate and Univariate logistic regression with separating by strain 
######################################################################################################################

# L5
expMat_L5 <- expMat[grep("L5", rownames(expMat)), ]
expMat_L5 <- as.data.frame(scale(expMat_L5))
expMat_L5$host <- factor(sapply(rownames(expMat_L5), function(x) strsplit(x, "_")[[1]][2]))

pcBiplot(prcomp(expMat_L5[, 1:4]))

# All variables at the same time 
expMat_L5_fit <- doLogit(".", data = expMat_L5)

summary(expMat_L5_fit)

# Separated fits for each variable (gene)
expMat_L5_Rv2349c_logit <- doLogit("Rv2349c", data = expMat_L5)

summary(expMat_L5_Rv2349c_logit)

L5_Rv2349c_logit <- plotLogit(expMat_L5_Rv2349c_logit)


expMat_L5_Rv2350c_logit <- doLogit("Rv2350c", data = expMat_L5)

summary(expMat_L5_Rv2350c_logit)

L5_Rv2350c_logit <- plotLogit(expMat_L5_Rv2350c_logit)


expMat_L5_Rv2351c_logit <- doLogit("Rv2351c", data = expMat_L5)

summary(expMat_L5_Rv2351c_logit)

L5_Rv2351c_logit <- plotLogit(expMat_L5_Rv2351c_logit)


expMat_L5_Mb1784c_logit <- doLogit("Mb1784c", data = expMat_L5)

summary(expMat_L5_Mb1784c_logit)

L5_Mb1784c_logit <- plotLogit(expMat_L5_Mb1784c_logit)

ggarrange(L5_Rv2351c_logit,
          L5_Rv2350c_logit + 
                  rremove("ylab"), 
          L5_Rv2349c_logit, 
          L5_Mb1784c_logit + 
                  rremove("ylab"),
          legend = "bottom",
          common.legend = T)

ggsave(paste0(plotDir, "plcABCD_Host_L5_logit.png"), height = 9, width = 10)
ggsave(paste0(plotDir, "plcABCD_Host_L5_logit.pdf"), height = 9, width = 10)

# L6
expMat_L6 <- expMat[grep("L6", rownames(expMat)), ]
expMat_L6 <- as.data.frame(scale(expMat_L6))
expMat_L6$host <- factor(sapply(rownames(expMat_L6), function(x) strsplit(x, "_")[[1]][2]))

pcBiplot(prcomp(expMat_L6[, 1:4]))

# All variables at the same time 
expMat_L6_fit <- doLogit(".", data = expMat_L6)

summary(expMat_L6_fit)

# Separated fits for each variable (gene)
expMat_L6_Rv2349c_logit <- doLogit("Rv2349c", data = expMat_L6)

summary(expMat_L6_Rv2349c_logit)

L6_Rv2349c_logit <- plotLogit(expMat_L6_Rv2349c_logit)


expMat_L6_Rv2350c_logit <- doLogit("Rv2350c", data = expMat_L6)

summary(expMat_L6_Rv2350c_logit)

L6_Rv2350c_logit <- plotLogit(expMat_L6_Rv2350c_logit)


expMat_L6_Rv2351c_logit <- doLogit("Rv2351c", data = expMat_L6)

summary(expMat_L6_Rv2351c_logit)

L6_Rv2351c_logit <- plotLogit(expMat_L6_Rv2351c_logit)


expMat_L6_Mb1784c_logit <- doLogit("Mb1784c", data = expMat_L6)

summary(expMat_L6_Mb1784c_logit)

L6_Mb1784c_logit <- plotLogit(expMat_L6_Mb1784c_logit)

ggarrange(L6_Rv2351c_logit,
          L6_Rv2350c_logit + 
                  rremove("ylab"), 
          L6_Rv2349c_logit, 
          L6_Mb1784c_logit + 
                  rremove("ylab"),
          legend = "bottom",
          common.legend = T)

ggsave(paste0(plotDir, "plcABCD_Host_L6_logit.png"), height = 9, width = 10)
ggsave(paste0(plotDir, "plcABCD_Host_L6_logit.pdf"), height = 9, width = 10)

# bovis
expMat_bovis <- expMat[grep("Bovis", rownames(expMat)), ]
expMat_bovis <- as.data.frame(scale(expMat_bovis))
expMat_bovis$host <- factor(sapply(rownames(expMat_bovis), function(x) strsplit(x, "_")[[1]][2]))

pcBiplot(prcomp(expMat_bovis[, c(1, 3, 4)]))

# All variables at the same time 
expMat_bovis_fit <- doLogit(".", data = expMat_bovis)

summary(expMat_bovis_fit)

# Separated fits for each variable (gene)
expMat_bovis_Rv2349c_logit <- doLogit("Rv2349c", data = expMat_bovis)

summary(expMat_bovis_Rv2349c_logit)

bov_Rv2349c_logit <- plotLogit(expMat_bovis_Rv2349c_logit)


expMat_bovis_Rv2350c_logit <- doLogit("Rv2350c", data = expMat_bovis)

summary(expMat_bovis_Rv2350c_logit)

bov_Rv2350c_logit <- plotLogit(expMat_bovis_Rv2350c_logit)


expMat_bovis_Rv2351c_logit <- doLogit("Rv2351c", data = expMat_bovis)

summary(expMat_bovis_Rv2351c_logit)

bov_Rv2351c_logit <- plotLogit(expMat_bovis_Rv2351c_logit)


expMat_bovis_Mb1784c_logit <- doLogit("Mb1784c", data = expMat_bovis)

summary(expMat_bovis_Mb1784c_logit)

bov_Mb1784c_logit <- plotLogit(expMat_bovis_Mb1784c_logit)

ggarrange(bov_Rv2351c_logit,
          bov_Rv2350c_logit + 
                  rremove("ylab"), 
          #bov_Rv2349c_logit, 
          bov_Mb1784c_logit + 
                  rremove("ylab"),
          legend = "bottom",
          common.legend = T)

ggsave(paste0(plotDir, "plcABCD_Host_bov_logit.png"), height = 9, width = 10)
ggsave(paste0(plotDir, "plcABCD_Host_bov_logit.pdf"), height = 9, width = 10)

# chimp
expMat_chimp <- expMat[grep("Chimp", rownames(expMat)), ]
expMat_chimp <- as.data.frame(scale(expMat_chimp))
expMat_chimp$host <- factor(sapply(rownames(expMat_chimp), function(x) strsplit(x, "_")[[1]][2]))

pcBiplot(prcomp(expMat_chimp[, 1:4]))

# All variables at the same time
expMat_chimp_fit <- doLogit(".", data = expMat_chimp)

summary(expMat_chimp_fit)

# Separated fits for each variable (gene)
expMat_chimp_Rv2349c_logit <- doLogit("Rv2349c", data = expMat_chimp)

summary(expMat_chimp_Rv2349c_logit)

cmp_Rv2349c_logit <- plotLogit(expMat_chimp_Rv2349c_logit)


expMat_chimp_Rv2350c_logit <- doLogit("Rv2350c", data = expMat_chimp)

summary(expMat_chimp_Rv2350c_logit)

cmp_Rv2350c_logit <- plotLogit(expMat_chimp_Rv2350c_logit)


expMat_chimp_Rv2351c_logit <- doLogit("Rv2351c", data = expMat_chimp)

summary(expMat_chimp_Rv2351c_logit)

cmp_Rv2351c_logit <- plotLogit(expMat_chimp_Rv2351c_logit)


expMat_chimp_Mb1784c_logit <- doLogit("Mb1784c", data = expMat_chimp)

summary(expMat_chimp_Mb1784c_logit)

cmp_Mb1784c_logit <- plotLogit(expMat_chimp_Mb1784c_logit)

ggarrange(cmp_Rv2351c_logit,
          cmp_Rv2350c_logit + 
                  rremove("ylab"), 
          cmp_Rv2349c_logit, 
          cmp_Mb1784c_logit + 
                  rremove("ylab"),
          legend = "bottom",
          common.legend = T)

ggsave(paste0(plotDir, "plcABCD_Host_cmp_logit.png"), height = 9, width = 10)
ggsave(paste0(plotDir, "plcABCD_Host_cmp_logit.pdf"), height = 9, width = 10)


######################################################################################################################
#                                                                                                                    #
# Francisco's models                                                                                                 #
#                                                                                                                    #
######################################################################################################################

# Create dataframe with host adaptation variables

franciscoDF <- expCitMat

franciscoDF$lin <- factor(sapply(rownames(franciscoDF), function(x) strsplit(x, "_")[[1]][1]))

franciscoDF$host <- factor(sapply(rownames(franciscoDF), function(x) strsplit(x, "_")[[1]][2]))

franciscoDF$adapt <- rep(NA, nrow(franciscoDF))

adapted <- c("L5_THP-1", "Bovis_BOMAC")
midAdapted <- c("L6_THP-1", "L6_BOMAC", "Chimp_THP-1", "Chimp_BOMAC")
notAdapted <- c("L5_BOMAC", "Bovis_THP-1")
for(i in 1:nrow(franciscoDF)){
        
        case <- paste(franciscoDF$lin[i],
                      franciscoDF$host[i],
                      sep = "_")
        if(case %in% adapted){
                franciscoDF$adapt[i] <- "adapted"
        }else if(case %in% midAdapted){
                franciscoDF$adapt[i] <- "midAdapted"
        }else if(case %in% notAdapted){
                franciscoDF$adapt[i] <- "notAdapted"
        }
}
franciscoDF$adapt <- factor(franciscoDF$adapt, levels = sort(unique(franciscoDF$adapt), decreasing = T))

franciscoDF$totalPlc <- rowSums(franciscoDF[, 1:4])

summary(glm(totalPlc ~ lin + host + adapt, data = franciscoDF))

summary(glm(Rv2349c ~ lin + host + adapt, data = franciscoDF))

summary(glm(Rv2350c ~ lin + host + adapt, data = franciscoDF))

summary(glm(Rv2351c ~ lin + host + adapt, data = franciscoDF))

summary(glm(Mb1784c ~ lin + host + adapt, data = franciscoDF))

summary(glm(virulence ~ lin + host + Rv2349c + Rv2350c + Rv2351c + Mb1784c, data = franciscoDF))

summary(glm(virulence ~ adapt + Rv2349c + Rv2350c + Rv2351c + Mb1784c, data = franciscoDF))

summary(glm(adapt ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c, data = franciscoDF))


franciscoDF$adapt <- relevel(franciscoDF$adapt, ref = "midAdapted")
summary(multinom(adapt ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c, data = franciscoDF))

######################################################################################################################
#                                                                                                                    #
# Linear models of day 3 post infection expression data Vs day 6 post infection citometry data                       #
#                                                                                                                    #
######################################################################################################################

# Scale and center the data
expCitMat_D3D6_scaled <- as.data.frame(scale(expCitMat_d3_d6))

# Fit linear models to each type of death
lateApFit_D3D6 <- lm(formula("late_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_D3D6_scaled)
summary(lateApFit_D3D6)
viableFit_D3D6 <- lm(formula("viable ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_D3D6_scaled)
summary(viableFit_D3D6)
earlApFit_D3D6 <- lm(formula("early_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_D3D6_scaled)
summary(earlApFit_D3D6)
necrotFit_D3D6 <- lm(formula("necrotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_D3D6_scaled)
summary(necrotFit_D3D6)

# Do added variable plots for the models of day3 exp vs day 6 cito.

lateApAll_D3D6 <- doMergedAddVarPlots(lateApFit_D3D6)

lateApAll_D3D6

ggsave(paste0(plotDir, "lateAp_addVar_D3D6.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "lateAp_addVar_D3D6.pdf"), height = 7.5, width = 6.75)


earlApAll_D3D6 <- doMergedAddVarPlots(earlApFit_D3D6)

earlApAll_D3D6

ggsave(paste0(plotDir, "earlAp_addVar_D3D6.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "earlAp_addVar_D3D6.pdf"), height = 7.5, width = 6.75)


viableAll_D3D6 <- doMergedAddVarPlots(viableFit_D3D6)

viableAll_D3D6

ggsave(paste0(plotDir, "viable_addVar_D3D6.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "viable_addVar_D3D6.pdf"), height = 7.5, width = 6.75)


necrotAll_D3D6 <- doMergedAddVarPlots(necrotFit_D3D6)

necrotAll_D3D6

ggsave(paste0(plotDir, "necrot_addVar_D3D6.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "necrot_addVar_D3D6.pdf"), height = 7.5, width = 6.75)


#
#  Linear models of day 1 post infection expression data Vs day 6 post infection citometry data
######################################################################################################################

# Scale and center the data
expCitMat_d6$virulence <- getVirulence(expCitMat_d6$necrotic, expCitMat_d6$viable, rotAngle = 40)
#expCitMat_d6$virulence <- expCitMat_d6$necrotic * 0.48 - expCitMat_d6$viable * 0.78
#expCitMat_d6$virulence <- expCitMat_d6$necrotic - expCitMat_d6$viable
expCitMat_d6_scaled <- as.data.frame(scale(expCitMat_d6))

# Create a host preference variable with the rotation data of the PCA at day 6 post infection, to get the projections over the Y
# axis, as this axis separates the infections according to how well adapted is the bacteria to the host
expCitMat_d6_scaled$host_preference <- apply(expCitMat_d6_scaled[, 5:8], 1, function(x) sum(x * (citD6_pca$rotation[c(1:4), 2])))

#expCitMat_d6_scaled$virulence <- as.numeric(scale(apply(expCitMat_d6_scaled[, 5:8], 1, function(x) sum(x * c(-0.5, 0.5, 1, -1)))))
citD6_pca$rotation[, 2]

randExpCitMatLst_d1d6[[1]]

for(i in seq_along(randExpCitMatLst_d1d6)){
        randExpCitMatLst_d1d6[[i]]$virulence <- apply(scale(expCitMat_d6[, 5:8]), 1, function(x) sum(x * (citD6_pca$rotation[c(1:4), 2])))
        #randExpCitMatLst_d1d6[[i]]$virulence <- as.numeric(scale(apply(expCitMat_d6_scaled[, 5:8], 1, function(x) sum(x * c(-0.5, 0.5, 1, -1)))))
}

getRandRegsFreqs(randExpCitMatLst_d1d6, getVirulence = F)

lateApFit_d6 <- lm(formula("late_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(lateApFit_d6)
viableFit_d6 <- lm(formula("viable ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(viableFit_d6)
earlApFit_d6 <- lm(formula("early_apoptotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(earlApFit_d6)
necrotFit_d6 <- lm(formula("necrotic ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(necrotFit_d6)
hostPrFit_d6 <- lm(formula("host_preference ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(hostPrFit_d6)
hostPrFit_d6 <- lm(formula("virulence ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCitMat_d6_scaled)
summary(hostPrFit_d6)


expMat_wNuGenes_hostPr_d6 <- expMat_wNuGenes

expMat_wNuGenes_hostPr_d6$host_preference <- expCitMat_d6_scaled$host_preference[match(rownames(expMat_wNuGenes),
                                                                                       rownames(expCitMat_d6_scaled))]



hostPrFit_wNuGenes_d6 <- lm(formula("host_preference ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c + mcr7 + Rv1106 + Rv2379c + Rv2383c + Rv3545c"), 
                   data = data.frame(apply(expMat_wNuGenes_hostPr_d6, 2, function(x) (x - mean(x))/sd(x))))
summary(hostPrFit_wNuGenes_d6)

expMat_wNuGenes

summary(lm(formula("virulence ~ Rv2351c"), data = expCitMat_d6_scaled))

summary(lm(formula("virulence ~ Rv2350c"), data = expCitMat_d6_scaled))

summary(lm(formula("virulence ~ Rv2349c"), data = expCitMat_d6_scaled))


doRegPlot(lm(formula("virulence ~ Rv2351c"), data = expCitMat_d6_scaled), varY = "virulence", varX = "Rv2351c")

doRegPlot(lm(formula("virulence ~ Rv2350c"), data = expCitMat_d6_scaled), varY = "virulence", varX = "Rv2350c")

doRegPlot(lm(formula("virulence ~ Rv2349c"), data = expCitMat_d6_scaled), varY = "virulence", varX = "Rv2349c")


expCitMat_d6_scaled[order(expCitMat_d6_scaled$virulence), ]
expCitMat_d6_scaled[order(expCitMat_d6_scaled$host_preference), ]

expCitMat_d6_scaled$virulence


# Do added variable plots for the models.

lateApAll_d6 <- doMergedAddVarPlots(lateApFit_d6)

lateApAll_d6

ggsave(paste0(plotDir, "lateAp_d6_addVar.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "lateAp_d6_addVar.pdf"), height = 7.5, width = 6.75)


earlApAll_d6 <- doMergedAddVarPlots(earlApFit_d6)

earlApAll_d6

ggsave(paste0(plotDir, "earlAp_d6_addVar.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "earlAp_d6_addVar.pdf"), height = 7.5, width = 6.75)


viableAll_d6 <- doMergedAddVarPlots(viableFit_d6)

viableAll_d6

ggsave(paste0(plotDir, "viable_d6_addVar.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "viable_d6_addVar.pdf"), height = 7.5, width = 6.75)


necrotAll_d6 <- doMergedAddVarPlots(necrotFit_d6)

necrotAll_d6

ggsave(paste0(plotDir, "necrot_d6_addVar.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "necrot_d6_addVar.pdf"), height = 7.5, width = 6.75)


hostPrAll_d6 <- doMergedAddVarPlots(hostPrFit_d6)

hostPrAll_d6

ggsave(paste0(plotDir, "hostPr_d6_addVar.png"), height = 7.5, width = 6.75)
ggsave(paste0(plotDir, "hostPr_d6_addVar.pdf"), height = 7.5, width = 6.75)
