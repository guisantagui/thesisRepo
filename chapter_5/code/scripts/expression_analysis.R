# Packages
##################################################################################################################

if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(ggpubr)) install.packages("ggpubr")
library(ggpubr)
if(!require(ggpattern))  install_github("coolbutuseless/ggpattern")
library(ggpattern)
if(!require(nnet))  install.packages("nnet")
library(nnet)
if(!require(stringr))  install.packages("stringr")
library(stringr)

# Directory stuff
##################################################################################################################
rootDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcPCRs/"
dataDir <- paste0(rootDir, "data/")
pcrResDir <- paste0(dataDir, "qPCR/")
resDir <- paste0(rootDir, "results/")
if(!dir.exists(resDir)){
        dir.create(resDir)
}
citoDir <- resDir
plotDir <- paste0(rootDir, "plots/")
if(!dir.exists(plotDir)){
        dir.create(plotDir)
}
expressPlotDir <- paste0(plotDir, "expression/")
if(!dir.exists(expressPlotDir)){
        dir.create(expressPlotDir)
}

# Load functions
##################################################################################################################

source(paste0(rootDir, "scripts/expressionAnalFnctns.R"))

# Load and parse the data
##################################################################################################################

# Oligo info dataframe
oligoInfo <- as.data.frame(readxl::read_xlsx(paste0(dataDir, "RNA_info.xlsx")))
colnames(oligoInfo) <- make.names(colnames(oligoInfo))
oligoInfo <- oligoInfo[oligoInfo$RNA.FRACTION == "Prokaryote" | oligoInfo$RNA.FRACTION == "Total", ]
oligoInfo$INFECTION.DATE <- as.Date(oligoInfo$INFECTION.DATE, origin = "1899-12-30")
oligoInfo$SAMPLING.DATE <- as.Date(oligoInfo$SAMPLING.DATE, origin = "1899-12-30")

# plcABCD expression data
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
        scale_color_manual(values = c("green", "blue", "red", "purple", "orange"), 
                           label = c("*plcD* (Mb1784c)",
                                     "*plcC* (Rv2349c)",
                                     "*plcB* (Rv2350c)",
                                     "*plcA* (Rv2351c)", 
                                     "*sigA* (Rv2703)")) + 
        stat_regline_equation(label.x = c(2, 2, 2, 2, 2), 
                              label.y = c(15, 14.2, 13.4, 12.6, 11.8),
                              aes(label =  paste(..eq.label.., ..rr.label.., ..adj.rr.label.., sep = "~~~~"))) +
        labs(y = "Ct", color = "primer pair") +
        xlab(expression(log[10](copies~"in"~PCR))) +
        theme(title = ggtext::element_markdown(),
              axis.text.x = ggtext::element_markdown(),
              axis.title.y = ggtext::element_markdown(),
              legend.text = ggtext::element_markdown(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.grid.major = element_line(colour = "#d4d4d4"))

regPlot

ggsave(paste0(expressPlotDir, "standCurves.pdf"), height = 6, width = 7)
ggsave(paste0(expressPlotDir, "standCurves.png"), height = 6, width = 7)

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

results <- results[results$oligo %in% bestSamples, ]


# Substitute NAs in relative expression by zeros
results$relExp[is.na(results$relExp)] <- 0

# Calculate relative expression normalized to the value of sigA in positive control, to remove any variability 
# between plates
results$relExp_plateNorm <- results$relExp/effs$sigA$ampFact^(-results$CtNorm)


write.csv(results, file = paste0(resDir, "bestSamplesAll.csv"))

# Load preprocessed citometry data
citoNum_d2 <- read.csv(paste0(citoDir, "citoNum_FC_noOLs_d2.csv"), row.names = 1)
citoNum_d6 <- read.csv(paste0(citoDir, "citoNum_FC_noOLs_d6.csv"), row.names = 1)

# Do expression boxplots
##################################################################################################################

plcA_boxplot <- doExpBoxPlot(results, plateNorm = T, gene = "Rv2351c", stats = "t_test", formatTitle = T)
plcB_boxplot <- doExpBoxPlot(results, plateNorm = T, gene = "Rv2350c", stats = "t_test", formatTitle = T)
plcC_boxplot <- doExpBoxPlot(results, plateNorm = T, gene = "Rv2349c", stats = "t_test", formatTitle = T)
plcD_boxplot <- doExpBoxPlot(results, plateNorm = T, gene = "Mb1784c", stats = "t_test", formatTitle = T)

plotAll <- ggarrange(plcA_boxplot$plot,
                     plcB_boxplot$plot,
                     plcC_boxplot$plot,
                     plcD_boxplot$plot,
                     common.legend = T,
                     legend = "bottom")

plotAll

ggsave(paste0(expressPlotDir, "plcExpression_d1.pdf"),
       height = 9.2, width = 9)


expMat <- getExpMat(results, plateNorm = T)

exprMat <- expMat
citoMat <- citoNum_d2

combExprCitoMats <- function(exprMat, citoMat){
        rownames(exprMat) <- gsub("-", "", rownames(exprMat))
        cases <- unique(tolower(str_extract(rownames(exprMat), "[^_]*_[^_]*")))
        
        exprCitoMat <- data.frame(matrix(nrow = 0,
                                         ncol = ncol(exprMat) + ncol(citoMat),
                                         dimnames = list(NULL,
                                                         c(colnames(exprMat),
                                                           colnames(citoMat)))))
        for(c in cases){
                exprToBind <- exprMat[grep(c, tolower(rownames(exprMat))), ]
                toSamp <- citoMat[grep(c, tolower(rownames(citoMat))), ]
                n <- length(grep(c, tolower(rownames(exprMat))))
                if(nrow(toSamp) >= n){
                        replace <- F
                }else{
                        replace <- T
                }
                sampled <- toSamp[sample(1:nrow(toSamp), n, replace = replace), ]
                binded <- cbind.data.frame(exprToBind, sampled)
                exprCitoMat <- rbind.data.frame(exprCitoMat, binded)
        }
        return(exprCitoMat)
}

set.seed(666)
expCit_d2 <-combExprCitoMats(expMat, citoNum_d2)
set.seed(666)
expCit_d6 <-combExprCitoMats(expMat, citoNum_d6)

expCit_d2_std <- apply(expCit_d2, 2, function(x) (x - mean(x))/sd(x))
expCit_d6_std <- apply(expCit_d6, 2, function(x) (x - mean(x))/sd(x))

expCitD2_pca <- pcBiplot(prcomp(expCit_d2_std, scale. = F, center = F))
expCitD2_pca

ggsave(filename = paste0(expressPlotDir, "expCit_d2_pca.pdf"),
       width = 6,
       height = 5)


expCitD6_pca <- pcBiplot(prcomp(expCit_d6_std, scale. = F, center = F))
expCitD6_pca

ggsave(filename = paste0(expressPlotDir, "expCit_d6_pca.pdf"),
       width = 6,
       height = 5)

ggarrange(expCitD2_pca,
          expCitD6_pca,
          common.legend = T,
          legend = "bottom")

ggsave(filename = paste0(expressPlotDir, "expCit_d2d6_pca.pdf"),
       width = 9.5,
       height = 5)

# Get Y axis of the d6 cytometry for obtaining a measure for virulence
virulFormulaPCA <- prcomp(citoNum_d6, scale. = T, center = T)$rotation[, 2]

matrix(virulFormulaPCA, nrow = 1, dimnames = list(NULL, names(virulFormulaPCA)))
virulencePCA_d6 <- expCit_d6_std[, 5:8] %*% matrix(virulFormulaPCA, ncol = 1, dimnames = list(names(virulFormulaPCA), NULL))

sort(virulencePCA_d6[, 1])

max(expCit_d6_std[, "viable"])
min(expCit_d6_std[, "viable"])
max(expCit_d6_std[, "necrotic"])
min(expCit_d6_std[, "necrotic"])

virulenceNecrotViable_d6 <- expCit_d6_std[, c("viable", "necrotic")] %*% matrix(c(-1, 1), ncol = 1)


expCit_d6_4Mods <- expCit_d6_std

virulenceMat_d6 <- cbind(virulencePCA_d6,
                         virulenceNecrotViable_d6)

colnames(virulenceMat_d6) <- c("virulencePCA",
                               "virulenceNecrotViabl")

virulenceMat_d6 <- apply(virulenceMat_d6, 2, function(x) (x - mean(x))/sd(x))

expCit_d6_4Mods <- data.frame(cbind(expCit_d6_4Mods, virulenceMat_d6))

summary(lm(formula("virulenceNecrotViabl ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCit_d6_4Mods))
summary(lm(formula("virulencePCA ~ Rv2349c + Rv2350c + Rv2351c + Mb1784c"), data = expCit_d6_4Mods))



Rv2349c_VN_lm_noBov <- lm(formula("virulenceNecrotViabl ~ Rv2349c"),
                          data = data.frame(apply(expCit_d6_4Mods[!grepl("Bovis", rownames(expCit_d6_4Mods)), ],
                                                  2,
                                                  function(x) (x - mean(x))/sd(x))))
summary(Rv2349c_VN_lm_noBov)

Rv2349c_VN_lm <- lm(formula("virulenceNecrotViabl ~ Rv2349c"),
                    data = expCit_d6_4Mods)
summary(Rv2349c_VN_lm)

Rv2350c_VN_lm <- lm(formula("virulenceNecrotViabl ~ Rv2350c"),
                    data = expCit_d6_4Mods)
summary(Rv2350c_VN_lm)

Rv2351c_VN_lm <- lm(formula("virulenceNecrotViabl ~ Rv2351c"),
                    data = expCit_d6_4Mods)
summary(Rv2351c_VN_lm)

Mb1784c_VN_lm <- lm(formula("virulenceNecrotViabl ~ Mb1784c"),
                    data = expCit_d6_4Mods)
summary(Mb1784c_VN_lm)

summary(lm(formula("virulenceNecrotViabl ~ Rv2349c + Mb1784c"), data = expCit_d6_4Mods))

sumTst <- summary(lm(formula("virulencePCA ~ Rv2349c"), data = expCit_d6_4Mods))


Rv2349c_regPlot_noBov <- doRegPlot(Rv2349c_VN_lm_noBov,
                                   varY = "virulenceNecrotViabl",
                                   varX = "Rv2349c")

Rv2349c_regPlot_noBov

ggsave(paste0(expressPlotDir, "Rv2349c_regPlot_noBov.pdf"), height = 5, width = 6)

Rv2349c_regPlot <- doRegPlot(Rv2349c_VN_lm,
                             varY = "virulenceNecrotViabl",
                             varX = "Rv2349c")

Rv2349c_regPlot

ggsave(paste0(expressPlotDir, "Rv2349c_regPlot.pdf"), height = 5, width = 6)


Rv2350c_regPlot <- doRegPlot(Rv2350c_VN_lm,
                             varY = "virulenceNecrotViabl",
                             varX = "Rv2350c")

Rv2350c_regPlot

ggsave(paste0(expressPlotDir, "Rv2350c_regPlot.pdf"), height = 5, width = 6)


Rv2351c_regPlot <- doRegPlot(Rv2351c_VN_lm,
                             varY = "virulenceNecrotViabl",
                             varX = "Rv2351c")

Rv2351c_regPlot

ggsave(paste0(expressPlotDir, "Rv2351c_regPlot.pdf"), height = 5, width = 6)


Mb1784c_regPlot <- doRegPlot(Mb1784c_VN_lm,
                             varY = "virulenceNecrotViabl",
                             varX = "Mb1784c")

Mb1784c_regPlot

ggsave(paste0(expressPlotDir, "Mb1784c_regPlot.pdf"), height = 5, width = 6)

allUnivRegs_virulence <- ggarrange(Rv2351c_regPlot,
                                   Rv2350c_regPlot,
                                   Rv2349c_regPlot,
                                   Mb1784c_regPlot, 
                                   common.legend = T,
                                   legend = "bottom")

allUnivRegs_virulence

ggsave(paste0(expressPlotDir, "plcABCD_univ_virulence_regPlot.pdf"), height = 10.5, width = 10)

checkExpCitMerge <- function(exprMat,
                             citoMat,
                             nPerm = 2000,
                             alpha = 0.05,
                             predictors = c("Rv2349c", "Rv2350c", "Rv2351c", "Mb1784c")){
        pValsMat <- matrix(nrow = 0, ncol = length(predictors), dimnames = list(NULL, predictors))
        for(i in 1:nPerm){
                exprCitMat <- combExprCitoMats(expMat, citoNum_d6)
                exprCitMat <- apply(exprCitMat, 2, function(x) (x - mean(x))/sd(x))
                virulence <- exprCitMat[, c("viable", "necrotic")] %*% matrix(c(-1, 1), ncol = 1)
                exprCitMat <- cbind(exprCitMat, virulence)
                exprCitMat <- apply(exprCitMat, 2, function(x) (x - mean(x))/sd(x))
                colnames(exprCitMat)[ncol(exprCitMat)] <- "virulence"
                predPVals <- c()
                for(p in predictors){
                        linearMod <- lm(formula(sprintf("virulence ~ %s", p)),
                                        data = data.frame(exprCitMat))
                        fStat <- summary(linearMod)$fstatistic
                        pVal <- pf(fStat[1], fStat[2], fStat[3], lower.tail = F)
                        predPVals <- c(predPVals, pVal)
                }
                names(predPVals) <- predictors
                pVals2Bind <- matrix(predPVals, nrow = 1)
                pValsMat <- rbind(pValsMat, pVals2Bind)
        }
        probs <- apply(pValsMat, 2, function(x) sum(x > alpha)/nrow(pValsMat))
        
        plotList <- list()
        for(p in predictors){
                dataHist <- data.frame(p_value = pValsMat[, p])
                plt <- ggplot(dataHist, mapping = aes(x = p_value)) +
                        geom_histogram(alpha = .6, fill = "black") +
                        geom_vline(xintercept = alpha, col = "red", linetype = "dashed") +
                        xlab("p-value") +
                        theme(axis.title.x = ggtext::element_markdown(),
                              axis.title.y = ggtext::element_markdown(),
                              panel.background = element_blank(),
                              panel.border = element_rect(colour = "black", fill=NA, size=1),
                              panel.grid.major = element_line(colour = "#d4d4d4"),
                              plot.caption=element_text(size = 8.5,
                                                        margin = margin(5,0,0,0)))
                plotList[[p]] <- plt
        }
        out <- list(probs = probs,
                    plots = plotList)
        return(out)
}

chkExpCitRes <- checkExpCitMerge(expMat, citoNum_d6)

chkExpCitRes$probs

chkExpCitRes$plots$Rv2349c
chkExpCitRes$plots$Rv2350c
chkExpCitRes$plots$Rv2351c
chkExpCitRes$plots$Mb1784c

chkExpCitRes_noBov <- checkExpCitMerge(expMat[!grepl("Bovis", rownames(expMat)), ], citoNum_d6)

chkExpCitRes_noBov$probs

chkExpCitRes_noBov$plots$Rv2349c
chkExpCitRes_noBov$plots$Rv2350c
chkExpCitRes_noBov$plots$Rv2351c
chkExpCitRes_noBov$plots$Mb1784c

ggarrange(chkExpCitRes$plots$Rv2349c, chkExpCitRes_noBov$plots$Rv2349c)

ggsave(paste0(expressPlotDir, "pValDistsPlcC.pdf"), width = 8, height = 4)
