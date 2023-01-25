
# Packages
##################################################################################################################

if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(ggsignif)) install.packages("ggsignif")
library(ggsignif)
if(!require(ggpubr)) install.packages("ggpubr")
library(ggpubr)
if(!require(ggpattern))  install_github("coolbutuseless/ggpattern")
library(ggpattern)
if(!require(tidyr)) install.packages("tidyr")
library(tidyr)


# Directory stuff
##################################################################################################################

rootDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcPCRs/"
citoDir <- paste0(rootDir, "data/")
plotDir <- paste0(rootDir, "plots/")
citoPlotDir <- paste0(plotDir, "cytometry/")
if(!dir.exists(citoPlotDir)){
        dir.create(citoPlotDir)
}
resDir <- paste0(rootDir, "results/")

if(!dir.exists(plotDir)){
        dir.create(plotDir)
}
if(!dir.exists(resDir)){
        dir.create(resDir)
}


# Load functions from script
##################################################################################################################
source("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcPCRs/scripts/pcrAnalFnctns.R")

parseRawCito <- function(citoRaw, returnControls = F){
        colnames(citoRaw) <- gsub("muestra",
                                  "samp_withinPlate",
                                  colnames(citoRaw))
        colnames(citoRaw) <- gsub("celula",
                                  "host",
                                  colnames(citoRaw))
        if(!"infeccion_numero" %in% colnames(citoRaw)){
                citoRaw$infeccion_numero <- as.numeric(factor(citoRaw$Infeccion_fecha,
                                                              levels = as.character(sort(unique(citoRaw$Infeccion_fecha)))))
        }
        citoRaw$Bacteria <- sapply(citoRaw$samp_withinPlate, 
                                   function(x) strsplit(x, split = "_")[[1]][1])
        citoRaw$valor <- as.numeric(citoRaw$valor)
        citoRaw$case <- paste(citoRaw$Bacteria, citoRaw$host, sep = "_")
        citoRaw$samp_withinPlate <- sapply(citoRaw$samp_withinPlate, function(x) strsplit(x, split = "_")[[1]][2])
        citoRaw$samp_withinPlate <- paste0("R", citoRaw$samp_withinPlate)
        colnames(citoRaw) <- gsub("samp_withinPlate", "rep_withinPlate", colnames(citoRaw))
        citoPars <- data.frame(matrix(nrow = 0,
                                      ncol = ncol(citoRaw) + 4,
                                      dimnames = list(NULL,
                                                      c(colnames(citoRaw),
                                                        "rep",
                                                        "sample",
                                                        "ctrlUsed",
                                                        "FC_ctrl"))))
        deathOrd <- unique(citoRaw$tipo_muerte)
        citoOrd <- data.frame(matrix(nrow = 0,
                                     ncol = ncol(citoRaw),
                                     dimnames = list(NULL,
                                                     colnames(citoRaw))))
        for(o in unique(paste(citoRaw$infeccion_numero,
                              citoRaw$tiempo,
                              citoRaw$case, 
                              citoRaw$rep_withinPlate,
                              sep = "_"))){
                citoSubOrd <- citoRaw[paste(citoRaw$infeccion_numero,
                                            citoRaw$tiempo,
                                            citoRaw$case, 
                                            citoRaw$rep_withinPlate,
                                            sep = "_") == o, ]
                citoSubOrd <- citoSubOrd[match(deathOrd, citoSubOrd$tipo_muerte), ]
                citoOrd <- rbind.data.frame(citoOrd,
                                            citoSubOrd)
        }
        rownames(citoOrd) <- 1:nrow(citoOrd)
        citoRaw <- citoOrd
        for(n in unique(citoRaw$infeccion_numero)){
                for(d in unique(citoRaw$tiempo)){
                        citoDay <- citoRaw[citoRaw$tiempo == d, ]
                        citoDay_wSampReps <- data.frame(matrix(nrow = 0,
                                                               ncol = ncol(citoDay) + 2,
                                                               dimnames = list(NULL,
                                                                               c(colnames(citoDay),
                                                                                 "rep",
                                                                                 "sample"))))
                        for(c in unique(citoDay$case)){
                                caseCitoDay <- citoDay[citoDay$case == c, ]
                                nReps <- nrow(caseCitoDay)/length(unique(citoRaw$tipo_muerte))
                                repVec <- rep(1:nReps, each = length(unique(citoRaw$tipo_muerte)))
                                repVec <- paste0("R",
                                                 as.character(repVec))
                                caseCitoDay$rep <- repVec
                                caseCitoDay$sample <- paste(caseCitoDay$case,
                                                            repVec,
                                                            sep = "_")
                                citoDay_wSampReps <- rbind.data.frame(citoDay_wSampReps, 
                                                                      caseCitoDay)
                                
                        }
                        for(h in unique(citoDay_wSampReps$host)){
                                subCit <- citoDay_wSampReps[citoDay_wSampReps$infeccion_numero == n & citoDay_wSampReps$tiempo == d & citoDay_wSampReps$host == h, ]
                                conts <- subCit[subCit$Bacteria == "control", ]
                                toDiv <- conts$valor[match(paste0(subCit$tipo_muerte, subCit$rep_withinPlate),
                                                           paste0(conts$tipo_muerte, conts$rep_withinPlate))]
                                ctrlUsed <- conts$sample[match(paste0(subCit$tipo_muerte, subCit$rep_withinPlate),
                                                               paste0(conts$tipo_muerte, conts$rep_withinPlate))]
                                # There are some plates with less controls than samples, so let's assign for these 
                                # samples without control with same number a random control
                                if(any(is.na(toDiv))){
                                        numNAs <- sum(is.na(toDiv))/4
                                        set.seed(123)
                                        randConts <- sample(unique(conts$sample),
                                                            numNAs,
                                                            replace = T)
                                        randContsDF <- conts[conts$sample %in% randConts, ]
                                        toDiv[is.na(toDiv)] <- randContsDF$valor
                                        ctrlUsed[is.na(ctrlUsed)] <- randContsDF$sample
                                        #print(numNAs)
                                        #print(randContsDF)
                                        #print(randContsDF$valor)
                                        #subCit$tipo_muerte[is.na(toDiv)]
                                }
                                subCit$ctrlUsed <- ctrlUsed
                                subCit$FC_ctrl <- subCit$valor / toDiv
                                citoPars <- rbind.data.frame(citoPars, subCit)
                        }
                }
        }        
        citoPars <- citoPars[order(citoPars$infeccion_numero), ]
        citoPars$valor <- as.numeric(citoPars$valor)
        if(returnControls){
                citoPars <- citoPars[citoPars$Bacteria == "control", ]
        }else{
                citoPars <- citoPars[citoPars$Bacteria != "control", ]
        }
        citoPars <- citoPars[!(is.infinite(citoPars$FC_ctrl) | is.na(citoPars$FC_ctrl)), ]
        # Remove samples that don't have values for all citometry deaths
        sampTab <- table(citoPars$sample)
        incomplete <- names(sampTab[sampTab %% length(deathOrd) != 0])
        citoPars <- citoPars[!citoPars$sample %in% incomplete, ]
        rownames(citoPars) <- 1:nrow(citoPars)
        return(citoPars)
}

# Do within or across groups RLA plot
doRLAPlot <- function(featData, whichRLA = NA, normMeth = ""){
        featData <- data.frame(featData)
        stop_quietly <- function() {
                opt <- options(show.error.messages = FALSE)
                on.exit(options(opt))
                stop()
        }
        whichRLA <- tolower(whichRLA)
        if(!(whichRLA == "wg" | whichRLA == "ag")){
                print("Please specify if RLA plot is within-group (WG) or across-group (AG)")
                stop_quietly()
        }
        if(whichRLA == "ag"){
                if(nchar(normMeth) > 0){
                        main <- sprintf("%s, across-group RLA plot", normMeth)
                }else{
                        main <- "Across-group RLA plot"
                }
                featData <- data.frame(apply(featData, 2, function(x) (x - median(x))/sd(x)))
        }
        abuDF <- data.frame(matrix(nrow = 0, 
                                   ncol = 4, 
                                   dimnames = list(NULL,
                                                   c("sample",
                                                     "strain",
                                                     "metabolite",
                                                     "abundance"))))
        for(i in 1:nrow(featData)){
                samp <- rownames(featData)[i]
                sampDF <- data.frame(t(featData[i, ]))
                colnames(sampDF) <- "abundance"
                sampDF$metabolite <- rownames(sampDF)
                sampDF$sample <- rep(samp, nrow(sampDF))
                sampDF$strain <- gsub("\\_.*", "", sampDF$sample)
                sampDF$strain <- str_extract(sampDF$sample, "[^_]*_[^_]*")
                sampDF <- sampDF[, c("sample", "strain", "metabolite", "abundance")]
                abuDF <- rbind.data.frame(abuDF, sampDF)
        }
        rownames(abuDF) <- 1:nrow(abuDF)
        if(whichRLA == "wg"){
                for(strain in unique(abuDF$strain)){
                        strainDF <- abuDF[abuDF$strain == strain, ]
                        for(met in unique(strainDF$metabolite)){
                                abStrainMet <- strainDF$abundance[strainDF$metabolite == met]
                                meanMet <- median(abStrainMet)
                                sdMet <- sd(abStrainMet)
                                abuDF$abundance[abuDF$strain == strain & abuDF$metabolite == met] <- sapply(abuDF$abundance[abuDF$strain == strain & abuDF$metabolite == met],
                                                                                                            function(x) (x - meanMet))
                        }
                }
                if(nchar(normMeth) > 1){
                        main <- sprintf("%s, within-group RLA plot", normMeth)
                }else{
                        main <- "Within-group RLA plot"
                }
        }
        strainCols <- rainbow(length(unique(abuDF$strain)))
        plot <- ggplot(data = abuDF, mapping = aes(x = sample, y = abundance, fill = strain)) +
                geom_boxplot() +
                ylim(c(-4, 4)) +
                scale_discrete_manual(aesthetics = "fill",
                                      values = strainCols) +
                labs(title = main, y = "relative log abundance") +
                theme(title = ggtext::element_markdown(),
                      axis.text.x = element_text(angle = 90),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "none")
        plot
}

# Load and parse the data
##################################################################################################################

infectTableFile <- paste0(citoDir,
                          "tabla_infecciones_completa_updated2.xlsx")

cito_raw_file <- paste0(citoDir,
                        "cito_backup_normalizado_bruto.xlsx")


cito_raw <- data.frame(readxl::read_xlsx(cito_raw_file, sheet = 1))

citoParsed <- parseRawCito(cito_raw)



# Remove contaminated bovis samples samples
citoParsed <- citoParsed[!grepl("contamination", citoParsed$Comments), ]

# Obtain dataframes of day 2 and day 6 post infection
cito_d2 <- citoParsed[citoParsed$tiempo == 2, ]
cito_d6 <- citoParsed[citoParsed$tiempo == 6, ]

# Get dataframes of only controls
controls <- parseRawCito(cito_raw, returnControls = T)
controls_d2 <- controls[controls$tiempo == 2, ]
controls_d6 <- controls[controls$tiempo == 6, ]


table(cito_d2$case)/4
table(cito_d6$case)/4

# Obtain numeric matrix of fold changes, calculed by dividing each sample 
# death value by the one of its respective control (uninfected equivalent
# macrophage cell in the same plate)
citoNum_d2 <- createNumMat(cito_d2, useFCs = T)
citoNum_d6 <- createNumMat(cito_d6, useFCs = T)


contNum_d2 <- createNumMat(controls_d2, useFCs = F)
contNum_d6 <- createNumMat(controls_d6, useFCs = F)


# Detect and remove outliers
##################################################################################################################

# It seems we have some outliers
pcBiplot(prcomp(citoNum_d2, scale. = T, center = T))
pcBiplot(prcomp(citoNum_d6, scale. = T, center = T))

# Detect outliers based on standard deviation method. sdThrshld refers to the standard
# deviation threshold a sample needs to deviate from the mean for bein considered an 
# outlier. If groupWise = T SD and the mean are calculated based on the replicates for
# each case. If it is false they are calculated based in all the samples
detectOLs_SD <- function(citoNum, sdThrshld = 2, groupWise = F){
        if(!require(stringr)) install.packages("stringr")
        library(stringr)
        if(groupWise){
                unCases <- unique(str_extract(rownames(citoNum), "[^_]*_[^_]*"))
                logicDF <- data.frame(matrix(nrow = 0,
                                             ncol = ncol(citoNum),
                                             dimnames = list(NULL,
                                                             colnames(citoNum))))
                for(c in unCases){
                        citCase <- citoNum[grep(c, rownames(citoNum)), ]
                        citCase <- apply(citCase, 2, function(x) (x - mean(x))/sd(x))
                        logCase <- citCase > sdThrshld | citCase < -sdThrshld
                        logicDF <- rbind.data.frame(logicDF, logCase)
                }
        }else{
                scaledCentered <- apply(citoNum, 2, function(x) (x - mean(x))/sd(x))
                logicDF <- scaledCentered > sdThrshld | scaledCentered < -sdThrshld
        }
        OLs <- apply(logicDF, 1, any)
        OLs <- names(OLs)[OLs]
        return(OLs)
}

detectOLs_SD(citoNum_d2, 3, groupWise = T)

# Detect outliers based on Interquartile range (IQR) method. k refers to the factor
# the IQR is multiplied to. groupWise = T refers to the IQR is calculated based 
# on the replicates. If it is false IQR is calculated based in all the samples
detectOLs_IQR <- function(citoNum, k = 1.5, groupWise = F, logDF = F){
        if(!require(stringr)) install.packages("stringr")
        library(stringr)
        getLogic <- function(citoDF, k){
                logicFun <- data.frame(matrix(nrow = nrow(citoDF),
                                              ncol = 0,
                                              dimnames = list(rownames(citoDF),
                                                              NULL)))
                for(i in 1:ncol(citoDF)){
                        name <- colnames(citoDF)[i]
                        quants <- quantile(citoDF[, i], probs = c(.75, .25))
                        iqr <- quants[1] - quants[2]
                        cutoff <- iqr * k
                        interval <- c(quants[2] - cutoff,
                                      quants[1] + cutoff)
                        log2Bind <- data.frame(matrix(citoDF[, i] < interval[1] | citoDF[, i] > interval[2],
                                                      nrow = nrow(citoDF),
                                                      ncol = 1,
                                                      dimnames = list(rownames(citoDF),
                                                                      name)))
                        logicFun <- cbind.data.frame(logicFun, log2Bind)
                }
                return(logicFun)
        }
        
        if(groupWise){
                logicDF <- data.frame(matrix(nrow = 0,
                                             ncol = ncol(citoNum),
                                             dimnames = list(NULL,
                                                             colnames(citoNum))))
                unCases <- unique(str_extract(rownames(citoNum), "[^_]*_[^_]*"))
                for(c in unCases){
                        citCase <- citoNum[grep(c, rownames(citoNum)), ]
                        subLog <- getLogic(citCase, k)
                        logicDF <- rbind.data.frame(logicDF,
                                                    subLog)
                }
        }else{
                logicDF <- getLogic(citoNum, k)
        }
        if(logDF){
                print(logicDF)
        }
        OLs <- apply(logicDF, 1, any)
        OLs <- names(OLs)[OLs]
        return(OLs)
}

detectOLs_IQR(citoNum_d2, groupWise = T)

# Remove outliers iteratively until all the samples are within the requirements.
# Maybe too conservative? 
removeOLs <- function(citoNum, sdThrshld, groupWise = F, OLmthd = "sd"){
        if(OLmthd == "sd"){
                OLFunc <- detectOLs_SD
        }else if(OLmthd == "iqr"){
                OLFunc <- detectOLs_IQR
        }
        noOLs <- citoNum
        OLs <- OLFunc(noOLs, sdThrshld)
        while(length(OLs) > 0){
                noOLs <- noOLs[!rownames(noOLs) %in% OLs, ]
                OLs <- OLFunc(noOLs, sdThrshld)
        }
        return(noOLs)
}

removeOLs(citoNum_d2, 3)

#citoNum_d2 <- removeOLs(citoNum_d2, 3)
#citoNum_d6 <- removeOLs(citoNum_d6, 3)

#citoNum_d2 <- citoNum_d2[!rownames(citoNum_d2) %in% detectOLs_SD(citoNum_d2, sdThrshld = 2, groupWise = T), ]
citoNum_d2_noOLs <- citoNum_d2[!rownames(citoNum_d2) %in% detectOLs_IQR(citoNum_d2, 1.5, groupWise = T), ]
citoNum_d6_noOLs <- citoNum_d6[!rownames(citoNum_d6) %in% detectOLs_IQR(citoNum_d6, 1.5, groupWise = T), ]

write.csv(citoNum_d2_noOLs, file = paste0(resDir, "citoNum_FC_noOLs_d2.csv"))
write.csv(citoNum_d6_noOLs, file = paste0(resDir, "citoNum_FC_noOLs_d6.csv"))

# Plot the raw values of viability after 6 days of infection, along controls, to see how virulent are the 
# strains

# Obtain unnormalized citonum
citoNum_d2_noFC <-  createNumMat(cito_d2, useFCs = F)
citoNum_d6_noFC <-  createNumMat(cito_d6, useFCs = F)

citoNum_d2_noFC[rownames(citoNum_d2_noFC) %in% rownames(citoNum_d2_noOLs), ]

citoNum_d6_noFC[rownames(citoNum_d6_noFC) %in% rownames(citoNum_d6_noOLs) & grepl("bomac", rownames(citoNum_d6_noFC)), ]
citoNum_d6_noFC[rownames(citoNum_d6_noFC) %in% rownames(citoNum_d6_noOLs) & grepl("thp1", rownames(citoNum_d6_noFC)), ]

doRawCitoBoxplots <- function(rawCitoNum, ctrlsNum, deathCat, daysPI){
        parsed <- citoParsed[citoParsed$tiempo == daysPI, ]
        ctrlsNum <- ctrlsNum[unique(parsed$ctrlUsed[match(rownames(rawCitoNum), parsed$sample)]), ]
  
        rawCitoNum <- rbind.data.frame(rawCitoNum, ctrlsNum)
        print(rawCitoNum)
        
        rawPlotDF <- data.frame(sample = rownames(rawCitoNum),
                                lin = sapply(rownames(rawCitoNum), function(x) strsplit(x, split = "_")[[1]][1]),
                                host = sapply(rownames(rawCitoNum), function(x) strsplit(x, split = "_")[[1]][2]),
                                value = rawCitoNum[, deathCat],
                                stringsAsFactors = F)
  
        rawPlotDF$host <- gsub("thp1", "THP1", rawPlotDF$host)
        rawPlotDF$host <- gsub("bomac", "BoMac", rawPlotDF$host)
  
  
        rawPlotDF$lin <- gsub("chimp", "Chimp", rawPlotDF$lin)
        rawPlotDF$lin <- gsub("bovis", "*M. bovis*", rawPlotDF$lin)
        rawPlotDF$lin <- gsub("control", "Uninfected", rawPlotDF$lin)
        rawPlotDF$lin <- factor(rawPlotDF$lin,
                                levels = c("Uninfected",
                                           "*M. bovis*",
                                           "Chimp",
                                           "L6",
                                           "L5",
                                           "ctrl")[c("Uninfected",
                                                     "*M. bovis*",
                                                     "Chimp",
                                                     "L6",
                                                     "L5",
                                                     "ctrl") %in% rawPlotDF$lin])
  
        deathCat <- gsub("viable", "viability", deathCat)
        deathCat <- gsub("necrotic", "necrosis", deathCat)
        deathCat <- gsub("apoptotic", "apoptosis", deathCat)
  
        deathCat <- paste0(toupper(substr(deathCat, 1, 1)),
                           substr(deathCat, 2, nchar(deathCat)))
        deathCat <- gsub("_", " ", deathCat)
  
        yLab <- paste(deathCat, "percentage", sep = " ")
  
        rawPlotDF$host <- factor(rawPlotDF$host,
                                 levels = sort(unique(rawPlotDF$host),
                                               decreasing = T))
  
        p <- ggplot(rawPlotDF, aes(y = value, x = lin, fill = lin)) + 
                geom_boxplot_pattern(pattern_color = "black",
                                     pattern_fill = "black",
                                     pattern_spacing = 0.015,
                                     aes(pattern = host),
                                     outlier.shape = NA) +
                scale_color_manual(values = c("#000000", "#010101")) +
                scale_pattern_manual(values = c("stripe", "none")) +
                labs(x = "Strain",
                     y = yLab) +
                scale_fill_manual(values = c("#808080", "#f279ce", "#c4bf62", "#24ad37", "#871414"), 
                                  labels = c("Uninfected", "*M. bovis*", "Chimp", "L6", "L5")) +
                geom_point(aes(col = host), position=position_jitterdodge(jitter.width = 0.005)) +
                theme(title = ggtext::element_markdown(),
                      axis.title.x = ggtext::element_markdown(),
                      axis.text.x = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      legend.text = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"))
        p
}

ggarrange(doRawCitoBoxplots(citoNum_d2_noFC[rownames(citoNum_d2_noFC) %in% rownames(citoNum_d2_noOLs), ],
                            contNum_d2, 
                            deathCat = "viable",
                            daysPI = 2),
          doRawCitoBoxplots(citoNum_d2_noFC[rownames(citoNum_d2_noFC) %in% rownames(citoNum_d2_noOLs), ],
                            contNum_d2, 
                            deathCat = "necrotic",
                            daysPI = 2),
          doRawCitoBoxplots(citoNum_d2_noFC[rownames(citoNum_d2_noFC) %in% rownames(citoNum_d2_noOLs), ],
                            contNum_d2, 
                            deathCat = "early_apoptotic",
                            daysPI = 2),
          doRawCitoBoxplots(citoNum_d2_noFC[rownames(citoNum_d2_noFC) %in% rownames(citoNum_d2_noOLs), ],
                            contNum_d6, 
                            deathCat = "late_apoptotic",
                            daysPI = 2),
          common.legend = T,
          legend = "bottom")

ggsave(paste0(citoPlotDir, "rawCitoBxpltsD2.pdf"),
       height = 9, width = 9)


ggarrange(doRawCitoBoxplots(citoNum_d6_noFC[rownames(citoNum_d6_noFC) %in% rownames(citoNum_d6_noOLs), ],
                            contNum_d6, 
                            deathCat = "viable",
                            daysPI = 6),
          doRawCitoBoxplots(citoNum_d6_noFC[rownames(citoNum_d6_noFC) %in% rownames(citoNum_d6_noOLs), ],
                            contNum_d6, 
                            deathCat = "necrotic",
                            daysPI = 6),
          doRawCitoBoxplots(citoNum_d6_noFC[rownames(citoNum_d6_noFC) %in% rownames(citoNum_d6_noOLs), ],
                            contNum_d6, 
                            deathCat = "early_apoptotic",
                            daysPI = 6),
          doRawCitoBoxplots(citoNum_d6_noFC[rownames(citoNum_d6_noFC) %in% rownames(citoNum_d6_noOLs), ],
                            contNum_d6, 
                            deathCat = "late_apoptotic",
                            daysPI = 6),
          common.legend = T,
          legend = "bottom")

ggsave(paste0(citoPlotDir, "rawCitoBxpltsD6.pdf"),
       height = 9, width = 9)

# day 6 has many zeros in necrosis, which appear in raw data--> suspicious, let's do biplots coloring according
# to date to see if there is something odd. These are the infections of 2021-08-03

# Coloring by strain
ggarrange(pcBiplot(prcomp(citoNum_d2_noOLs, scale. = T, center = T)),
          pcBiplot(prcomp(citoNum_d6_noOLs, scale. = T, center = T)),
          legend = "bottom",
          common.legend = T)

ggsave(paste0(citoPlotDir, "FC_PCA_noOLs.pdf"),
       height = 4, width = 7.5)

# Coloring by date
ggarrange(pcBiplot_date(prcomp(citoNum_d2_noOLs, scale. = T, center = T),
                        sampInfo = citoParsed, linNames = T),
          pcBiplot_date(prcomp(citoNum_d6_noOLs, scale. = T, center = T),
                        sampInfo = citoParsed, linNames = T),
          legend = "right",
          common.legend = T)
ggsave(paste0(citoPlotDir, "FC_PCA_noOLs_wDate.pdf"),
       height = 4, width = 9)

# Determine if there are differences between uninfected bomac and thp1 and if we're suffering from batch effect
##################################################################################################################

# Do univariate tests between bomac and thp1 uninfected cells to see if they
# have differences in the cell death phenotypes

# Do Shapiro test to see if we can do parametric tests
doShapsCitos <- function(citoNum){
        hosts <- sort(unique(sapply(rownames(citoNum),
                                    function(x) strsplit(x, split = "_")[[1]][2])))
        pValsDF <- data.frame(matrix(nrow = 0,
                                     ncol = ncol(citoNum),
                                     dimnames = list(NULL,
                                                     colnames(citoNum))))
        for(h in hosts){
                pVals <- apply(citoNum,
                               2,
                               function(x) shapiro.test(x[grep(h, rownames(citoNum))])$p.value)
                pVals <- data.frame(matrix(pVals,
                                           nrow = 1,
                                           ncol = length(pVals),
                                           dimnames = list(h,
                                                           names(pVals))))
                pValsDF <- rbind.data.frame(pValsDF, pVals)
        }
        return(pValsDF)
}

# This function does mann whitney test between bomac and thp1 cells
doMannWhit <- function(citoNum, whatComp = "hosts"){
        hosts <- sort(unique(sapply(rownames(citoNum),
                                    function(x) strsplit(x, split = "_")[[1]][2])))
        lins <- sort(unique(sapply(rownames(citoNum),
                                   function(x) strsplit(x, split = "_")[[1]][1])))
        pValsDF <- data.frame(matrix(nrow = length(lins),
                                     ncol = 0,
                                     dimnames = list(lins,
                                                     NULL)))
        for(i in 1:ncol(citoNum)){
                deathCat <- colnames(citoNum)[i]
                death <- citoNum[, i]
                pVals <- sapply(lins,
                                function(x) wilcox.test(death[grepl(x,
                                                                    rownames(citoNum)) & 
                                                                      grepl(hosts[1],
                                                                            rownames(citoNum))],
                                                        death[grepl(x,
                                                                    rownames(citoNum))
                                                              & grepl(hosts[2],
                                                                      rownames(citoNum))])$p.value)
                pVals <- data.frame(matrix(pVals,
                                           ncol = 1,
                                           nrow = length(lins),
                                           dimnames = list(names(pVals),
                                                           deathCat)))
                pValsDF <- cbind.data.frame(pValsDF, pVals)
        }
        return(pValsDF)
}

doMannWhit <- function(citoNum, whatComp = "hosts"){
        hosts <- sort(unique(sapply(rownames(citoNum),
                                    function(x) strsplit(x, split = "_")[[1]][2])))
        lins <- sort(unique(sapply(rownames(citoNum),
                                   function(x) strsplit(x, split = "_")[[1]][1])))
        pValsDF <- data.frame(matrix(nrow = length(lins),
                                     ncol = 0,
                                     dimnames = list(lins,
                                                     NULL)))
        outList <- list()
        for(i in 1:ncol(citoNum)){
                deathCat <- colnames(citoNum)[i]
                death <- citoNum[, i]
                names(death) <- rownames(citoNum)
                if(whatComp == "hosts"){
                        pVals <- sapply(lins,
                                        function(x) wilcox.test(death[grepl(x,
                                                                            rownames(citoNum)) & 
                                                                              grepl(hosts[1],
                                                                                    rownames(citoNum))],
                                                                death[grepl(x,
                                                                            rownames(citoNum))
                                                                      & grepl(hosts[2],
                                                                              rownames(citoNum))])$p.value)
                        pVals <- data.frame(matrix(pVals,
                                                   ncol = 1,
                                                   nrow = length(lins),
                                                   dimnames = list(names(pVals),
                                                                   deathCat)))
                        pValsDF <- cbind.data.frame(pValsDF, pVals)
                }else if(whatComp == "strains"){
                        hostListPVals <- list()
                        for(h in hosts){
                                deathHost <- death[grep(h, names(death))]
                                deathHostPVals <- data.frame(matrix(nrow = length(lins) - 1,
                                                                    ncol = length(lins) - 1,
                                                                    dimnames = list(lins[1:(length(lins) - 1)],
                                                                                    lins[2:length(lins)])))
                                for(j in 1:(length(lins) - 1)){
                                        strainA <- lins[j]
                                        for(k in (j + 1):length(lins)){
                                                strainB <- lins[k]
                                                pVal <- wilcox.test(deathHost[grep(strainA, names(deathHost))],
                                                                    deathHost[grep(strainB, names(deathHost))])$p.value
                                                if(deathCat == "necrotic" & h == "thp1" & strainA == "chimp" & strainB == "L6"){
                                                        print(deathCat)
                                                        print(h)
                                                        print(strainA)
                                                        print(strainB)
                                                        print(pVal)
                                                        print(deathHost[grep(strainA, names(deathHost))])
                                                        print(deathHost[grep(strainB, names(deathHost))])
                                                        print(wilcox.test(deathHost[grep(strainA, names(deathHost))],
                                                                          deathHost[grep(strainB, names(deathHost))]))
                                                }
                                                deathHostPVals[strainA, strainB] <- pVal
                                        }
                                }
                                hostListPVals[[h]] <- deathHostPVals
                        }
                        outList[[deathCat]] <- hostListPVals
                }
        }
        if(whatComp == "hosts"){
                return(pValsDF)
        }else if(whatComp == "strains"){
                return(outList)
        }
}
for(i in 1:(4-1)){
        for(j in (i + 1):4){
                print(paste(i, j, sep = "_"))
        }
}




doShapsCitos(contNum_d2[!rownames(contNum_d2) %in% detectOLs_IQR(contNum_d2, groupWise = T), ])
doShapsCitos(contNum_d6[!rownames(contNum_d6) %in% detectOLs_IQR(contNum_d6, groupWise = T), ])

doShapsCitos(contNum_d2)
doShapsCitos(contNum_d6)
# In general it seems that data doesn't follow a normal distribution

# Do Mann Whitney U test between BOMAC and THP-1 for each type of death to
# see if they have significative differences in any of the types of cell death.
apply(contNum_d2,
      2,
      function(x) wilcox.test(x[grep("thp1",
                                     rownames(contNum_d2))],
                              x[grep("bomac",
                                     rownames(contNum_d2))])$p.value)
# At day 2 there are no significative differences in any one of the types of
# cell death

apply(contNum_d6,
      2,
      function(x) wilcox.test(x[grep("thp1",
                                     rownames(contNum_d6))],
                              x[grep("bomac",
                                     rownames(contNum_d6))])$p.value)

# In day 6 there are significative differences in late apoptosis and viability


# Obtain boxplots of the types of death of bomac and thp1 controls, to see if there are differences

# Obtain boxplots of the types of death of bomac and thp1 controls, to see if there are differences
doCitoBoxplots <- function(citoNum, xAxis = "", showDate = F, sampInfo = NULL, daysPI = NULL, deathCat = NULL, plotStats = F){
        scaleFUN = function(x) sprintf("%.1f", x)
        citoNum <- data.frame(citoNum, stringsAsFactors = F)
        if(plotStats){
                statsDF <- doMannWhit(citoNum)
                rownames(statsDF) <- gsub("bovis", "*M. bovis*", rownames(statsDF))
                rownames(statsDF) <- gsub("chimp", "Chimp. B.", rownames(statsDF))
                statSymbDF <- data.frame(value = c(1, 0.1, 0.05, 0.01, 0.001),
                                         symbol = c("ns", ".", "*", "**", "***"))
        }
        stop_quietly <- function(){
                opt <- options(show.error.messages = FALSE)
                on.exit(options(opt))
                stop()
        }
        if(xAxis == "strain"){
                showDate <- F
                sampInfo <- NULL
                daysPI <- NULL
                if(is.null(deathCat)){
                        print("Introduce the death category")
                        stop_quietly()
                }
        }else{
                deathCat <- NULL
        }
        samps <- rownames(citoNum)
        dataPlot <- data.frame(matrix(nrow = 0,
                                      ncol = 4,
                                      dimnames = list(NULL,
                                                      c("sample",
                                                        "host",
                                                        "type_of_death",
                                                        "value"))))
        for(s in samps){
                toBind <- data.frame(sample = rep(s, ncol(citoNum)),
                                     stringsAsFactors = F)
                
                toBind$host <- sapply(toBind$sample,
                                      function(x) strsplit(x, split = "_")[[1]][2])
                toBind$host <- gsub("thp1", "THP1", toBind$host)
                toBind$host <- gsub("bomac", "BoMac", toBind$host)
                subDeath <- citoNum[rownames(citoNum) == s, ]
                toBind$type_of_death <- colnames(subDeath)
                toBind$hostDeath <- paste(toBind$host,
                                          toBind$type_of_death, 
                                          sep = "_")
                toBind$value <- as.numeric(subDeath)
                dataPlot <- rbind.data.frame(dataPlot, toBind)
        }
        dataPlot$strain <- sapply(dataPlot$sample,
                                  function(x) strsplit(x, split = "_")[[1]][1])
        dataPlot$strain <- gsub("chimp", "Chimp. B.", dataPlot$strain)
        dataPlot$strain <- gsub("bovis", "*M. bovis*", dataPlot$strain)
        dataPlot$strain <- factor(dataPlot$strain,
                                  levels = c("*M. bovis*",
                                             "Chimp. B.",
                                             "L6",
                                             "L5",
                                             "ctrl")[c("*M. bovis*",
                                                       "Chimp. B.",
                                                       "L6",
                                                       "L5",
                                                       "ctrl") %in% dataPlot$strain])
        
        dataPlot$host <- factor(dataPlot$host,
                                levels = sort(unique(dataPlot$host),
                                              decreasing = T))
        if(!is.null(deathCat)){
                dataPlot <- dataPlot[dataPlot$type_of_death == deathCat, ]
        }
        if(showDate == T){
                if(is.null(sampInfo) | is.null(daysPI)){
                        print("For showing dates it's necessary to provide data info datframe and days post infection")
                        stop_quietly()
                }else{
                        sampInfo <- sampInfo[sampInfo$tiempo == daysPI, ]
                        dataPlot$date <- as.character(sampInfo$Infeccion_fecha[match(dataPlot$sample, sampInfo$ctrlUsed)])
                        plot <- ggplot(dataPlot, aes(y = value, x = type_of_death, col = date)) +
                                geom_boxplot_pattern(color = "gray",
                                                     pattern_color = "gray",
                                                     pattern_fill = "gray",
                                                     pattern_spacing = 0.015,
                                                     aes(pattern = host),
                                                     outlier.shape = NA) + 
                                scale_pattern_manual(values = c("stripe", "none")) +
                                geom_point(aes(shape = host), position=position_jitterdodge()) +
                                labs(x = "death category",
                                     y = "percentage") +
                          scale_y_continuous(labels = scaleFUN)
                }
        }else if(xAxis == "strain"){
                #dataPlot$type_of_death <- gsub("viable", "viability",
                #                               dataPlot$type_of_death)
                #dataPlot$type_of_death <- gsub("necrotic", "necrosis",
                #                               dataPlot$type_of_death)
                yLab <- unique(dataPlot$type_of_death)
                yLab <- gsub("viable", "viability",
                             yLab)
                yLab <- gsub("necrotic", "necrosis",
                             yLab)
                yLab <- gsub("apoptotic", "apoptosis",
                             yLab)
                yLab <- gsub("_", " ", yLab)
                yLab <- paste(yLab, "fold change", sep = " ")
                plot <- ggplot(dataPlot, aes_string(y = "value", x = xAxis, fill = xAxis)) +
                        geom_boxplot_pattern(pattern_color = "black",
                                             pattern_fill = "black",
                                             pattern_spacing = 0.015,
                                             aes(pattern = host),
                                             outlier.shape = NA) + 
                        scale_color_manual(values = c("#000000", "#010101")) +
                        scale_pattern_manual(values = c("stripe", "none")) +
                        labs(x = "Strain",
                             y = yLab) +
                        scale_fill_manual(values = c("#f279ce", "#c4bf62", "#24ad37", "#871414"), 
                                          labels = c("*M. bovis*", "Chimp. B.", "L6", "L5")) +
                        geom_point(aes(col = host), position=position_jitterdodge(jitter.width = 0.005)) +
                        scale_y_continuous(labels = scaleFUN)
                
        }
        else{
                plot <- ggplot(dataPlot, aes(y = value, x = type_of_death)) +
                        scale_color_manual(values = c("#000000", "#010101")) +
                        geom_boxplot_pattern(color = "gray",
                                             pattern_color = "gray",
                                             pattern_fill = "gray",
                                             pattern_spacing = 0.015,
                                             aes(pattern = host),
                                             outlier.shape = NA) + 
                        scale_pattern_manual(values = c("none", "stripe")) +
                        geom_point(aes(col = host), position=position_jitterdodge()) +
                        labs(x = "death category",
                             y = "percentage")  +
                        scale_y_continuous(labels = scaleFUN)
        }
        if(plotStats){
                dataPlot$type_of_death <- factor(dataPlot$type_of_death)
                # Change order of stats DF to make them fit to factors
                
                trueDimChange <- unlist(lapply(dimnames(statsDF),
                                               function (x) all(x %in% levels(dataPlot[, colnames(dataPlot) == xAxis]))))
                
                ordVec <- match(levels(dataPlot[, colnames(dataPlot) == xAxis]),
                                dimnames(statsDF)[trueDimChange][[1]])
                if(trueDimChange[1]){
                        statsDF <- statsDF[ordVec, ]
                }else{
                        statsDF <- statsDF[, ordVec]
                }
                # Obtain the y position of significance bars, proportional to maximum value of each
                # x category
                maxYs <- sapply(unique(dataPlot[, colnames(dataPlot) == xAxis]),
                                function(x) max(dataPlot$value[dataPlot[, colnames(dataPlot) == xAxis] == x]))
                
                #maxYs <- sapply(levels(dataPlot[, colnames(dataPlot) == xAxis]),
                #                function(x) max(dataPlot$value[dataPlot[, colnames(dataPlot) == xAxis] == x]))
                names(maxYs) <- unique(dataPlot[, colnames(dataPlot) == xAxis])
                maxYs <- maxYs[match(levels(dataPlot[, colnames(dataPlot) == xAxis]),
                                     names(maxYs))]
                
                # Map the pValues of the stat dataframe to the appropriate symbol
                pVals <- statsDF[, colnames(statsDF) %in% unique(dataPlot$type_of_death)]
                pVals <- as.numeric(pVals)
                names(pVals) <- rownames(statsDF)
                names(pVals) <- levels(dataPlot[, colnames(dataPlot) == xAxis])
                pVals <- pVals[match(levels(dataPlot[, colnames(dataPlot) == xAxis]),
                                     names(pVals))]
                
                
                symbols <- sapply(pVals, function(x) statSymbDF$symbol[max(which(x < statSymbDF$value))])
                maxYs <- maxYs + max(maxYs) * 0.025
                signifDF <- data.frame(x = 1:(length(levels(dataPlot[, colnames(dataPlot) == xAxis]))) - 0.2,
                                       xend = 1:(length(levels(dataPlot[, colnames(dataPlot) == xAxis]))) + 0.2,
                                       y = maxYs,
                                       annotation = symbols,
                                       group = 1:length(maxYs))
                plot <- plot +
                        geom_signif(stat = "identity",
                                    data = signifDF,
                                    aes(x=x, xend=xend, y=y, yend=y, fill = NULL, annotation = annotation, group = group, col = NULL))
        }
        plot <- plot +
                theme(title = ggtext::element_markdown(),
                      axis.title.x = ggtext::element_markdown(),
                      axis.text.x = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      legend.text = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"))
        
        if(plotStats){
                out <- list(stats = statsDF,
                            plot = plot)
                return(out)
        }else{
                return(plot)
        }
}

ctrlsBxplt_d2 <- doCitoBoxplots(contNum_d2, xAxis = "type_of_death", plotStats = T)
ctrlsBxplt_d6 <- doCitoBoxplots(contNum_d6, xAxis = "type_of_death", plotStats = T)

ggarrange(ctrlsBxplt_d2$plot,
          ctrlsBxplt_d6$plot +
                  rremove("ylab"),
          common.legend = T,
          legend = "bottom")

ggsave(paste0(citoPlotDir, "ctrls_bxplts.pdf"),
       height = 4.5, width = 9)

ctrlsBxplt_d2_wDate <- doCitoBoxplots(contNum_d2,
                                      showDate = T,
                                      citoParsed,
                                      daysPI = 2,
                                      xAxis = "type_of_death",
                                      plotStats = T)

ctrlsBxplt_d6_wDate <- doCitoBoxplots(contNum_d6,
                                      showDate = T,
                                      citoParsed,
                                      daysPI = 6,
                                      xAxis = "type_of_death",
                                      plotStats = T)

ggarrange(ctrlsBxplt_d2_wDate$plot,
          ctrlsBxplt_d6_wDate$plot +
                  rremove("ylab"),
          common.legend = T,
          legend = "right")

ggsave(paste0(citoPlotDir, "ctrls_bxplts_wDate.pdf"),
       height = 4.5, width = 10)

# The controls of the older infections (2020-12-01 and 2020-12-15) have higher values 
# of viablility, while there have lower values of early and late apoptosis. There is 
# Something odd about the controls of 2021-08-03.


# Obtain biplots of control samples to asess if BOMAC and THP-1 have intrinsic differences in the 
# cell death profiles, and to determine if there is batch variation
ggarrange(pcBiplot(prcomp(contNum_d2, scale. = T, center = T)),
          pcBiplot(prcomp(contNum_d6, scale. = T, center = T)),
          legend = "bottom",
          common.legend = T)

ggsave(paste0(citoPlotDir, "ctrls_PCA.pdf"),
       height = 4, width = 7.5)

ggarrange(pcBiplot_date(prcomp(contNum_d2, scale. = T, center = T),
                        sampInfo = controls_d2),
          pcBiplot_date(prcomp(contNum_d6, scale. = T, center = T),
                        sampInfo = controls_d6),
          legend = "right",
          common.legend = T)

ggsave(paste0(citoPlotDir, "ctrlsPCA_wDate.pdf"),
       height = 4, width = 9)
# It seems that there is some separation by date, being older controls separated from the other ones,
# specially infections from 2020-12-01 and 2020-12-15, which have higher viability.



# Do PCAs of unnormalized data along controls, coloring by date, to see if there is any 
# trend related to when the experiments were ran. 

# Do PCAs of unnormalized data
ggarrange(pcBiplot_date(prcomp(citoNum_d2_noFC, scale. = T, center = T),
                        sampInfo = citoParsed, linNames = T),
          pcBiplot_date(prcomp(citoNum_d6_noFC, scale. = T, center = T),
                        sampInfo = citoParsed, linNames = T),
          common.legend = T,
          legend = "right")

ggsave(paste0(citoPlotDir, "unNormPCA_wDate.pdf"),
       height = 4, width = 9)

# Do PCAs of unnormalized data along controls, coloring by date, to see if there is any pattern
ggarrange(pcBiplot_date(prcomp(rbind.data.frame(citoNum_d2_noFC, contNum_d2), scale. = T, center = T),
                        sampInfo = rbind.data.frame(citoParsed, controls), linNames = T),
          pcBiplot_date(prcomp(rbind.data.frame(citoNum_d6_noFC, contNum_d6), scale. = T, center = T),
                        sampInfo = rbind.data.frame(citoParsed, controls), linNames = T),
          common.legend = T,
          legend = "right")

ggsave(paste0(citoPlotDir, "unNormCtrlsPCA_wDate.pdf"),
       height = 4, width = 9)

ggsave(paste0(citoPlotDir, "unNormCtrlsPCA_wDate.png"),
       height = 4, width = 9)

ggarrange(pcBiplot(prcomp(rbind.data.frame(citoNum_d2_noFC, contNum_d2), scale. = T, center = T)),
          pcBiplot(prcomp(rbind.data.frame(citoNum_d6_noFC, contNum_d6), scale. = T, center = T)),
          common.legend = T,
          legend = "right")

ggsave(paste0(citoPlotDir, "unNormCtrlsPCA.pdf"),
       height = 4, width = 9)

ggsave(paste0(citoPlotDir, "unNormCtrlsPCA.png"),
       height = 4, width = 9)

# These PCAs show that the infections are aggregated according to date and that they cluster along 
# their controls, so this aggregation is not due to in the same date were analyzed strains with 
# similar behaviour, so as a consequence biological and experimental variation overlap. The fact
# that controls aggregate with the samples they were analyzed on the same date instead of with the 
# other controls suggests that WE ARE SUFFERING FROM BATCH EFFECT.


# Remove batch effect with fancy normalization methods
##################################################################################################################

# Let's yry to use some NOMISish and CCMNish methods for normalizing the data, using the controls 
# (BOMAC and THP1) between runs. We are gonna add columns with the cell death values of the 
# corresponding control to each row. We will use these values as internal standards

# Remove outliers with IQR
#citoNum_d2_noFC <- citoNum_d2_noFC[!rownames(citoNum_d2_noFC) %in% detectOLs_IQR(citoNum_d2_noFC, groupWise = T), ]
#citoNum_d6_noFC <- citoNum_d6_noFC[!rownames(citoNum_d6_noFC) %in% detectOLs_IQR(citoNum_d6_noFC, groupWise = T), ]

# Add the corresponding controls to each row
citoNum_d2_noFC_wConts <- cbind.data.frame(citoNum_d2_noFC,
                                           contNum_d2[citoParsed$ctrlUsed[match(rownames(citoNum_d2_noFC),
                                                                                citoParsed$sample)], ])

colnames(citoNum_d2_noFC_wConts)[5:8] <- paste(colnames(citoNum_d2_noFC_wConts)[5:8],
                                               "ctrl",
                                               sep = "_")

citoNum_d6_noFC_wConts <- cbind.data.frame(citoNum_d6_noFC,
                                           contNum_d6[citoParsed$ctrlUsed[match(rownames(citoNum_d6_noFC),
                                                                                citoParsed$sample)], ])

colnames(citoNum_d6_noFC_wConts)[5:8] <- paste(colnames(citoNum_d6_noFC_wConts)[5:8],
                                               "ctrl",
                                               sep = "_")

#citoNum_d2_noFC_wConts <- citoNum_d2_noFC_wConts[rownames(citoNum_d2_noFC_wConts) %in% detectOLs_IQR(citoNum_d2_noFC_wConts, groupWise = T), ]
#citoNum_d6_noFC_wConts <- citoNum_d6_noFC_wConts[rownames(citoNum_d6_noFC_wConts) %in% detectOLs_IQR(citoNum_d6_noFC_wConts, groupWise = T), ]

detectOLs_IQR(citoNum_d6_noFC_wConts, groupWise = T)
detectOLs_IQR(citoNum_d6, groupWise = T)

detectOLs_IQR(citoNum_d6, groupWise = T) %in% detectOLs_IQR(citoNum_d6_noFC_wConts, groupWise = T)



# NOMIS
citoNomis <- function(citoNum_ctrls){
        isMat <- citoNum_ctrls[, grepl("ctrl", colnames(citoNum_ctrls))]
        anMat <- citoNum_ctrls[, !grepl("ctrl", colnames(citoNum_ctrls))]
        isMat <- cbind(rep(1, nrow(isMat)),
                       isMat)
        isMat <- as.matrix(isMat)
        anMat <- as.matrix(anMat)
        beta <- solve(t(isMat) %*% isMat) %*% (t(isMat) %*% anMat)
        #print(beta)
        anMatNorm <- anMat - isMat %*% beta + t(replicate(nrow(anMat), 
                                                          apply(anMat, 2, mean)))
        return(anMatNorm)
}


citoNum_d2_nomis <- citoNomis(citoNum_d2_noFC_wConts[!rownames(citoNum_d2_noFC_wConts) %in% detectOLs_IQR(citoNum_d2_noFC_wConts,
                                                                                                          groupWise = T), ])
citoNum_d6_nomis <- citoNomis(citoNum_d2_noFC_wConts[!rownames(citoNum_d6_noFC_wConts) %in% detectOLs_IQR(citoNum_d6_noFC_wConts,
                                                                                                          groupWise = T), ])

# Do biplots of the normalized data
ggarrange(pcBiplot(prcomp(apply(citoNum_d2_nomis,
                                2,
                                function(x) (x - mean(x))/sd(x)),
                          scale. = F,
                          center = F)),
          pcBiplot(prcomp(apply(citoNum_d6_nomis,
                                2,
                                function(x) (x - mean(x))/sd(x)),
                          scale. = F,
                          center = F)),
          legend = "bottom",
          common.legend = T)

ggsave(paste0(citoPlotDir, "nomisPCA.pdf"),
       height = 4, width = 7.5)

ggarrange(pcBiplot_date(prcomp(citoNum_d2_nomis, scale. = T, center = T),
                        sampInfo = citoParsed, linNames = T),
          pcBiplot_date(prcomp(citoNum_d6_nomis, scale. = T, center = T),
                        sampInfo = citoParsed, linNames = T),
          legend = "right",
          common.legend = T)

ggsave(paste0(citoPlotDir, "nomisPCA_wDate.pdf"),
       height = 4, width = 9)

# CCMN
library(NormalizeMets)
ccmnNorm_package <- NormQcmets(citoNum_d2_noFC_wConts, method = "ccmn", 
                               qcmets = grep("ctrl", colnames(citoNum_d2_noFC_wConts)),
                               factors=str_extract(rownames(citoNum_d2_noFC_wConts), "[^_]*_[^_]*"))$featuredata
pcBiplot_date(prcomp(ccmnNorm_package,
                     scale. = T,
                     center = T),
              sampInfo = citoParsed, linNames = T)



citoNum_ctrls <- citoNum_d2_noFC_wConts
# Obtain the mean and sd matrixes in order to desescaling the data in the last step of
# the normalization
meanMat <- matrix(rep(1, nrow(citoNum_ctrls))) %*% t(matrix(apply(citoNum_ctrls, 2, mean)))
dimnames(meanMat) <- dimnames(citoNum_ctrls)
meanMat <- meanMat[, !grepl("ctrl", colnames(meanMat))]
sdMat <- diag(apply(citoNum_ctrls, 2, sd))
dimnames(sdMat) <- list(colnames(citoNum_ctrls), colnames(citoNum_ctrls))
sdMat <- sdMat[!grepl("ctrl", colnames(sdMat)),
               !grepl("ctrl", colnames(sdMat))]
citoNum_ctrls <- apply(citoNum_ctrls, 2, function(x) (x - mean(x))/sd(x))
isMat <- citoNum_ctrls[, grepl("ctrl", colnames(citoNum_ctrls))]
anMat <- citoNum_ctrls[, !grepl("ctrl", colnames(citoNum_ctrls))]

factMat <- sapply(unique(str_extract(rownames(citoNum_ctrls), "[^_]*_[^_]*")),
                  function(x) as.numeric(grepl(x, rownames(citoNum_ctrls))))
factMat <- cbind(rep(1,
                     nrow(factMat)),
                 factMat)
rownames(factMat) <- rownames(citoNum_ctrls)
factMat <- factMat[, 1:(ncol(factMat) - 1)]
# Do MLR between isMat and factMat to determine the effect the grouping 
# has on each one of the variables (cell death categories)
beta_ccmn1 <- solve(t(factMat) %*% factMat) %*% (t(factMat) %*% isMat)
# Remove the effect the grouping has on the controls from the controls
isMat_prima <- isMat - factMat %*% beta_ccmn1
# Use PCA to estimate the structured variance of the controls
isPCA <- prcomp(isMat_prima, scale. = F, center = F)
# Plot the PCA to see if the controls are still aggrupated according to date
# (sampInfo is citoParsed DF because the rownames of the pca correspond to the 
# samples each control is associated to, but date is the same so citoParsed works
# for this)
pcBiplot_date(isPCA, sampInfo = citoParsed)
pcLoads <- isPCA$rotation
pcScors <- isPCA$x
# Regress the scores of the PCA on the citometry values of the samples, to
# determine how much the variance of the controls (without grouping influence)
# contributes to variance of the samples.
beta_ccmn2 <- solve(t(pcScors) %*% pcScors) %*% (t(pcScors) %*% anMat)

# Substract the correlated part of the data obtained in past regression 
# from the sample data to obtain normalized values
ccmnNormed <- anMat - pcScors %*% beta_ccmn2
# Unscale the normalized result by adding the means and multiplying by
# the SD
ccmnNormed <- (ccmnNormed %*% sdMat + meanMat) 

citoCcmn <- function(citoNum_ctrls){
        # Obtain the mean and sd matrixes in order to desescaling the data in the last step of
        # the normalization
        meanMat <- matrix(rep(1, nrow(citoNum_ctrls))) %*% t(matrix(apply(citoNum_ctrls, 2, mean)))
        dimnames(meanMat) <- dimnames(citoNum_ctrls)
        meanMat <- meanMat[, !grepl("ctrl", colnames(meanMat))]
        sdMat <- diag(apply(citoNum_ctrls, 2, sd))
        dimnames(sdMat) <- list(colnames(citoNum_ctrls), colnames(citoNum_ctrls))
        sdMat <- sdMat[!grepl("ctrl", colnames(sdMat)),
                       !grepl("ctrl", colnames(sdMat))]
        citoNum_ctrls <- apply(citoNum_ctrls, 2, function(x) (x - mean(x))/sd(x))
        isMat <- citoNum_ctrls[, grepl("ctrl", colnames(citoNum_ctrls))]
        anMat <- citoNum_ctrls[, !grepl("ctrl", colnames(citoNum_ctrls))]
        
        factMat <- sapply(unique(str_extract(rownames(citoNum_ctrls), "[^_]*_[^_]*")),
                          function(x) as.numeric(grepl(x, rownames(citoNum_ctrls))))
        factMat <- cbind(rep(1,
                             nrow(factMat)),
                         factMat)
        rownames(factMat) <- rownames(citoNum_ctrls)
        factMat <- factMat[, 1:(ncol(factMat) - 1)]
        # Do MLR between isMat and factMat to determine the effect the grouping 
        # has on each one of the variables (cell death categories)
        beta_ccmn1 <- solve(t(factMat) %*% factMat) %*% (t(factMat) %*% isMat)
        # Remove the effect the grouping has on the controls from the controls
        isMat_prima <- isMat - factMat %*% beta_ccmn1
        # Use PCA to estimate the structured variance of the controls
        isPCA <- prcomp(isMat_prima, scale. = F, center = F)
        pcLoads <- isPCA$rotation
        pcScors <- isPCA$x
        # Regress the scores of the PCA on the citometry values of the samples, to
        # determine how much the variance of the controls (without grouping influence)
        # contributes to variance of the samples.
        beta_ccmn2 <- solve(t(pcScors) %*% pcScors) %*% (t(pcScors) %*% anMat)
        
        # Substract the correlated part of the data obtained in past regression 
        # from the sample data to obtain normalized values
        ccmnNormed <- anMat - pcScors %*% beta_ccmn2
        # Unscale the normalized result by multiplying by the SD and adding the 
        # means
        ccmnNormed <- (ccmnNormed %*% sdMat + meanMat)
        return(ccmnNormed)
}

detectOLs_IQR(citoNum_d2_noFC_wConts, groupWise = T)

citoNum_d2_ccmn <- citoCcmn(citoNum_d2_noFC_wConts[!rownames(citoNum_d2_noFC_wConts) %in% detectOLs_IQR(citoNum_d2_noFC_wConts,
                                                                                                        groupWise = T), ])
citoNum_d6_ccmn <- citoCcmn(citoNum_d6_noFC_wConts[!rownames(citoNum_d6_noFC_wConts) %in% detectOLs_IQR(citoNum_d6_noFC_wConts,
                                                                                                        groupWise = T), ])

ggarrange(pcBiplot_date(prcomp(citoNum_d2_ccmn, scale. = T, center = T),
                        sampInfo = citoParsed, linNames = T),
          pcBiplot_date(prcomp(citoNum_d6_ccmn, scale. = T, center = T),
                        sampInfo = citoParsed, linNames = T),
          legend = "right",
          common.legend = T)

ggsave(paste0(citoPlotDir, "ccmnPCA_wDate.pdf"),
       height = 4, width = 9)

ggarrange(doRLAPlot(citoNum_d2_noFC, "WG", normMeth = "Unnormalized") +
                  rremove("xlab") +
                  rremove("x.text"), 
          doRLAPlot(citoNum_d2_noFC, "AG", normMeth = "Unnormalized") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"), 
          doRLAPlot(citoNum_d6_noFC, "WG", normMeth = "Unnormalized") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"), 
          doRLAPlot(citoNum_d6_noFC, "AG", normMeth = "Unnormalized") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"),
          
          doRLAPlot(citoNum_d2, "WG", normMeth = "Fold Change") +
                  rremove("xlab") +
                  rremove("x.text"),
          doRLAPlot(citoNum_d2, "AG", normMeth = "Fold Change") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"),
          doRLAPlot(citoNum_d6, "WG", normMeth = "Fold Change") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"),
          doRLAPlot(citoNum_d6, "AG", normMeth = "Fold Change") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"),
          
          doRLAPlot(citoNum_d2_nomis, "WG", normMeth = "NOMIS") +
                  rremove("xlab") +
                  rremove("x.text"), 
          doRLAPlot(citoNum_d2_nomis, "AG", normMeth = "NOMIS") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"),
          doRLAPlot(citoNum_d6_nomis, "WG", normMeth = "NOMIS") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"), 
          doRLAPlot(citoNum_d6_nomis, "AG", normMeth = "NOMIS") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"),
          
          doRLAPlot(citoNum_d2_ccmn, "WG", normMeth = "CCMN") +
                  rremove("xlab") +
                  rremove("x.text"), 
          doRLAPlot(citoNum_d2_ccmn, "AG", normMeth = "CCMN") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"),
          doRLAPlot(citoNum_d6_ccmn, "WG", normMeth = "CCMN") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"), 
          doRLAPlot(citoNum_d6_ccmn, "AG", normMeth = "NOMIS") +
                  rremove("xlab") +
                  rremove("x.text") +
                  rremove("ylab"),
          ncol = 4, nrow = 4)



ggsave(paste0(citoPlotDir, "RLAPlots.pdf"),
       height = 10, width = 20)




# Do univariate tests
##################################################################################################################

# Plot the viability fold change

print(doMannWhit(citoNum_d2_noOLs))

print(doMannWhit(citoNum_d6_noOLs))


# Viability
viableBxplt_FC_d2 <- doCitoBoxplots(citoNum_d2_noOLs, xAxis = "strain", deathCat = "viable", plotStats = T)
viableBxplt_FC_d2

ggsave(paste0(citoPlotDir, "viableBxplt_FC_d2.pdf"),
       height = 4.5, width = 5)

viableBxplt_FC_d6 <- doCitoBoxplots(citoNum_d6_noOLs, xAxis = "strain", deathCat = "viable", plotStats = T)
viableBxplt_FC_d6

ggsave(paste0(citoPlotDir, "viableBxplt_FC_d6.pdf"),
       height = 4.5, width = 5)


# Necrosis
necrotBxplt_FC_d2 <- doCitoBoxplots(citoNum_d2_noOLs, xAxis = "strain", deathCat = "necrotic", plotStats = T)
necrotBxplt_FC_d2

ggsave(paste0(citoPlotDir, "necrotBxplt_FC_d2.pdf"),
       height = 4.5, width = 5)

necrotBxplt_FC_d6 <- doCitoBoxplots(citoNum_d6_noOLs, xAxis = "strain", deathCat = "necrotic", plotStats = T)
necrotBxplt_FC_d6

ggsave(paste0(citoPlotDir, "necrotBxplt_FC_d6.pdf"),
       height = 4.5, width = 5)

# Early apoptosis
earlApBxplt_FC_d2 <- doCitoBoxplots(citoNum_d2_noOLs, xAxis = "strain", deathCat = "early_apoptotic", plotStats = T)
earlApBxplt_FC_d2

ggsave(paste0(citoPlotDir, "earlApBxplt_FC_d2.pdf"),
       height = 4.5, width = 5)

earlApBxplt_FC_d6 <- doCitoBoxplots(citoNum_d6_noOLs, xAxis = "strain", deathCat = "early_apoptotic", plotStats = T)
earlApBxplt_FC_d6

ggsave(paste0(citoPlotDir, "earlApBxplt_FC_d6.pdf"),
       height = 4.5, width = 5)

# Late apoptosis
lateApBxplt_FC_d2 <- doCitoBoxplots(citoNum_d2_noOLs, xAxis = "strain", deathCat = "late_apoptotic", plotStats = T)
lateApBxplt_FC_d2

ggsave(paste0(citoPlotDir, "lateApBxplt_FC_d2.pdf"),
       height = 4.5, width = 5)

lateApBxplt_FC_d6 <- doCitoBoxplots(citoNum_d6_noOLs, xAxis = "strain", deathCat = "late_apoptotic", plotStats = T)
lateApBxplt_FC_d6

ggsave(paste0(citoPlotDir, "lateApBxplt_FC_d6.pdf"),
       height = 4.5, width = 5)


# Do combined boxplots by days

ggarrange(viableBxplt_FC_d2$plot,
          necrotBxplt_FC_d2$plot,
          legend = "bottom",
          common.legend = T)

ggsave(paste0(citoPlotDir, "necrotViablBxplts_FC_d2.pdf"),
       height = 5, width = 9)


ggarrange(viableBxplt_FC_d6$plot,
          necrotBxplt_FC_d6$plot,
          legend = "bottom",
          common.legend = T)

ggsave(paste0(citoPlotDir, "necrotViablBxplts_FC_d6.pdf"),
       height = 5, width = 9)


ggarrange(earlApBxplt_FC_d2$plot,
          lateApBxplt_FC_d2$plot,
          legend = "bottom",
          common.legend = T)

ggsave(paste0(citoPlotDir, "apoptBxplts_FC_d2.pdf"),
       height = 5, width = 9)

ggarrange(earlApBxplt_FC_d6$plot,
          lateApBxplt_FC_d6$plot,
          legend = "bottom",
          common.legend = T)

ggsave(paste0(citoPlotDir, "apoptBxplts_FC_d6.pdf"),
       height = 5, width = 9)

ggarrange(earlApBxplt_FC_d2$plot,
          lateApBxplt_FC_d2$plot,
          earlApBxplt_FC_d6$plot,
          lateApBxplt_FC_d6$plot,
          legend = "bottom",
          common.legend = T)

ggsave(paste0(citoPlotDir, "apoptBxplts_FC_d26.pdf"),
       height = 9.2, width = 9)

# Obtain tile plot of p-values

doMannWhit(citoNum_d2_noOLs, whatComp = "strains")

getBtwStrainsStatTilePlot <- function(citoNum,
                                      deathCat,
                                      orientation = "horizontal",
                                      color = T){
        statDeath <- doMannWhit(citoNum = citoNum,
                                whatComp = "strains")[[deathCat]]
        statSymbDF <- data.frame(value = c(1, 0.1, 0.05, 0.01, 0.001),
                                 symbol = c("ns", ".", "*", "**", "***"),
                                 color = c("#FF0000", "#FFA500", "#7CFC00", "#228B22", "#32CD32"),
                                 stringsAsFactors = F)
        plotLst <- list()
        for(h in names(statDeath)){
                mwDF <- statDeath[[h]]
                dimnames(mwDF) <- lapply(dimnames(mwDF),
                                         function(x) gsub("bovis", "*M. bovis*", x))
                dimnames(mwDF) <- lapply(dimnames(mwDF),
                                         function(x) gsub("chimp", "Chimp. B.", x))
                mwDF[is.nan(as.matrix(mwDF))] <- 1
                mwDF$strainA <- rownames(mwDF)
                mwDF <- mwDF[, c("strainA",
                                 colnames(mwDF)[1:(ncol(mwDF) - 1)])]
                tileDF <- as.data.frame(pivot_longer(mwDF,
                                                     cols = colnames(mwDF)[2:ncol(mwDF)],
                                                     names_to = "strainB"))
                tileDF <- tileDF[!is.na(tileDF$value), ]
                tileDF$color <- sapply(tileDF$value,
                                       function(x) statSymbDF$color[max(which(x <= statSymbDF$value))])
                tileDF$symbol <- sapply(tileDF$value, function(x) statSymbDF$symbol[max(which(x <= statSymbDF$value))])
                plotCols <- statSymbDF$color[match(levels(factor(tileDF$symbol)),
                                                   statSymbDF$symbol)]
                plotTitle <- h
                plotTitle <- gsub("thp1", "THP1", plotTitle)
                plotTitle <- gsub("bomac", "BoMac", plotTitle)
                #plotTitle <- sprintf("**%s**", plotTitle)
                tileDF$label <- formatC(tileDF$value,
                                        format = "e",
                                        digits = 2)
                tileDF$face <- rep("plain", nrow(tileDF))
                tileDF$face[tileDF$value <= 0.05] <- "bold"
                if(color){
                        p <- ggplot(tileDF,
                                    aes(x = strainB,
                                        y = factor(strainA,
                                                   levels = sort(unique(strainA),
                                                                 decreasing = T)),
                                        fill=symbol,
                                        label=label)) +
                                scale_fill_manual(values = plotCols) +
                                geom_tile(color="black") +
                                geom_text(aes(label=label),
                                          size=5,
                                          color="white")
                }else{
                        p <- ggplot(tileDF,
                                    aes(x = strainB,
                                        y = factor(strainA,
                                                   levels = sort(unique(strainA),
                                                                 decreasing = T)),
                                        fill = "none",
                                        label=label)) +
                                scale_fill_manual(values = "white") +
                                geom_tile(color="black") +
                                geom_text(aes(label=label,
                                              fontface = face),
                                          size=5,
                                          color="black")
                }
                p <- p +

                        scale_x_discrete(position = "top") +
                        labs(title =plotTitle) +
                        theme(title = element_text(size = 16, face = "bold"),
                              axis.text.x = ggtext::element_markdown(size = 15),
                              axis.text.y = ggtext::element_markdown(size = 15),
                              axis.title.x = element_blank(),
                              axis.title.y = element_blank(),
                              panel.background = element_blank(),
                              #panel.border = element_rect(colour = "black", fill=NA, size=1),
                              #panel.grid.major = element_line(colour = "#d4d4d4"),
                              legend.position = "none")
                plotLst[[h]] <- p
        }
        if(orientation == "horizontal"){
                plotGlob <- ggarrange(plotLst[[2]],
                                      plotLst[[1]])
        }else if(orientation == "vertical"){
                plotGlob <- ggarrange(plotLst[[2]],
                                      plotLst[[1]],
                                      ncol = 1,
                                      nrow = 2)
        }
        
        plotGlob
}

getBtwStrainsStatTilePlot(citoNum_d6_noOLs,
                          deathCat = "necrotic")

ggsave(paste0(citoPlotDir, "mannWStrainsNecrotD6_tilePlot.pdf"),
       width = 9, height = 2)

getBtwStrainsStatTilePlot(citoNum_d6_noOLs,
                          deathCat = "necrotic",
                          color = F)

ggsave(paste0(citoPlotDir, "mannWStrainsNecrotD6_tilePlot_noCol.pdf"),
       width = 9, height = 2)



getBtwStrainsStatTilePlot(citoNum_d6_noOLs,
                          deathCat = "late_apoptotic",
                          orientation = "vertical")

ggsave(paste0(citoPlotDir, "mannWStrainslateApD6_tilePlot.pdf"),
       width = 4.5, height = 4)

getBtwStrainsStatTilePlot(citoNum_d6_noOLs,
                          deathCat = "late_apoptotic",
                          orientation = "vertical",
                          color = F)

ggsave(paste0(citoPlotDir, "mannWStrainslateApD6_tilePlot_noCol.pdf"),
       width = 4.5, height = 4)



getBtwStrainsStatTilePlot(citoNum_d6_noOLs,
                          deathCat = "late_apoptotic")

ggsave(paste0(citoPlotDir, "mannWStrainslateApD6_tilePlot_horiz.pdf"),
       width = 9, height = 2)

getBtwStrainsStatTilePlot(citoNum_d6_noOLs,
                          deathCat = "late_apoptotic",
                          color = F)

ggsave(paste0(citoPlotDir, "mannWStrainslateApD6_tilePlot_horiz_noCol.pdf"),
       width = 9, height = 2)



getBtwStrainsStatTilePlot(citoNum_d6_noOLs,
                          deathCat = "early_apoptotic")

ggsave(paste0(citoPlotDir, "mannWStrainsEarlyApD6_tilePlot_horiz.pdf"),
       width = 9, height = 2)

getBtwStrainsStatTilePlot(citoNum_d6_noOLs,
                          deathCat = "early_apoptotic",
                          color = F)

ggsave(paste0(citoPlotDir, "mannWStrainsEarlyApD6_tilePlot_horiz_noCol.pdf"),
       width = 9, height = 2)