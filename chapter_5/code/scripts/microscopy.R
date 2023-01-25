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
if(!require(readxl)) install.packages("readxl")
library(readxl)
if(!require(ggpattern)) install.packages("ggpattern")
library(ggpattern)


# Directory stuff
##################################################################################################################
rootDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcPCRs/"


dataDir <- paste0(rootDir, "data/")
plotDir <- paste0(rootDir, "plots/")
microPlotDir <- paste0(plotDir, "micro/")
if(!dir.exists(microPlotDir)){
        dir.create(microPlotDir)
}

# Load and parse the data
##################################################################################################################
microData <- data.frame(readxl::read_xlsx(paste0(dataDir,
                                                 "tabla_infecciones_completa_updated3.xlsx"),
                                          sheet = 1))

microData <- microData[, colnames(microData) %in% c("Cepas", "Célula", "Día_micro", "Valor_micro", "Fecha_INF")]

colnames(microData) <- c("date", "strain", "host", "dpi", "infection_ratio")

microData <- microData[!is.na(microData$infection_ratio), ]

dateNaIdxs <- c(which(!is.na(microData$date)),
                (length(microData$date) + 1))
dateVec <- c()
for(i in 1:(length(dateNaIdxs) - 1)){
        n <- dateNaIdxs[i]
        date <- as.character(microData$date[n])
        idx1 <- n + 1
        idx2 <- dateNaIdxs[i + 1]
        print(idx1)
        print(idx2)
        dateVec <- c(dateVec, date, rep(date, idx2 - idx1))
}

microData$date <- dateVec

# At day 3 in L5 infections there is just one data point because all macrophages were dead: there was 
# nothing at sight. We can consider this as 100% macrophages infected, so let's add this to the 
# dataframe.
#L5D3 <- data.frame(matrix(c("2021-07-06", "L5", "THP1",
#                            "2021-07-06", "L5", "THP1",
#                            "2021-07-06", "L5", "THP1",
#                            "2021-07-06", "L5", "THP1"),
#                          byrow = T,
#                          nrow = 4,
#                          ncol = 3,
#                          dimnames = list(NULL,
#                                          colnames(microData)[1:3])))

#L5D3$dpi <- rep(3, nrow(L5D3))
#L5D3$infection_ratio <- rep(100, nrow(L5D3))

#microData <- rbind.data.frame(microData,
#                              L5D3)

# Add a column with unique sample names
casesMicro <- paste(microData$strain, microData$host, microData$dpi, sep = "_")
sampVec <- rep(NA, length(casesMicro))
for(c in unique(casesMicro)){
        sampVec[casesMicro == c] <- paste0(paste(c, "R", sep = "_"),
                                           1:sum(casesMicro == c))
}

microData$sample <- sampVec
microData <- microData[, c("date", "strain", "host", "sample", "dpi", "infection_ratio")]

# Functions
##################################################################################################################
# Detect outliers based on Interquartile range (IQR) method. k refers to the factor
# the IQR is multiplied to. groupWise = T refers to the IQR is calculated based 
# on the replicates. If it is false IQR is calculated based in all the samples
detectOLs_IQR <- function(data, k = 1.5, groupWise = F, logDF = F){
        if(!require(stringr)) install.packages("stringr")
        library(stringr)
        data$case <- paste(data$strain, data$host, as.character(data$dpi), sep =  "_")
        dataNum <- data.frame(infection_ratio = data$infection_ratio)
        rownames(dataNum) <- data$sample
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
                                             ncol = ncol(dataNum),
                                             dimnames = list(NULL,
                                                             colnames(dataNum))))
                unCases <- unique(data$case)
                for(c in unCases){
                        numCase <- dataNum[grepl(c, rownames(dataNum)), , drop = F]
                        subLog <- getLogic(numCase, k)
                        logicDF <- rbind.data.frame(logicDF,
                                                    subLog)
                }
        }else{
                logicDF <- getLogic(dataNum, k)
        }
        if(logDF){
                print(logicDF)
        }
        OLs <- apply(logicDF, 1, any)
        OLs <- names(OLs)[OLs]
        return(OLs)
}

# Do boxplot of the infection ratio data with statistics (Mann-Whitney U test). Input should be microscopy data
# of a particular day post infection.

doMicroBxplt <- function(dataPlot, col = "strain"){
        statSymbDF <- data.frame(value = c(1, 0.1, 0.05, 0.01, 0.001),
                                 symbol = c("ns", ".", "*", "**", "***"))
        dataPlot$strain <- gsub("Bovis", "*M. bovis*", dataPlot$strain)
        dataPlot$strain <- gsub("Chimp", "Chimp. B.", dataPlot$strain)
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
        yLab <- "% of infected macrophages"
        plot <- ggplot(dataPlot, aes(y = infection_ratio, x = strain, fill = strain)) +
                geom_boxplot_pattern(pattern_color = "black",
                                     pattern_fill = "black",
                                     pattern_spacing = 0.015,
                                     aes(pattern = host),
                                     outlier.shape = NA) + 
                scale_pattern_manual(values = c("stripe", "none")) +
                labs(x = "Strain",
                     y = yLab)
        if(col == "strain"){
                plot <- plot +
                        geom_point(aes(col = host), position=position_jitterdodge(jitter.width = 0.005)) +
                        scale_color_manual(values = c("#000000", "#010101")) +
                        scale_fill_manual(values = c("#f279ce", "#c4bf62", "#24ad37", "#871414"), 
                                          labels = c("*M. bovis*", "Chimp. B.", "L6", "L5")) +
                        theme(title = ggtext::element_markdown(),
                              axis.title.x = ggtext::element_markdown(),
                              axis.text.x = ggtext::element_markdown(),
                              legend.text = ggtext::element_markdown(size = 8),
                              legend.title = ggtext::element_markdown(size = 10),
                              legend.key.height = unit(14, "points"),
                              legend.key.width = unit(14, "points"),
                              panel.background = element_blank(),
                              panel.border = element_rect(colour = "black", fill=NA, size=1),
                              panel.grid.major = element_line(colour = "#d4d4d4"),
                              legend.position = "bottom")
        }else if(col == "date"){
                plot <- plot +
                        geom_point(aes(shape = host, col = date), position=position_jitterdodge(jitter.width = 0.0001)) +
                        scale_fill_manual(values = alpha(c("#f279ce", "#c4bf62", "#24ad37", "#871414"), .3), 
                                          labels = c("*M. bovis*", "Chimp. B.", "L6", "L5")) +
                        theme(title = ggtext::element_markdown(),
                              axis.title.x = ggtext::element_markdown(),
                              axis.text.x = ggtext::element_markdown(),
                              legend.text = ggtext::element_markdown(size = 8),
                              legend.title = ggtext::element_markdown(size = 10),
                              legend.key.height = unit(14, "points"),
                              legend.key.width = unit(14, "points"),
                              panel.background = element_blank(),
                              panel.border = element_rect(colour = "black", fill=NA, size=1),
                              panel.grid.major = element_line(colour = "#d4d4d4"),
                              legend.position = "right")
        }
        
        maxYs <- sapply(levels(dataPlot$strain),
                        function(x) max(dataPlot$infection_ratio[dataPlot$strain == x]))
        maxYs <- maxYs + max(maxYs) * 0.025
        pVals <- sapply(levels(dataPlot$strain),
                        function(x) wilcox.test(dataPlot$infection_ratio[dataPlot$strain == x & dataPlot$host == "THP1"],
                                                dataPlot$infection_ratio[dataPlot$strain == x & dataPlot$host == "BoMac"])$p.value)
        symbols <- sapply(pVals, function(x) statSymbDF$symbol[max(which(x <= statSymbDF$value))])
        
        statDF <- data.frame(matrix(c(pVals,
                                      sapply(levels(dataPlot$strain),
                                             function(x) median(dataPlot$infection_ratio[dataPlot$strain == x & dataPlot$host == "THP1"])),
                                      sapply(levels(dataPlot$strain),
                                             function(x) median(dataPlot$infection_ratio[dataPlot$strain == x & dataPlot$host == "BoMac"])),
                                      sapply(levels(dataPlot$strain),
                                             function(x) mean(dataPlot$infection_ratio[dataPlot$strain == x & dataPlot$host == "THP1"])),
                                      sapply(levels(dataPlot$strain),
                                             function(x) mean(dataPlot$infection_ratio[dataPlot$strain == x & dataPlot$host == "BoMac"]))),
                                    byrow = T,
                                    nrow = 5, 
                                    ncol = length(pVals),
                                    dimnames = list(c("p-values_MW", "THP-1_median", "BoMac_median", "THP-1_mean", "BoMac_mean"),
                                                    names(pVals))))
        
        signifDF <- data.frame(x = 1:(length(levels(dataPlot$strain))) - 0.2,
                               xend = 1:(length(levels(dataPlot$strain))) + 0.2,
                               y = maxYs,
                               annotation = symbols,
                               group = 1:length(maxYs))
        
        plot <- plot + 
                geom_signif(stat = "identity",
                            data = signifDF,
                            aes(x=x,
                                xend=xend,
                                y=y,
                                yend=y,
                                fill = NULL,
                                annotation = annotation,
                                group = group,
                                col = NULL))
        
        out <- list(stats = statDF, bxplt = plot)
        return(out)
}

# Do univariate (MW tests) between strains infecting the same macrophage and plot the 
# resulting p-values as tile plots.
doInfTilePlots <- function(dataPlot, revX = T, orientation = "vertical", color = T){
        statSymbDF <- data.frame(value = c(1, 0.1, 0.05, 0.01, 0.001),
                                 symbol = c("ns", ".", "*", "**", "***"),
                                 color = c("#FF0000", "#FFA500", "#7CFC00", "#228B22", "#32CD32"),
                                 stringsAsFactors = F)
        
        dataPlot$host <- factor(dataPlot$host,
                                levels = sort(unique(dataPlot$host),
                                              decreasing = T))
        dataPlot$strain <- gsub("Bovis", "*M. bovis*", dataPlot$strain)
        dataPlot$strain <- gsub("Chimp", "Chimp. B.", dataPlot$strain)
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
        hosts <- levels(dataPlot$host)
        lins <- levels(dataPlot$strain)
        plotLst <- list()
        for(h in hosts){
                infHost <- dataPlot[dataPlot$host == h, ]
                infHostPVals <- data.frame(matrix(nrow = length(lins) - 1,
                                                  ncol = length(lins) - 1,
                                                  dimnames = list(lins[1:(length(lins) - 1)],
                                                                  lins[2:length(lins)])))
                for(j in 1:(length(lins) - 1)){
                        strainA <- lins[j]
                        for(k in (j + 1):length(lins)){
                                strainB <- lins[k]
                                pVal <- wilcox.test(infHost$infection_ratio[infHost$strain == strainA],
                                                    infHost$infection_ratio[infHost$strain == strainB])$p.value
                                infHostPVals[strainA, strainB] <- pVal
                        }
                }
                infHostPVals$strainA <- rownames(infHostPVals)
                infHostPVals <- infHostPVals[, c("strainA",
                                                 colnames(infHostPVals)[1:(ncol(infHostPVals) - 1)])]
                tileDF <- as.data.frame(pivot_longer(infHostPVals,
                                                     cols = colnames(infHostPVals)[2:ncol(infHostPVals)],
                                                     names_to = "strainB"))
                tileDF <- tileDF[!is.na(tileDF$value), ]
                tileDF$color <- sapply(tileDF$value,
                                       function(x) statSymbDF$color[max(which(x <= statSymbDF$value))])
                tileDF$symbol <- sapply(tileDF$value, function(x) statSymbDF$symbol[max(which(x <= statSymbDF$value))])
                tileDF$label <- formatC(tileDF$value,
                                        format = "e",
                                        digits = 2)
                tileDF$face <- rep("plain", nrow(tileDF))
                tileDF$face[tileDF$value <= 0.05] <- "bold"
                targStrains <- c("*M. bovis*",
                                 "Chimp. B.",
                                 "L6",
                                 "L5",
                                 "ctrl")
                if(revX){
                        targStrains <- targStrains[length(targStrains):1]
                }
                tileDF$strainB <- factor(tileDF$strainB,
                                         levels = targStrains[targStrains %in% tileDF$strainB])
                
                plotCols <- statSymbDF$color[match(levels(factor(tileDF$symbol)),
                                                   statSymbDF$symbol)]
                plotTitle <- h
                plotTitle <- gsub("thp1", "THP1", plotTitle)
                plotTitle <- gsub("bomac", "BoMac", plotTitle)
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
                              legend.position = "none")
                plotLst[[h]] <- p
        }
        
        if(orientation == "horizontal"){
                plotGlob <- ggarrange(plotLst[[1]],
                                      plotLst[[2]])
        }else if(orientation == "vertical"){
                plotGlob <- ggarrange(plotLst[[1]],
                                      plotLst[[2]],
                                      ncol = 1,
                                      nrow = 2)
        }
        return(plotGlob)
}

# Plots time series of infection ratio per strain in a given host
doTimeInfPlot <- function(inData, host, CI = 0.95, plotSamps = T, yLims = c(10, 100)){
        plotData <- inData[inData$host == host, ]
        meansDF <- data.frame(matrix(nrow = 0,
                                     ncol = 6,
                                     dimnames = list(NULL,
                                                     c(colnames(plotData)[!colnames(plotData) %in% c("date", "sample")],
                                                       "upper_bound",
                                                       "lower_bound"))))
        dpis <- unique(plotData$dpi)
        strains <- unique(plotData$strain)
        
        for(dpi in dpis){
                for(strain in strains){
                        strDpi <- plotData$infection_ratio[plotData$strain == strain & plotData$dpi == dpi]
                        strDpiMean <- mean(strDpi)
                        strDpiSD <- sd(strDpi)
                        CI_dev <- CI*(strDpiSD/sqrt(length(strDpi)))
                        upLim <- strDpiMean + CI_dev
                        lwLim <- strDpiMean - CI_dev
                        numToBind <- data.frame(matrix(c(dpi, strDpiMean, upLim, lwLim),
                                                       ncol = 4,
                                                       nrow = 1,
                                                       dimnames = list(NULL,
                                                                       colnames(meansDF)[3:6])))
                        strToBind <- data.frame(matrix(c(strain, host),
                                                       ncol = 2,
                                                       nrow = 1,
                                                       dimnames = list(NULL,
                                                                       colnames(meansDF)[1:2])))
                        toBindDF <- cbind.data.frame(strToBind,
                                                     numToBind)
                        meansDF <- rbind.data.frame(meansDF, toBindDF)
                }
        }
        
        meansDF$strain <- gsub("Bovis", "*M. bovis*", meansDF$strain)
        meansDF$strain <- gsub("Chimp", "Chimp. B.", meansDF$strain)
        meansDF$strain <- factor(meansDF$strain,
                                 levels = c("*M. bovis*",
                                            "Chimp. B.",
                                            "L6",
                                            "L5",
                                            "ctrl")[c("*M. bovis*",
                                                      "Chimp. B.",
                                                      "L6",
                                                      "L5",
                                                      "ctrl") %in% meansDF$strain])
        
        plotData$strain <- gsub("Bovis", "*M. bovis*", plotData$strain)
        plotData$strain <- gsub("Chimp", "Chimp. B.", plotData$strain)
        plotData$strain <- factor(plotData$strain,
                                  levels = c("*M. bovis*",
                                             "Chimp. B.",
                                             "L6",
                                             "L5",
                                             "ctrl")[c("*M. bovis*",
                                                       "Chimp. B.",
                                                       "L6",
                                                       "L5",
                                                       "ctrl") %in% plotData$strain])
        
        yLab <- sprintf("of infected %s cells", host)
        yLab <- paste0("% ", yLab)
        xLab <- "Days post-infection"
        plt <- ggplot(meansDF, mapping = aes(x = dpi, y = infection_ratio, col = strain)) +
                geom_line() +
                geom_ribbon(aes(ymin=lower_bound, ymax=upper_bound, x=dpi, fill = strain), alpha = 0.3, colour = NA) +
                scale_color_manual(values = c("#f279ce", "#c4bf62", "#24ad37", "#871414"), 
                                   labels = c("*M. bovis*", "Chimp. B.", "L6", "L5")) +
                scale_fill_manual(values = c("#f279ce", "#c4bf62", "#24ad37", "#871414"), 
                                  labels = c("*M. bovis*", "Chimp. B.", "L6", "L5")) +
                
                labs(x = xLab,
                     y = yLab) +
                ylim(yLims) +
                theme(title = ggtext::element_markdown(),
                      axis.title.x = ggtext::element_markdown(),
                      axis.text.x = ggtext::element_markdown(),
                      legend.text = ggtext::element_markdown(),
                      legend.title = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "bottom")
        
        if(plotSamps){
                plt <- plt +
                        geom_point(plotData, mapping = aes(x = dpi, y = infection_ratio, col = strain))
        }else{
                plt <- plt + 
                        geom_point()
        }
        return(plt)
}

# Remove outliers from microscopy data
##################################################################################################################

OLs <- detectOLs_IQR(microData, groupWise = T, k = 1.5)

microData <- microData[!microData$sample %in% OLs, ]

# Do boxplots of 0 dpi
##################################################################################################################

# Get day 0 infection ratios
d0_micro <- microData[microData$dpi == 0, ]

microD0Plot <- doMicroBxplt(d0_micro)

microD0Plot$stats
microD0Plot$bxplt

ggsave(paste0(microPlotDir, "d0_infRatio.pdf"),
       height = 5, width = 5)

microD0Plot_byDate <- doMicroBxplt(d0_micro, col = "date")

microD0Plot_byDate$bxplt

ggsave(paste0(microPlotDir, "d0_infRatio_byDate.pdf"),
       height = 5, width = 5)

# Get day 1 infection ratios
d1_micro <- microData[microData$dpi == 1, ]

microD1Plot <- doMicroBxplt(d1_micro)

microD1Plot$stats
microD1Plot$bxplt

ggsave(paste0(microPlotDir, "d1_infRatio.pdf"),
       height = 5, width = 5)

microD1Plot_byDate <- doMicroBxplt(d1_micro, col = "date")

microD1Plot_byDate$bxplt

ggsave(paste0(microPlotDir, "d1_infRatio_byDate.pdf"),
       height = 5, width = 5)

# Get day 3 infection ratios
d3_micro <- microData[microData$dpi == 3, ]

microD3Plot <- doMicroBxplt(d3_micro)

microD3Plot$stats
microD3Plot$bxplt

ggsave(paste0(microPlotDir, "d3_infRatio.pdf"),
       height = 5, width = 5)

microD3Plot_byDate <- doMicroBxplt(d3_micro, col = "date")

microD3Plot_byDate$bxplt

ggsave(paste0(microPlotDir, "d3_infRatio_byDate.pdf"),
       height = 5, width = 5)

ggarrange(microD0Plot_byDate$bxplt,
          microD1Plot_byDate$bxplt,
          microD3Plot_byDate$bxplt,
          nrow = 1,
          common.legend = T,
          legend = "bottom")

ggsave(paste0(microPlotDir, "d013_infRatio_byDate.pdf"),
       height = 5, width = 10)

# Do tileplots of 0 dpi
##################################################################################################################
doInfTilePlots(d0_micro, revX = T, orientation = "vertical")
ggsave(paste0(microPlotDir, "d0_tilePlot.pdf"),
       height = 4, width = 4)

doInfTilePlots(d0_micro, revX = T, orientation = "vertical", color = F)
ggsave(paste0(microPlotDir, "d0_tilePlot_noCol.pdf"),
       height = 4, width = 4)

# Do timeseries of infection ratio
##################################################################################################################
ggarrange(doTimeInfPlot(microData, "THP1", CI = 0.95, plotSamps = F, yLims = c(0, 100)),
          doTimeInfPlot(microData, "BoMac", CI = 0.95, plotSamps = F, yLims = c(0, 100)),
          common.legend = T,
          legend = "bottom")
ggsave(paste0(microPlotDir, "infTimeProg.pdf"),
       height = 4, width = 7)

ggarrange(doTimeInfPlot(microData, "THP1", CI = 0.95, plotSamps = T, yLims = c(0, 100)),
          doTimeInfPlot(microData, "BoMac", CI = 0.95, plotSamps = T, yLims = c(0, 100)),
          common.legend = T,
          legend = "bottom")
ggsave(paste0(microPlotDir, "infTimeProg_wSamps.pdf"),
       height = 4, width = 7)

ggsave(paste0(microPlotDir, "infTimeProg_wSamps.png"),
       height = 4, width = 7)
