if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(ggpattern)) install.packages("ggpattern")
library(ggpattern)
if(!require(ggpubr)) install.packages("ggpubr")
library(ggpubr)

# Directory stuff
##################################################################################################################

rootDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcPCRs/"


dataDir <- paste0(rootDir, "data/")
plotDir <- paste0(rootDir, "plots/")
moiPlotDir <- paste0(plotDir, "moi/")

if(!dir.exists(moiPlotDir)){
        dir.create(moiPlotDir)
}

# Load and parse the data
##################################################################################################################
moiData <- read.csv(paste0(dataDir,
                           "moiTests.csv"))


doMoiBxplt <- function(data, strain, time, title = T, color = T){
        strainCols <- data.frame(color = c("#f279ce", "#f279ce", "#c4bf62", "#24ad37", "#871414"),
                                 strain = c("BCG", "M. bovis", "Chimp. B.", "L6", "L5"),
                                 stringsAsFactors = F)
        colorSt <- strainCols$color[grep(tolower(strain),
                                         tolower(strainCols$strain))]
        strainN <- strainCols$strain[grep(tolower(strain),
                                          tolower(strainCols$strain))]
        data <- data[tolower(data$strain) == tolower(strain), ]
        data <- data[data$time == time, ]
        data$strain <- gsub("Bovis",
                            "*M. bovis*",
                            data$strain)
        data$strain <- gsub("chimp|Chimp",
                            "Chimp. B.",
                            data$strain)
        data$host <- factor(data$host,
                            levels = c("THP1", "BoMac"))
        data$moi <- factor(data$moi)
        plotTitle <- sprintf("%s, %s HPI", strainN, as.character(time))
        if(color){
                boxPlot <- ggplot(data, aes(x = moi, fill = strain, y = inf_perc)) +
                        scale_fill_manual(values = colorSt, 
                                          labels = strainN) 
        }else{
                boxPlot <- ggplot(data, aes(x = moi, y = inf_perc))
        }
        boxPlot <- boxPlot +
                geom_boxplot_pattern(pattern_color = "black",
                                     pattern_fill = "black",
                                     pattern_spacing = 0.015,
                                     aes(pattern = host)) + 
                scale_pattern_manual(values = c("stripe", "none"),
                                     labels = c("THP1", "BoMac")) +
                theme(title = ggtext::element_markdown(),
                      axis.text.x = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      legend.text = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"))
        if(title){
                boxPlot <- boxPlot + labs(title = plotTitle,
                                          x = "MOI",
                                          y = "percentage of infected cells")
        }else{
                boxPlot <- boxPlot + labs(x = "MOI",
                                          y = "percentage of infected cells")
        }
        boxPlot
}

BCG_1hpi <- doMoiBxplt(moiData, strain = "BCG", time = 1)

bovis_1hpi <- doMoiBxplt(moiData, strain = "bovis", time = 1)

chimp_1hpi <- doMoiBxplt(moiData, strain = "chimp", time = 1)

L6_1hpi <- doMoiBxplt(moiData, strain = "L6", time = 1)

L5_1hpi <- doMoiBxplt(moiData, strain = "L5", time = 1)


BCG_4hpi <- doMoiBxplt(moiData, strain = "BCG", time = 4)

bovis_4hpi <- doMoiBxplt(moiData, strain = "bovis", time = 4)

chimp_4hpi <- doMoiBxplt(moiData, strain = "chimp", time = 4)

L6_4hpi <- doMoiBxplt(moiData, strain = "L6", time = 4)

L5_4hpi <- doMoiBxplt(moiData, strain = "L5", time = 4)


ggarrange(BCG_1hpi + 
                  rremove("xlab"),
          bovis_1hpi + 
                  rremove("ylab") + 
                  rremove("xlab"),
          chimp_1hpi + 
                  rremove("ylab") + 
                  rremove("xlab"),
          L6_1hpi + 
                  rremove("ylab") + 
                  rremove("xlab"),
          L5_1hpi + 
                  rremove("ylab") + 
                  rremove("xlab"),
          BCG_4hpi,
          bovis_4hpi + 
                  rremove("ylab"),
          chimp_4hpi + 
                  rremove("ylab"),
          L6_4hpi + 
                  rremove("ylab"),
          L5_4hpi + 
                  rremove("ylab"),
          ncol = 5,
          nrow = 2,
          common.legend = T)

ggsave(paste0(moiPlotDir, "allStrainsMoiTsts.pdf"), height = 6, width = 12)

BCG_1hpi_noTitle_noCol <- doMoiBxplt(moiData, strain = "BCG", time = 1, title = F, color = F)
BCG_4hpi_noTitle_noCol <- doMoiBxplt(moiData, strain = "BCG", time = 4, title = F, color = F)

ggarrange(BCG_1hpi_noTitle_noCol,
          BCG_4hpi_noTitle_noCol,
          common.legend = T,
          legend = "bottom")

ggsave(paste0(moiPlotDir, "BCGMoiTsts.pdf"), height = 3.5, width = 6)

BCG_4hpi_noTitle_noCol

ggsave(paste0(moiPlotDir, "BCGMoiTsts_4H.pdf"), height = 3.5, width = 4)
