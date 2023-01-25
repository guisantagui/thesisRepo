if(!require(factoextra)) install.packages("factoextra")
library("factoextra")
if(!require(ropls)) install.packages("ropls")
library("ropls")
if(!require(caret)) install.packages("caret")
library("caret")
if(!require(gplots)) install.packages("gplots")
library("gplots")
if(!require(ggrepel)) install.packages("ggrepel")
library("ggrepel")

# Input variables
####################################################################################################
whatMod <- "delsAllSNPs"
whatMed <- "mi7H9OADCMed"

# Directory paths
####################################################################################################
inpDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/results/simulations/"
resDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/results/"

fbaDir <- sprintf("%s%s/FBA/%s/", inpDir, whatMod, whatMed)
        
medName <- gsub("mi|Med", "", whatMed)

plotDir <- paste0(resDir, "plots/")

if(!dir.exists(plotDir)){
        dir.create(plotDir)
}

# Load data 
####################################################################################################
# All flux distributions
fba <- read.csv(sprintf("%sfbaFlxs_%s.csv", fbaDir, medName), row.names = 1)
# Flux distributions with the reactions directly altered when building models removed
fba_noRemRXNs <- read.csv(sprintf("%sfbaFlxs_%s_noRemRXNs.csv", fbaDir, medName), row.names = 1)
# Flux distributions when h2o2 is present
fbaH2O2 <- read.csv(sprintf("%sfbaFlxs_%s_h2o2.csv", fbaDir, medName), row.names = 1)
# Flux distributions when Fe3 is scarce
fbaLowFe <- read.csv(sprintf("%sfbaFlxs_%s_lowFe.csv", fbaDir, medName), row.names = 1)
# Flux distributions in starvation
fbaStarv <- read.csv(sprintf("%sfbaFlxs_%s_starv.csv", fbaDir, medName), row.names = 1)

# Create vector of lineage colors
lineages <- sort(rownames(fba))

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

names(linCols) <- lineages

# Do PCA
fba_PCA <- prcomp(fba)
fba_noRemRXNs_PCA <- prcomp(fba_noRemRXNs)
fba_h2o2_PCA <- prcomp(fbaH2O2)
fba_lowFe_PCA <- prcomp(fbaLowFe)
fba_starv_PCA <- prcomp(fbaStarv)

pcaPlot <- function(pca, x="PC1", y="PC2"){
        data <- data.frame(obsnames=row.names(pca$x), pca$x)
        data <- data[, c("obsnames", x, y)]
        data <- data[order(data$obsnames), ]
        data$obsnames <- factor(data$obsnames, levels = data$obsnames)
        data$values <- linCols[match(data$obsnames, names(linCols))]
        propVar <- summary(pca)$importance[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)
        #return(data)
        plot <- ggplot(data, aes(data[, 2], 
                                 data[, 3], 
                                 label = obsnames, 
                                 color = obsnames)) + 
                scale_discrete_manual("Lineage",
                                      aesthetics = "colour",
                                      values = data$values) +
                geom_hline(yintercept = 0, alpha = 0.6) +
                geom_vline(xintercept = 0, alpha = 0.6) +
                geom_point() + 
                xlab(sprintf("%s (%s%%)", x, propX)) +
                ylab(sprintf("%s (%s%%)", y, propY)) +
                geom_text_repel() +
                theme_minimal()
        plot
}


pcaPlotPath <- sprintf("%sfba_pca_%s.pdf", plotDir, medName)

pdf(pcaPlotPath, height = 4, width = 4.5)
pcaPlot(fba_PCA)
dev.off()


pcaPlotPath_noRemRXNs <- sprintf("%sfba_pca_%s_noRemRXNs.pdf", plotDir, medName)

pdf(pcaPlotPath_noRemRXNs, height = 4, width = 4.5)
pcaPlot(fba_noRemRXNs_PCA)
dev.off()


pcaPlotPath_h2o2 <- sprintf("%sfba_pca_%s_h2o2.pdf", plotDir, medName)

pdf(pcaPlotPath_h2o2, height = 4, width = 4.5)
pcaPlot(fba_h2o2_PCA)
dev.off()


pcaPlotPath_lowFe <- sprintf("%sfba_pca_%s_lowFe.pdf", plotDir, medName)

pdf(pcaPlotPath_lowFe, height = 4, width = 4.5)
pcaPlot(fba_lowFe_PCA)
dev.off()


pcaPlotPath_starv <- sprintf("%sfba_pca_%s_starv.pdf", plotDir, medName)

pdf(pcaPlotPath_starv, height = 4, width = 4.5)
pcaPlot(fba_starv_PCA)
dev.off()


fbaLowFe$linGroup <- gsub("[[:digit:]]", "", rownames(fbaLowFe))

ropls::opls(fbaLowFe[, 1:(nrow(fbaLowFe)-1)], fbaLowFe$linGroup, orthoI = NA, predI = 1)
