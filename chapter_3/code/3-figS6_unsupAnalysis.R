####################################################################################################
####################################################################################################
#                                                                                                  #
# Do the unsupervised analysises (PCA and HCA) of the CCMN normalized metabolomic data.            #
#                                                                                                  #
####################################################################################################
####################################################################################################
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(gplots)) install.packages("gplots")
library(gplots)
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(ggrepel)) install.packages("ggrepel")
library(ggrepel)
if(!require(dendextend)) install.packages("dendextend")
library(dendextend)

# Directory paths
####################################################################################################
rootDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/PA_paper_scripts/" ## Change this 
                                                                                      ## path 
                                                                                      ## according 
                                                                                      ## user's 
                                                                                      ## direcory
                                                                                      ## location.
dataDir <- paste0(rootDir, "data/")
outDir <- paste0(rootDir, "results/")
plotDir <- paste0(rootDir, "plots/")
figDir <- paste0(plotDir, "figS6/")

# Create directories if necessary
if(!dir.exists(plotDir)){
        dir.create(plotDir)
}
if(!dir.exists(outDir)){
        dir.create(outDir)
}
if(!dir.exists(figDir)){
        dir.create(figDir)
}

####################################################################################################
#                                                                                                  #
# Load the data.                                                                                   #
#                                                                                                  #
####################################################################################################

# Dictionary of metabolite names between different metabolite sets
dictionary <- read.csv(file = paste0(outDir, "dictionary_solvedPut.csv"), 
                       row.names = 1, 
                       stringsAsFactors = F)

# CCMN-normalized metabolomic data (the .tsv file, to avoid data loss).
mets_ccmn <- read.table(file = paste0(outDir, "mets_ccmn.tsv"), sep = "\t")
# When tsv is imported R changes the column names. Correct it with dictionary
colnames(mets_ccmn) <- dictionary$Consensus[match(colnames(mets_ccmn),
                                                  make.names(dictionary$Consensus))]

# Rhamnolipid production data
rhamnMat <- read.csv(file = paste0(dataDir, "rhamnMat.csv"),
                     row.names = 1, 
                     stringsAsFactors = F)

####################################################################################################
#                                                                                                  #
# Do HCA.                                                                                          #
#                                                                                                  #
####################################################################################################
# Remove sample 25 (H47921), as it's very different to the other replicates. 
mets_ccmn <- mets_ccmn[rownames(mets_ccmn) != "H47921_25", ]

rhamn3CatCols <- rhamnMat$rhamn3catsCols[match(gsub("\\_.*", "", rownames(mets_ccmn)), 
                                               rhamnMat$strains)]


# Obtain a dendrogram with putative metabolites to be able to extract the column dendrograme
# and plot just the non putative metabolites without altering the clustering

# Scale and center the data before clustering 
mets_ccmn_scaled <- apply(mets_ccmn, 2, function(x) (x-mean(x))/sd(x))
 
allMetsHeatMap <- heatmap.2(as.matrix(t(mets_ccmn_scaled)), 
                            distfun = function(x) dist(x, method = "euclidean"), 
                            density.info = "none", 
                            hclust = function(x) hclust(x, method = "ward.D"), 
                            dendrogram = "both", 
                            col = redgreen(75), breaks = 76, ColSideColors = rhamn3CatCols, 
                            notecol = NULL, trace = "none", xlab = "Strains", 
                            ylab = "Metabolites", 
                            margins = c(10, 16), 
                            cex.main = 20,
                            keysize = 0.7,
                            cexRow = 0.7,
                            cexCol = 1.2,
                            scale = "none",
                            notecex = 0.7,
                            key.xtickfun=function() {
                                    cex <- par("cex")*par("cex.axis")
                                    side <- 1
                                    line <- 0
                                    col <- par("col.axis")
                                    font <- par("font.axis")
                                    mtext("low", side=side, at=0, adj=0,
                                          line=line, cex=cex, col=col, font=font)
                                    mtext("high", side=side, at=1, adj=1,
                                          line=line, cex=cex, col=col, font=font)
                                    return(list(labels=FALSE, tick=FALSE))
                            })

# Extract dendrograms from past heatmap to keep branch line thick and to keep same 
# clustering for figure, which won't have putative compounds shown. 
colThick <- dendextend::set(allMetsHeatMap$colDendrogram, "branches_lwd", 1)
rowThick <- dendextend::set(allMetsHeatMap$rowDendrogram, "branches_lwd", 1)

# Remove putative compounds and change names to the Kyu Rhee ones. 

mets_ccmn_scaled_noPut <- mets_ccmn_scaled[, -grep("_?",
                                                   colnames(mets_ccmn_scaled), 
                                                   fixed = T)]

colnames(mets_ccmn_scaled_noPut) <- dictionary$KyuRheeNames[match(colnames(mets_ccmn_scaled_noPut), 
                                                                  dictionary$Consensus)]

pdf(paste0(figDir, "figS6_A.pdf"), width = 13, height = 11, pointsize = 2.5)
heatmap.2(as.matrix(t(mets_ccmn_scaled_noPut)), 
          distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), 
          dendrogram = "both", 
          col = redgreen(75), 
          breaks = 76, 
          ColSideColors = rhamn3CatCols, 
          notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", 
          cex.axis = 300,
          Colv = colThick,
          margins = c(13, 42), 
          cex.main = 20,
          keysize = 0.5,
          cexRow = 2.4,
          cexCol = 2.1,
          key = T,
          key.title = "Metabolite Levels",
          key.par=list(mgp=c(2.0, 0.5, 0),
                       mar=c(3.5, 0.5, 8, 2)),
          key.xtickfun=function() {
                  cex <- par("cex")*2*par("cex.axis")
                  side <- 1
                  line <- 1
                  col <- par("col.axis")
                  font <- par("font.axis")
                  mtext("low", side=side, at=0, adj=0,
                        line=line, cex=cex, col=col, font=font)
                  mtext("high", side=side, at=1, adj=1,
                        line=line, cex=cex, col=col, font=font)
                  return(list(labels=F, tick=F))
          }
)
dev.off()

####################################################################################################
#                                                                                                  #
# Do PCA.                                                                                          #
#                                                                                                  #
####################################################################################################

# Compute median of the strain replicates.

strains <- unique(gsub("\\_.*|(PA14).*", rep = "\\1", rownames(mets_ccmn)))

mets_ccmn_med <- data.frame(matrix(ncol = ncol(mets_ccmn),
                                   nrow = 0,
                                   dimnames = list(NULL,
                                                   colnames(mets_ccmn))))

for(i in seq_along(strains)){
        strain <- strains[i]
        strainMat <- mets_ccmn[grep(strain, rownames(mets_ccmn)), ]
        toRBind <- apply(strainMat, 2, median)
        mets_ccmn_med <- rbind.data.frame(mets_ccmn_med, toRBind)
}
rownames(mets_ccmn_med) <- strains
colnames(mets_ccmn_med) <- colnames(mets_ccmn)

# Remove putative metabolites and change names to definitive ones
mets_ccmn_med <- mets_ccmn_med[, -grep("_?", colnames(mets_ccmn_med), fixed = T)]
colnames(mets_ccmn_med) <- dictionary$KyuRheeNames[match(colnames(mets_ccmn_med), 
                                                         dictionary$Consensus)]


mets_ccmn_med_pca <- prcomp(apply(mets_ccmn_med, 2, function(x) (x - mean(x))/sd(x)),
                            scale. = F,
                            center = F)

# Define a function for obtaining the top N contributor metabolites
# to the specified components (in our case PC1 and PC2)
getTopContrib <- function(PC, topN = 12, x = "PC1", y = "PC2"){
        contrib <- facto_summarize(PC, 
                                   "var", 
                                   axes = c(as.numeric(gsub("PC", "", x)),
                                            as.numeric(gsub("PC", "", y))))
        contrib <- contrib[order(contrib$contrib, decreasing = T), 
                           c("name", "contrib")]
        topContrib <- as.character(contrib$name[1:topN])
        return(topContrib)
}

# Get top contributor compounds
topCont_ccmn_med_PC1_PC2 <- getTopContrib(mets_ccmn_med_pca,
                                          topN = 11,
                                          x = "PC1",
                                          y = "PC2")

# Define a function for obtaining the plot
pcBiplot <- function(PC, x="PC1", y="PC2", varPlotFilt = NULL){
        x <- "PC1"
        y <- "PC2"
        data <- data.frame(obsnames=row.names(PC$x), PC$x)
        data <- data[, c("obsnames", x, y)]
        data$rhamn3Cats <- rhamnMat$rhamn3cats[match(rownames(data), 
                                                     gsub("(PA14).*", 
                                                          rep ="\\1",
                                                          rhamnMat$strains))]
        rhamn3CatsDict <- c(0:2)
        names(rhamn3CatsDict) <- c("Rhamnolipid -", 
                                   "Rhamnolipid +/-",
                                   "Rhamnolipid +")
        data$rhamn3Cats <- names(rhamn3CatsDict)[match(data$rhamn3Cats, 
                                                       rhamn3CatsDict)]
        data$rhamn3Cats <- factor(data$rhamn3Cats, 
                                  levels = c("Rhamnolipid +", 
                                             "Rhamnolipid +/-",
                                              "Rhamnolipid -"))
        propVar <- summary(PC)$importance[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)

        plot <- ggplot(data, aes(x = data[, x], 
                                 y = data[, y], 
                                 label = obsnames, 
                                 color = rhamn3Cats)) + 
                scale_discrete_manual("Rhamnolipid production",
                                      aesthetics = "colour",
                                      values = unique(rhamnMat$rhamn3catsCols[order(rhamnMat$rhamn3cats,
                                                                                    decreasing = T)])) +
                geom_hline(yintercept = 0, alpha = 0.6) +
                geom_vline(xintercept = 0, alpha = 0.6) +
                geom_point() + 
                xlab(sprintf("%s (%s%%)", x, propX)) +
                ylab(sprintf("%s (%s%%)", y, propY)) +
                geom_text_repel() +
                #theme_minimal()
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        
        datapc <- data.frame(varnames=rownames(PC$rotation), 
                             PC$rotation)
        mult <- min(
                (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
                (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
        )
        datapc <- transform(datapc,
                            v1 = .7 * mult * (get(x)),
                            v2 = .7 * mult * (get(y))
        )
        datapc$x0 <- rep(0, nrow(datapc))
        datapc$y0 <- rep(0, nrow(datapc))
        if(!is.null(varPlotFilt)){
                datapc <- datapc[datapc$varnames %in% varPlotFilt, ]
        }
        plot <- plot +
                geom_text_repel(data=datapc, 
                                aes(x=v1, y=v2, label=varnames), 
                                color = "black", 
                                size = 3) + 
                geom_segment(data = datapc, aes(x=x0, 
                                                y=y0, 
                                                xend=v1, 
                                                yend=v2, 
                                                label = varnames),
                             arrow = arrow(length=unit(0.2,"cm"),
                                           type = "closed",
                                           angle = 20), 
                             alpha=0.75, 
                             color="black", 
                             size = 0.5)
        plot
}


pcBiplot(mets_ccmn_med_pca, varPlotFilt = topCont_ccmn_med_PC1_PC2)
ggsave(filename = paste0(figDir, "figS6_B.pdf"), height = 5, width = 7)


# Write a .tsv file of the ccmn normalized metabolites with putative compounds removed 
# and final names (Kyu Rhee's) in it for further analysis. 
####################################################################################################

mets_ccmn_noPut <- mets_ccmn[, -grep("_?", colnames(mets_ccmn), fixed = T)]


colnames(mets_ccmn_noPut) <- dictionary$KyuRheeNames[match(colnames(mets_ccmn_noPut), 
                                                           dictionary$Consensus)]

write.table(format(mets_ccmn_noPut, digits = 22), 
            file = paste0(outDir, "mets_ccmn_noPut.tsv"), 
            quote = F, 
            row.names = T, 
            sep = "\t")
