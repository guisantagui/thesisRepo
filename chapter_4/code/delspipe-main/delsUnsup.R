#########################################################################################################
#                                                                                                       #
# This script runs several unsupervised analyses on the input deletion data.                            #
#                                                                                                       #
#########################################################################################################

if(!require(dplyr)) install.packages("dplyr", repos='http://cran.us.r-project.org')
library(dplyr)
if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(stringr)) install.packages("stringr", repos='http://cran.us.r-project.org')
library(stringr)
if(!require(factoextra)) install.packages("factoextra", repos = "https://cran.r-project.org")
library(factoextra)
if(!require(ggplot2)) install.packages("ggplot2", repos = "https://cran.r-project.org")
library(ggplot2)
if(!require(ape)) install.packages("ape", repos = "https://cran.r-project.org")
library(ape)

# Create parser, which accepts arguments from terminal and pass to the script
parser <- arg_parser("This script generates runs Principal Component Analysis and Hierarchical Clustering Analysis on the input deletion data and, if indicated, on the set of enzyme-coding genes.")

parser <- add_argument(parser = parser,
                       arg = c("input", 
                               "--output"),
                       help = c("Input deletion percentage file",
                                "Output directory where results of the analysis will be stored"),
                       flag = c(F, F),
                       default = list("input" = "/storage/PGO/results/mtb/deletions_allGenes/allDels.csv",
                                      "--output" = "/storage/PGO/results/mtb/deletions_analysis/"))
                                      
parser <- add_argument(parser = parser, 
                       arg = c("--enzs", "--plots"),
                       help = c("Do the analyses on a reduced gene set of just enzymatic genes", 
                                "Plot the result of the analyses"),
                       flag = c(T, T))
                                      
parsed <- parse_args(parser)

# Load functions
source("functions.R")

# Directory stuff
#########################################################################################################
input <- parsed$input
outDir <- parsed$output
outName <- gsub(".csv", "", basename(input))

# Load the data
#########################################################################################################
print("Loading deletion data...")
allDels <- read.csv(input)
#colnames(allDels)[1] <- "gene"
print(sprintf("%s.csv imported", outName))

print(head(colnames(allDels), 20))

# Run analyses.
#########################################################################################################
outPathPCA <- sprintf("%s%sPCA.RData", outDir, outName)
allDelsNum <- t(allDels[, sapply(allDels, class) == "numeric"])
if(!file.exists(outPathPCA)){
        print("Running PCA of deletion data...")
        allDelsPCA <- prcomp(allDelsNum)
        save(allDelsPCA, file = outPathPCA)
        print(sprintf("%sPCA.RData saved at %s.", outName, outDir))
}else{
        print(sprintf("PCA already computed, loading %sPCA.RData from %s...", outName, outDir))
        load(outPathPCA)
}

#outPathHCA <- sprintf("%s%sHCA.RData", outDir, outName)
#if(!file.exists(outPathHCA)){
#        print("Running HCA of deletion data...")
#        allDelsDist <- dist(allDelsNum, method = "euclidean")
#        allDelsHCA <- hclust(allDelsDist, method = "ward.D")
#        save(allDelsHCA, file = outPathHCA)
#        print(sprintf("%sHCA.RData saved at %s.", outName, outDir))
#}else{
#        print(sprintf("HCA already computed, loading %sHCA.RData from %s...", outName, outDir))
#        load(outPathHCA)
#}

# Plot stuff.
#########################################################################################################

# Do colorscale 
lineages <- sapply(rownames(allDelsPCA$x), function(x) strsplit(x, "_")[[1]][2])
linGroup <- gsub("[[:digit:]]+", "", lineages)
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

print(linCols)

if(parsed$plots == T){
        print("Plotting results...")
        # Plot PCA of all deletions 
        #pdf(file = sprintf("%s%sPCA.pdf", outDir, outName))
        fviz_pca_ind(allDelsPCA, 
                     col.ind = lineages,
                     habillage = lineages,
                     geom = "point") + 
        scale_color_manual(name = "Lineages", 
                           labels = lineages,
                           values = linCols) +
        scale_shape_manual(name = "Lineages", 
                           values = c(rep(2, 4),
                                      rep(19, 9)),
                           labels = sapply(linGroup, function(x) if(x == "A") x <- "Animal" else x <- "Human"))
        #dev.off()
        ggsave(sprintf("%s%sPCA.pdf", outDir, outName))
        print(sprintf("%sPCA.pdf saved at %s.", outName, outDir))
        # Plot dendrogram of all deletions 
        #pdf(sprintf("%s%sHCA.pdf", outDir, outName))
        #plot(as.phylo(allDelsHCA), 
        #     type = "fan",
        #     tip.color = dendCols,
        #     cex = 0.05,
        #     label.offset = 1)
        #dev.off()
        #print(sprintf("%sHCA.pdf saved at %s.", outName, outDir))
}

if(parsed$enzs == T){
        print("Obtaining enzymatic deletion data...")
        enzDels <- filtEnzDels(allDels)
        enzOutName <- paste0("enz", str_match(outName, "Dels")[, 1])
        enzOutPathPCA <- sprintf("%s%sPCA.RData", outDir, enzOutName)
        if(!file.exists(enzOutPathPCA)){
        print("Running PCA of enzymatic deletion data...")
                enzDelsNum <- t(enzDels[, sapply(enzDels, class) == "numeric"])
                enzDelsPCA <- prcomp(enzDelsNum)
                save(enzDelsPCA, file = enzOutPathPCA)
                print(sprintf("%sPCA.RData saved at %s.", enzOutName, outDir))
        }else{
                print(sprintf("PCA of enzymatic genes already computed, loading %sPCA.RData from %s...", enzOutName, outDir))
                load(enzOutPathPCA)
        }
        enzOutPathHCA <- sprintf("%s%sHCA.RData", outDir, enzOutName)
        #if(!file.exists(enzOutPathHCA)){
        #        print("Running HCA of enzymatic deletion data...")
        #        enzDelsDist <- dist(enzDelsNum, method = "euclidean")
        #        enzDelsHCA <- hclust(enzDelsDist, method = "ward.D")
        #        save(enzDelsHCA, file = enzOutPathHCA)
        #        print(sprintf("%sHCA.RData saved at %s.", enzOutName, outDir))
        #}else{
        #        print(sprintf("HCA of enzymatic genes already computed, loading %sPCA.RData from %s...", enzOutName, outDir))
        #        load(enzOutPathHCA)
        #}
        if(parsed$plots == T){
                print("Plotting results of enzymatic genes...")
                # Plot PCA of enzymatic deletions 
                #pdf(file = sprintf("%s%sPCA.pdf", outDir, enzOutName))
                fviz_pca_ind(enzDelsPCA, 
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
                #dev.off()
                ggsave(sprintf("%s%sPCA.pdf", outDir, enzOutName))
                print(sprintf("%sPCA.pdf saved at %s.", enzOutName, outDir))
                # Plot dendrogram of enzymatic deletions 
                #pdf(sprintf("%s%sHCA.pdf", outDir, enzOutName))
                #plot(as.phylo(enzDelsHCA), 
                #     type = "fan",
                #     tip.color = dendCols,
                #     cex = 0.05,
                #     label.offset = 1)
                #dev.off()
                #print(sprintf("%sHCA.pdf saved at %s.", enzOutName, outDir))
        }
}

print("PCA and HCA of non-binarized deletion data done!")