#########################################################################################################
#                                                                                                       #
# Run Correspondence Analysis on binaryzed deletion data.                                               #
#                                                                                                       #
#########################################################################################################

if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(FactoMineR)) install.packages("FactoMineR")
library(FactoMineR)
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(stringr)) install.packages("stringr", repos='http://cran.us.r-project.org')
library(stringr)

# Create parser, which accepts arguments from terminal and pass to the script
parser <- arg_parser("This script runs Correspondence analysis on the binarized set of the input deletion data and, if selected, on the sets of just enzymatic genes and/or genes included in iEK1011 Genome Scale Metabolic Model.")

parser <- add_argument(parser = parser,
                       arg = c("input", 
                               "--output",
                               "--iEKGenes"),
                       help = c("Input deletion file",
                                "Output directory where results of the analysis will be stored",
                                "Path to .txt file with list of genes included in iEK1011 Genome Scale Metabolic Model"),
                       flag = c(F, F, F),
                       default = list("input" = "/storage/PGO/results/mtb/deletions_allGenes/allDels.csv",
                                      "--output" = "/storage/PGO/results/mtb/deletions_analysis/",
                                      "--iEKGenes" = "/storage/PGO/data/mtb/GSMMs/iEK1011_2.0/iEK1011_2.0_genes.txt"))
                                      
parser <- add_argument(parser = parser, 
                       arg = c("--enzs", "--iEK", "--plots", "--geneLabl"),
                       help = c("Do the Correspondence Analysis on a reduced gene set of just enzymatic genes", 
                                "Do the Correspondence Analysis on a reduced gene set including just genes in iEK1011 model",
                                "Plot the result of the analyses", 
                                "Write gene labels in plots"),
                       flag = c(T, T, T, T))
                                      
parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################

# Use as input allDels.csv and then change to the Bin name because in the pipeline this file will be 
# the input for the whole pipeline. 
inputBin <- paste0(gsub(".csv", "", parsed$input), "Bin")
inputDir <- parsed$output
inputBin <- basename(inputBin)
outDir <- parsed$output

# Load and parse the data
#########################################################################################################

# Binarized deletion data
print(sprintf("Loading %s.csv from %s...", inputBin, inputDir))
allDelsBin <- as.data.frame(read.csv(sprintf("%s%s.csv", inputDir, inputBin), row.names = 1))
rownames(allDelsBin) <- allDelsBin$gene
allDelsBin <- allDelsBin[, (sapply(allDelsBin, class) == "numeric" | sapply(allDelsBin, class) == "integer") & colnames(allDelsBin) != "ncbi_geneID"]
#allDelsBin <- allDelsBin[, colSums(allDelsBin) <= nrow(allDelsBin)]
for(i in 1:ncol(allDelsBin)){
        allDelsBin[, i][is.na(allDelsBin[, i])] <- 0
}
allDelsBin <- allDelsBin[, (colSums(allDelsBin) <= nrow(allDelsBin))]
print(sprintf("%s.csv loaded!", inputBin))

if(parsed$enzs == T){
        enzInputBin <- paste0("enz", 
                              str_match(inputBin, 
                                        "DelsBin")[, 1])
        print(sprintf("Loading %s.csv from %s...", enzInputBin, inputDir))
        enzDelsBin <- as.data.frame(read.csv(sprintf("%s%s.csv", inputDir, enzInputBin), row.names  = 1))
        rownames(enzDelsBin) <- enzDelsBin$gene
        enzDelsBin <- enzDelsBin[, sapply(enzDelsBin, class) == "integer"]
        for(i in 1:ncol(allDelsBin)){
                allDelsBin[, i][is.na(allDelsBin[, i])] <- 0
        }
        enzDelsBin <- enzDelsBin[, colSums(enzDelsBin) <= nrow(enzDelsBin)]
        print(sprintf("%s.csv loaded!", enzInputBin))
        #enzDelsBin[, 8:ncol(enzDelsBin)] <- apply(enzDelsBin[, 8:ncol(enzDelsBin)], 2, as.numeric)
}


#allDelsBin <- allDelsBin[, 2:ncol(allDelsBin)]
#enzDelsBin <- enzDelsBin[, 2:ncol(enzDelsBin)]

# load iEK1011 genes. 
iEK1011_genes <- read.table(parsed$iEKGenes, sep = "\t")$V1

print(head(allDelsBin$gene))
print(head(colnames(allDelsBin)))

# Parse data into counts dataframes for doing Correspondence Analysis)
mkCADF <- function(delsBin){
        genes <- rownames(delsBin)
        print(head(genes))
        #delsBin <- delsBin[, (sapply(delsBin, class) == "numeric" | sapply(delsBin, class) == "integer") & colnames(delsBin) != "ncbi_geneID"]
        lins <- unique(sapply(colnames(delsBin), function(x) strsplit(x, "_")[[1]][2]))
        CADF <- data.frame(matrix(ncol = length(lins),
                                  nrow = nrow(delsBin),
                                  dimnames = list(rownames(delsBin),
                                                  lins)))
        for(l in lins){
                CADF[, l] <- rowSums(delsBin[, grep(paste0("_", l), 
                                                    colnames(delsBin))])
        }
        rownames(CADF) <- genes
        return(CADF)
}

allDelsCorrDF <- mkCADF(allDelsBin)
print(head(allDelsCorrDF))
write.csv(allDelsCorrDF, file = paste0(outDir, "allDelsCorrDF.csv"))

# Obtain enzymatic deletion count matrix in case the flag is present.                                                          
if(parsed$enzs == T){
        enzDelsCorrDF <- mkCADF(enzDelsBin)
        print(head(enzDelsCorrDF))
}

# Filter to get just iEK1011 genes in case the flag is present.
if(parsed$iEK == T){
        #iEKDelsCorrDF <- allDelsCorrDF[match(iEK1011_genes, allDelsBin$gene), ]
        iEKDelsCorrDF <- allDelsCorrDF[rownames(allDelsCorrDF) %in% iEK1011_genes, ]
        iEKDelsCorrDF <- iEKDelsCorrDF[!is.na(iEKDelsCorrDF$L1), ]
        print((iEK1011_genes))
        print(head(rownames(allDelsCorrDF)))
        print(head(iEKDelsCorrDF))
}

# Run Correspondence analysis and do plots, if flag present.
#########################################################################################################
outName <- gsub("Bin", "", inputBin)
print(sprintf("Running Correspondence Analysis of %s...", inputBin))
allDelsCA <- CA(allDelsCorrDF, graph = FALSE)
save(allDelsCA, file = sprintf("%s%sCA.RData", outDir, outName))
print(sprintf("%s saved at %s.", outName, outDir))
if(parsed$plots == T){
        pdf(file = sprintf("%s%sCA.pdf", outDir, outName))
        fviz_ca_biplot(allDelsCA, repel = TRUE, label = "col")
        dev.off()
        print(sprintf("%sCA.pdf saved at %s.", outName, outDir))
}
if(parsed$plots == T & parsed$geneLabl == T){
        pdf(file = sprintf("%s%sCA_wGenLabl.pdf", outDir, outName))
        fviz_ca_biplot(allDelsCA, repel = FALSE, label = "all")
        dev.off()
        print(sprintf("%sCA_wGenLabl.pdf saved at %s.", outName, outDir))
}
if(parsed$enzs == T){
        enzOutName <- gsub("Bin", "", enzInputBin)
        print(sprintf("Running Correspondence Analysis of %s...", enzInputBin))
        enzDelsCA <- CA(enzDelsCorrDF, graph = FALSE)
        save(enzDelsCA, file = sprintf("%s%sCA.RData", outDir, enzOutName))
        print(sprintf("%s saved at %s.", enzOutName, outDir))
        if(parsed$plots == T){
                pdf(file = sprintf("%s%sCA.pdf", outDir, enzOutName))
                fviz_ca_biplot(enzDelsCA, repel = TRUE, label = "col")
                dev.off()
                print(sprintf("%sCA.pdf saved at %s.", enzOutName, outDir))
        }
        if(parsed$plots == T & parsed$geneLabl){
                pdf(file = sprintf("%s%sCA_wGenLabl.pdf", outDir, enzOutName))
                fviz_ca_biplot(enzDelsCA, repel = FALSE, label = "all")
                dev.off()
                print(sprintf("%sCA_wGenLabl.pdf saved at %s.", enzOutName, outDir))
        }
}
if(parsed$iEK == T){
        iEKOutName <- paste0("iEK", str_match(outName, "Dels"))
        iEKDelsCA <- CA(iEKDelsCorrDF, graph = FALSE)
        save(iEKDelsCA, file = sprintf("%s%sCA.RData", outDir, iEKOutName))
        print(sprintf("%s saved at %s.", iEKOutName, outDir))
        if(parsed$plots == T){
                pdf(file = sprintf("%s%sCA.pdf", outDir, iEKOutName))
                fviz_ca_biplot(iEKDelsCA, repel = TRUE, label = "col")
                dev.off()
                print(sprintf("%s.pdf saved at %s.", iEKOutName, outDir))
        }
        if(parsed$plots == T & parsed$geneLabl == T){
                pdf(file = sprintf("%s%sCA_wGenLabl.pdf", outDir, iEKOutName))
                fviz_ca_biplot(iEKDelsCA, repel = FALSE, label = "all")
                dev.off()
                print(sprintf("%sCA_wGenLabl.pdf saved at %s.", iEKOutName, outDir))
        }
}

print("All Correspondence Analyses done!")