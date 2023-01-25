#########################################################################################################
#                                                                                                       #
# Take as input a vector of genes and gives as output dataframes with the proportion of strains within  #
# each lineage that have these genes deleted with above two given percentages.                          #
#                                                                                                       #
#########################################################################################################

if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)

# Create parsed, which would accept arguments from terminal and pass to the script
parser <- arg_parser("This script takes as input a txt file of genes and gives as output dataframes with the proportion of strains within each lineage that have these genes deleted with above a given percentage")

parser <- add_argument(parser = parser,
                       arg = c("input", 
                               "--output", 
                               "--low",
                               "--high",
                               "--plots"),
                       help = c("Input file",
                                "Output directory where output files will be stored", 
                                "Low threshold percentage of deletion to test",
                                "Low threshold percentage of deletion to test",
                                "Plot heatmap"),
                       flag = c(F, F, F, F, T))
                       
parser <- add_argument(parser = parser, 
                       arg = "--delFile",
                       help = "CSV file with deletion percentages. Default: allDels.csv",
                       default = "/storage/PGO/results/mtb/deletions_allGenes/allDels.csv")
                       
parser <- add_argument(parser = parser, 
                       arg = "--lineages",
                       help = "Select what lineages are desired to see the deleted proportion",
                       default = c("L1,L2,L3,L4,L5,L6,L7,L8,L9,A1,A2,A3,A4"))
                       
parser <- add_argument(parser = parser, 
                       arg = "--rowReord",
                       help = "Determine a order of the rows (lineages) different than default",
                       default = c("L1,L2,L3,L4,L5,L6,L7,L8,L9,A1,A2,A3,A4"))


parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################
resultsDir <- parsed$output

# Load data
#########################################################################################################
allDels <- read.csv(parsed$delFile)
colnames(allDels)[1] <- "gene"

# Obtain dataframe of >--low and >high of deletion proportion of the input genes
#########################################################################################################
getDelProp <- function(delObjkt, geneVec, linVec, thrhld = 5){
        #geneVec <- factor(geneVec)
        propVecList <- list()
        for(i in 1:length(linVec)){
                lin <- linVec[i]
                genePropVec <- c()
                for(j in 1:length(geneVec)){
                        gene <- geneVec[j]
                        geneProp <- sum(delObjkt[gene==delObjkt$gene, grep(lin, colnames(delObjkt))] > thrhld)/length(delObjkt[gene==delObjkt$gene, grep(lin, colnames(delObjkt))])
                        genePropVec <- c(genePropVec, geneProp)
                }
                names(genePropVec) <- geneVec
                propVecList[[i]] <-genePropVec
        }
        names(propVecList) <- linVec
        propVecDF <- data.frame(propVecList)
        rownames(propVecDF) <- geneVec
        return(propVecDF)
}

genes <- as.character(read.table(parsed$input)$V1)
#print(genes)
genes <- genes[genes %in% allDels$gene]

linVec <- strsplit(parsed$lineages, ",")[[1]]
linOrdVec <- strsplit(parsed$rowReord, ",")[[1]]

# print(linOrdVec)

low <- as.numeric(parsed$low)
high <- as.numeric(parsed$high)

delPropLow <- t(getDelProp(allDels, genes, linVec = linVec, thrhld = low))
delPropHigh <- t(getDelProp(allDels, genes, linVec = linVec, thrhld = high))

#print(delPropLow)
#print(delPropLow)

ordVec <- c()
for(i in 1:length(genes)){
        ordVec <- c(ordVec, i, (length(genes) + i))
}

delProp <- cbind.data.frame(delPropLow, delPropHigh)[, ordVec]

colnames(delProp) <- paste(gsub(".1", "", colnames(delProp), fixed = T), as.character(c(low, high)), sep = "_")


delProp <- delProp * 100


delProp <- delProp[match(linOrdVec, rownames(delProp)), ]

inputName <- gsub(".txt", "", basename(parsed$input))
delPropFile <- paste(inputName, "DelProp.csv", sep = "")

write.csv(delProp, file = paste(resultsDir, delPropFile, sep = ""))
print(paste(delPropFile, "saved at", resultsDir))

# Plot heatmaps. 
#########################################################################################################
if(parsed$plots == T){
        plotFile <- paste(inputName, "DelPropHM.pdf", sep = "")
        pdf(file = paste(resultsDir, plotFile, sep = ""))
                heatmap(as.matrix(delProp), #[nrow(delProp):1, ]), 
                        Rowv = NA, 
                        Colv = NA, 
                        scale = "none", 
                        col = rev(gplots::redgreen(nrow(delProp))),
                        labCol = paste(rep(genes, each = 2), 
                                       rep(c(paste(as.character(low), "%", sep = ""), 
                                             paste(as.character(high), "%", sep = "")),
                                           length(genes))),
                        margins = c(5, 5))
        dev.off()
        print(paste(plotFile, "saved at", resultsDir))
}