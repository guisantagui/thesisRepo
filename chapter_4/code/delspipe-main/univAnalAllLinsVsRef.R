#########################################################################################################
#                                                                                                       #
# This script does univariate comparisons ef each lineage to reference lineage (default is L4, as H37Rv #
# belongs to it) with the aim of identifying what are the differences of each lineage to H37Rv (L4). It #
# stores the genes with a significant adjusted p-value (default < 0.05) in .txt files in a subdirectory #
# within the output folder.                                                                             #
#                                                                                                       #
#########################################################################################################
if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)

# Create parser, which accepts arguments from terminal and pass to the script
parser <- arg_parser("This script does univariate comparisons ef each lineage to reference lineage (default is L4, as H37Rv belongs to it) with the aim of identifying what are the differences of each lineage to reference strain. It stores the genes with a significant adjusted p-value (default < 0.05) in .txt files in a subdirectory within the output folder.")

parser <- add_argument(parser = parser,
                       arg = c("input", 
                               "--output",
                               "--ref",
                               "--alpha"),
                       help = c("Input deletion file",
                                "Output directory where results of the analysis will be stored",
                                "Reference lineage to do the comparisons versus it", 
                                "Alpha value for filtering the significative p-values"),
                       flag = c(F, F, F, F),
                       default = list("input" = "/storage/PGO/results/mtb/deletions_allGenes/allDels.csv",
                                      "--output" = "/storage/PGO/results/mtb/deletions_analysis/",
                                      "--ref" = "L4",
                                      "--alpha" = 0.05))
                                      
parser <- add_argument(parser = parser, 
                       arg = "--chckIfUnivDone",
                       help = "Check before running univariate analysis if it's already done",
                       flag = T)
                                      
parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################
ref <- parsed$ref
outDir <- parsed$output
# Create a subdirectory, if it doesn't exist, for storing the results of the filtered p-values.
filtOutDir <- paste0(parsed$output, sprintf("linVs%sDiffGenes/", ref))
if(!dir.exists(filtOutDir)){
        dir.create(filtOutDir)
}
resuDir <- "/storage/PGO/results/mtb/deletions_analysis/"

# Load the data and do the dataframe of reference lineage to do the comparisons to 
#########################################################################################################
allDels <- read.csv(parsed$input)
allDelsBin <- read.csv(paste0(outDir, gsub(".csv", "Bin.csv", basename(parsed$input))))

colnames(allDels)[1] <- "gene"
colnames(allDelsBin)[1] <- "gene"

allDels_ref <- allDels[, grep(ref, colnames(allDels))]
allDelsBin_ref <- allDelsBin[, grep(ref, colnames(allDelsBin))]

rownames(allDels_ref) <- allDels$gene
rownames(allDelsBin_ref) <- allDelsBin$gene

# Do the univariate tests
#########################################################################################################
linsToCompare <- unique(gsub(".*_", "", colnames(allDels[, sapply(allDels, class) == "numeric"])))
linsToCompare <- linsToCompare[linsToCompare  != ref]


outName <- sprintf("linVs%sPValsDF.csv", ref)

if(file.exists(sprintf("%s%s", outDir, outName)) & parsed$chckIfUnivDone == T){
        print(sprintf("Loading %s from %s...", outName, outDir))
        linVsRefPValsDF <- read.csv(file = sprintf("%s%s", outDir, outName), row.names = 1)
        #colnames(linVsRefPValsDF) <- paste0(colnames(linVsRefPValsDF), 
        #                                    c("manWPVal", 
        #                                      "manWPValAdj", 
        #                                      "fishPVal", 
        #                                      "fishPValAdj", 
        #                                      "chiSPVal", 
        #                                      "chiSPValAdj")) # Just this time
}else{
        print(sprintf("Doing univariate tests of %s versus the rest of the lineages...", ref))
        linVsRefPVals <- list()
        for(i in seq_along(linsToCompare)){
                lin <- linsToCompare[i]
                print(sprintf("Computing univariate tests of %s versus %s (reference)...", lin, ref))
                MWPVals <- c()
                fiPVals <- c()
                CSPVals <- c()
                for(j in seq_along(allDels$gene)){
                        allDels_refVec <- as.numeric(allDels_ref[j, ])
                        allDelsBin_refVec <-  as.numeric(allDelsBin_ref[j, ])
                        allDels_linVec <- as.numeric(allDels[j, grep(lin, colnames(allDels))])
                        allDelsBin_linVec <- as.numeric(allDelsBin[j, grep(lin, colnames(allDelsBin))])
                        contMat <- matrix(c(sum(allDelsBin_refVec), 
                                            (length(allDelsBin_refVec) - sum(allDelsBin_refVec)),
                                            sum(allDelsBin_linVec), 
                                            (length(allDelsBin_linVec) - sum(allDelsBin_linVec))),
                                          ncol = 2,
                                          nrow = 2, 
                                          dimnames = list(c("mut", "int"), 
                                                          c(ref, lin))) 
                        manWPVal <- wilcox.test(allDels_refVec, allDels_linVec)$p.value
                        fishPVal <- fisher.test(contMat)$p.value
                        chiSPVal <- chisq.test(contMat)$p.value
                        MWPVals <- c(MWPVals, manWPVal)
                        fiPVals <- c(fiPVals, fishPVal)
                        CSPVals <- c(CSPVals, chiSPVal)
                }
                pValsDF <- data.frame(manWPVal = MWPVals, 
                                      manWPValAdj = p.adjust(MWPVals, "BH"),
                                      fishPVal = fiPVals, 
                                      fishPValAdj = p.adjust(fiPVals, "BH"),
                                      chiSPVal = CSPVals,
                                      chiSPValAdj = p.adjust(CSPVals, "BH"))
                rownames(pValsDF) <- allDels$gene
                linVsRefPVals[[i]] <- pValsDF
        }
        names(linVsRefPVals) <- linsToCompare

        linVsRefPValsDF <- do.call("cbind.data.frame", linVsRefPVals)

        colnames(linVsRefPValsDF) <- paste(rep(names(linVsRefPVals),
                                               each = 6),
                                           rep(colnames(linVsRefPVals[[1]]), 6),
                                           sep = "_")
                                   
        rownames(linVsRefPValsDF) <- allDels$gene
}

# Export the results.
#########################################################################################################
write.csv(linVsRefPValsDF, file = paste0(outDir, outName))
print(sprintf("%s saved at %s.", outName, outDir))

# Filter the significative p-values and store them in the filtOutDir.
#########################################################################################################
lins <- linsToCompare
alpha <- parsed$alpha
signGeneLins <- list()
for(i in seq_along(lins)){
        lin <- lins[i]
        subMat <- linVsRefPValsDF[, grep(lin, colnames(linVsRefPValsDF))]
        subMat <- subMat[, grep("Adj", colnames(subMat))]
        sgnfMW <- rownames(linVsRefPValsDF)[subMat[, grep("manWPValAdj", colnames(subMat))] < alpha]
        sgnfMW <- sgnfMW[!is.na(sgnfMW)]
        sgnfFi <- rownames(linVsRefPValsDF)[subMat[, grep("fishPValAdj", colnames(subMat))] < alpha]
        sgnfFi <- sgnfFi[!is.na(sgnfFi)]
        sgnfCS <- rownames(linVsRefPValsDF)[subMat[, grep("chiSPValAdj", colnames(subMat))] < alpha]
        sgnfCS <- sgnfCS[!is.na(sgnfCS)]
        signLst <- list(MW = sgnfMW,
                        Fi = sgnfFi,
                        CS = sgnfCS)
        signGeneLins[[i]] <- signLst
        for(j in seq_along(signLst)){
                fileName <- paste(filtOutDir, 
                                  lin,
                                  names(signLst)[j],
                                  "Sign.txt",
                                  sep = "")
                write.table(signLst[[j]],
                            file = fileName,
                            col.names = F,
                            row.names = F,
                            quote = F)
                print(sprintf("%s saved at %s.", basename(fileName), dirname(fileName)))
        }
}
names(signGeneLins) <- lins
