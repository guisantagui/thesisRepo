#########################################################################################################
#                                                                                                       #
# Plot the results of individual SNPs random forest.                                                    #
#                                                                                                       #
#########################################################################################################

if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(ggplot2)) install.packages("ggplot2", repos='http://cran.us.r-project.org')
library(ggplot2)
if(!require(gplots)) install.packages("gplots", repos='http://cran.us.r-project.org')
library(gplots)

# Create parsed, which would accept arguments from terminal and pass to the script
parser <- arg_parser("This script does random forest separating animal and human lineages")

parser <- add_argument(parser = parser,
                       arg = c("--deletSNPsMat", 
                               "--varImpFile", 
                               "--top",
                               "--provFile",
                               "--output"),
                       help = c("CSV file with proportion of strains within each lineage with each potentially deletereous SNP",
                                "Variable importance file from RF", 
                                "Number of top variables to plot",
                                "File with provean scores for each potentially deletereous SNP",
                                "Output directory to store the plots"),
                       flag = c(F, F, F, F, F),
                       default = list("--deletSNPsMat" = "/storage/PGO/results/mtb/snpsPipe/deletSNPsMat.csv",
                                      "--varImpFile" = "/storage/PGO/results/mtb/snpsPipe/RF/byStrain/varImpdeletSNPs_byStrainDF.csv",
                                      "--top" = 40,
                                      "--provFile" = "/storage/PGO/results/mtb/snpsPipe/provean/proveanSign.csv",
                                      "--output" = "/storage/PGO/results/mtb/snpsPipe/RF/byStrain/"))

parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################
delSNPsMatFile <- parsed$deletSNPsMat
varImpFile <- parsed$varImpFile
provFile <- parsed$provFile
outDir <- parsed$output
top <- parsed$top

# Load and parse the data
#########################################################################################################

deletSNPs <- read.csv(delSNPsMatFile, row.names = 1, check.names = F)
varImp <- read.csv(varImpFile, row.names = 1)
prov <- read.csv(provFile, row.names = 1)

# Prepare dataframe for heatmap
#########################################################################################################

varImp <- varImp[order(varImp$Overall, decreasing = T), ]
varImp <- varImp[1:top, ]
                                    
HM_DF <- data.frame(matrix(ncol = nrow(deletSNPs),
                           nrow = 0,
                           dimnames = list(NULL, 
                                           rownames(deletSNPs))))
                                    
print("Generating dataframe for building heatmap.")

varImpSNPs <- c()
for(i in seq_along(varImp$gene)){
        gene <- varImp$gene[i]
        snpsInGene <- colnames(deletSNPs)[gsub("\\..*", "", colnames(deletSNPs)) == gene]
        subMat <- deletSNPs[, snpsInGene]
        if(is.null(dim(subMat))){
                subMat <- data.frame(matrix(subMat, 
                                            ncol = length(subMat),
                                            nrow = 1,
                                            dimnames = list(snpsInGene,
                                                            colnames(HM_DF))))
        }else{
                subMat <- as.data.frame(t(subMat))
                colnames(subMat) <- colnames(HM_DF)
                rownames(subMat) <- snpsInGene
        }
        HM_DF <- rbind.data.frame(HM_DF, subMat)
}

write.csv(HM_DF, file = sprintf("%stop%sSNPProp_preFilt.csv", outDir, top))

HM_DF <- HM_DF[, c("L8", 
                   "L3",
                   "L2", 
                   "L4", 
                   "L7", 
                   "L1", 
                   "A3", 
                   "A4", 
                   "A2", 
                   "L6", 
                   "L9", 
                   "A1", 
                   "L5")]

# Get rid of SNPs that are almost absent (below 5% of the isolates within lineage) in all the lineages
filtVec <- apply(HM_DF[, colnames(HM_DF) != "gene"] , 1, function(x) any(x > 0.2))

HM_DF <- HM_DF[filtVec, ] 


write.csv(HM_DF, file = sprintf("%stop%sSNPProp.csv", outDir, top))

pdf(file = sprintf("%stop%sRF_indStrains_deletSNPs_hm.pdf", outDir, as.character(top)), height = 11, width = 4)
heatmap.2(as.matrix(HM_DF), 
          dendrogram = "none",
          Rowv = F, 
          Colv = F, 
          col = greenred(75), 
          breaks = 76, 
          trace = "none", 
          key = F,
          margins = c(4, 7),
          cexRow = 0.85,
          cexCol = 1)
dev.off()

print(sprintf("top%sRF_indStrains_deletSNPs_hm.pdf saved at %s.", as.character(top), outDir))

# Do variable importance and provean score barplots
#########################################################################################################
doBarPlots <- function(HMDF, vImpDF, provDF){
        HMDF <- HMDF[1:top, ]
        HMDF_genes <- gsub("\\..*", "", rownames(HMDF))
        variant <- gsub(".*\\.", "", rownames(HMDF))
        plotDF <- data.frame(gene = HMDF_genes,
                             variant = variant,
                             importance = vImpDF$Overall[match(HMDF_genes, vImpDF$gene)],
                             provean = provDF$PROVEAN[match(paste(HMDF_genes, 
                                                                  variant, sep = "_"),
                                                            paste(provDF$LOCUS, 
                                                                  provDF$AA_one, sep = "_"))])
        plotDF$provean[is.na(plotDF$provean)] <- 0
        plotDF <- plotDF[nrow(plotDF):1, ]
        
        plotVImp <- ggplot(data = plotDF, mapping = aes(x = importance, 
                                            y = factor(paste(gene, variant, sep = "_"), 
                                                       levels = paste(gene, variant, sep = "_")))) +
                geom_bar(stat = "identity", fill = "steelblue") + 
                labs(x = "Variable importance", y = "Variant") + 
                theme_minimal() + 
                theme(axis.text.y = element_text(hjust = 0))
        plotProv <- ggplot(data = plotDF, mapping = aes(x = provean, 
                                                        y = factor(paste(gene, variant, sep = "_"), 
                                                                   levels = paste(gene, variant, sep = "_")))) +
                geom_bar(stat = "identity", fill = "steelblue") + 
                labs(x = "PROVEAN score", y = "Variant") + 
                theme(axis.text.y = element_text(hjust = 0))
        out <- list(vImpPlt = plotVImp, 
                    provPlt = plotProv)
        return(out)
}

barPlots <- doBarPlots(HM_DF, varImp, prov)

barPlots$vImpPlt
ggsave(paste0(outDir, "deletSnp_byStrain_vImp.pdf"), height = 10, width = 3)
print(sprintf("deletSnp_byStrain_vImp.pdf saved at %s.", outDir))

barPlots$provPlt
ggsave(paste0(outDir, "deletSnp_byStrain_prov.pdf"), height = 10, width = 3)
print(sprintf("deletSnp_byStrain_prov.pdf saved at %s.", outDir))