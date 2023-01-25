#########################################################################################################
#                                                                                                       #
# Plot the results of SNPs random forest                                                                #
#                                                                                                       #
#########################################################################################################

if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(ggplot2)) install.packages("ggplot2", repos='http://cran.us.r-project.org')
library(ggplot2)
if(!require(gplots)) install.packages("gplots", repos='http://cran.us.r-project.org')
library(gplots)
if(!require(rtracklayer)) install.packages("rtracklayer", repos='http://cran.us.r-project.org')
library(rtracklayer)

# Create parsed, which accepts arguments from terminal and pass to the script
parser <- arg_parser("This script plots results from Random Forest")

parser <- add_argument(parser = parser,
                       arg = c("--rfDir",
                               "--provean",
                               "--annot",
                               "--snpPath",
                               "--top",
                               "--output"),
                       help = c("Directory where RF results are located",
                                "Path to file with provean scores of each SNP.",
                                "Path to annotation GFF file",
                                "Path to deletereous SNP matrix file",
                                "Number of top variables to show in RF plots",
                                "Output directory where output files will be stored"),
                       flag = c(F, F, F, F, F, F))
                       
parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################

snpFile <- paste0(parsed$snpPath)
rfDir <- parsed$rfDir
provFile <- paste0(parsed$provean, "proveanScores.txt")
annoFile <- paste0(parsed$annot, "genes.gff")
outDir <- parsed$output


# Load data
#########################################################################################################

print(sprintf("Loading %s deletereous SNPs file...", basename(snpFile)))
deletSNPs <- read.csv(snpFile, row.names = 1, stringsAsFactors = T)
vImp <- read.csv(paste0(rfDir, "varImpdeletSNPsMatNoBinUpDF.csv"),
                 row.names = 1)
                 
provScores <- read.table(provFile, sep = "\t", header = T, row.names = 1)
anno <- as.data.frame(readGFF(annoFile))

# Do heatmap
#########################################################################################################
top <- parsed$top

vImp <- vImp[1:top, ]

HM_DF <- data.frame(matrix(ncol = nrow(deletSNPs),
                    nrow = 0,
                    dimnames = list(NULL, 
                                    rownames(deletSNPs))))

print("Generating dataframe for building heatmap.")
for(i in seq_along(vImp$gene)){
        gene <- vImp$gene
        subMat <- deletSNPs[, grep(gene, colnames(deletSNPs))]
        snpsInGene <- colnames(subMat)
        lins <- rownames(subMat)
        subMat <- as.data.frame(t(subMat))
        colnames(subMat) <- lins
        rownames(subMat) <- snpsInGene
        HM_DF <- rbind.data.frame(HM_DF, subMat)
}

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

write.csv(HM_DF, file = sprintf("%stop%sSNPProp.csv", outDir, top))
         
pdf(file = sprintf("%stop%sRFDeletSNPs_hm.pdf", outDir, as.character(top)), height = 10, width = 4)
heatmap.2(as.matrix(HM_DF), 
          Rowv = F, 
          Colv = F, 
          col = redgreen(75), 
          breaks = 76, 
          trace = "none", 
          key = F,
          margins = c(4, 4))
dev.off()

print(sprintf("top%sRFDeletSNPs_hm.pdf saved at %s.", as.character(top), outDir))