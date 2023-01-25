#########################################################################################################
#                                                                                                       #
# Run Correspondence Analysis on SNP proportion data.                                                   #
#                                                                                                       #
#########################################################################################################

if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(FactoMineR)) install.packages("FactoMineR")
library(FactoMineR)
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)

# Create parser, which accepts arguments from terminal and pass to the script
parser <- arg_parser("This script runs Correspondence analysis on the binarized set of the input deletion data and, if selected, on the sets of just enzymatic genes and/or genes included in iEK1011 Genome Scale Metabolic Model.")

parser <- add_argument(parser = parser,
                       arg = c("input", 
                               "--output"),
                       help = c("Input deletereous SNP matrix",
                                "Output directory where results of the analysis will be stored"),
                       flag = c(F, F),
                       default = list("input" = "/storage/PGO/results/mtb/snpsPipe/deletSNPsMat.csv",
                                      "--output" = "/storage/PGO/results/mtb/snpsPipe/CA/"))
                                      
parser <- add_argument(parser = parser, 
                       arg = "--multSampSize",
                       help = "Multiply by the sampled size of either animal or human lineages if specified.",
                       flag = T)

parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################

inFile <- parsed$input
outDir <- parsed$output


# Load and parse the data
#########################################################################################################
print(sprintf("Loading %s potentially deletereous SNP matrix...", basename(inFile)))
deletSNPs <- read.csv(inFile, row.names = 1, check.names = F)

if(parsed$multSampSize == T){
        sampSizes <- c(400, 178)
        names(sampSizes) <- c("A", "L")
        for(i in 1:nrow(deletSNPs)){
                lin <- rownames(deletSNPs)[[i]]
                linGroup <- gsub("[[:digit:]]", "", lin)
                deletSNPs[i, ] <- deletSNPs[i, ] * sampSizes[names(sampSizes) == linGroup]
        }
}

genesWSNPs <- unique(gsub("\\..*", "", colnames(deletSNPs)))

genesWSNPs <- genesWSNPs[genesWSNPs != "NA"]

# Build a matrix accounting for the sums of the proportions within each lineage of all the SNPs that 
# affect a single gene
deletSNPs_sum <- data.frame(matrix(ncol = nrow(deletSNPs),
                                   nrow = 0,
                                   dimnames = list(NULL,
                                                   rownames(deletSNPs))))
                                                   
print(sprintf("Building SNPs sums matrix..."))
for(i in seq_along(genesWSNPs)){
        gene <- genesWSNPs[i]
        subMat <- deletSNPs[, grep(gene, colnames(deletSNPs))]
        if(is.null(dim(subMat))){
                sums <- subMat
        }else{
                sums <- rowSums(subMat)
        }
        toBind <- matrix(sums,
                         ncol = nrow(deletSNPs),
                         dimnames = list(gene,
                                         rownames(deletSNPs)))
        
        deletSNPs_sum <- rbind.data.frame(deletSNPs_sum, toBind)
}

# Write the summed deletereous SNP DF to a CSV file

outName <- gsub(".csv", "_sums.csv", basename(inFile))

write.csv(deletSNPs_sum, file = paste0(outDir, outName))

print(sprintf("%s saved at %s.", outName, outDir))

# Do Correspondence Analysis and plot the results. 
#########################################################################################################
deletSNPsCA <- CA(deletSNPs_sum, graph = F)

pdf(paste0(outDir, "deletSNPsCA.pdf"))
fviz_ca_biplot(deletSNPsCA, label = "col", repel = T)
dev.off()


print(sprintf("deletSNPsCA.pdf saved at %s.", outDir))










