#########################################################################################################
#                                                                                                       #
# Generate binarized versions of input deletion file with the option of filtering the enzymatic genes   #
# and exporting a binarized set of them.                                                                #
#                                                                                                       #
#########################################################################################################
if(!require(dplyr)) install.packages("dplyr", repos='http://cran.us.r-project.org')
library(dplyr)
if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(stringr)) install.packages("stringr", repos='http://cran.us.r-project.org')
library(stringr)

# Create parser, which accepts arguments from terminal and pass to the script
parser <- arg_parser("This script generates binarized versions of input deletion file, with the option of filtering the enzymatic genes and exporting a binarized set of them.")

parser <- add_argument(parser = parser,
                       arg = c("input", 
                               "--output", 
                               "--threshold"),
                       help = c("Input deletion percentage file",
                                "Output directory where binarized output deletion files will be stored", 
                                "Threshold percentage of deleted ORF used to determine the gene as nonfunctional"),
                       flag = c(F, F, F),
                       default = list("input" = "/storage/PGO/results/mtb/deletions_allGenes/allDels.csv",
                                      "--output" = "/storage/PGO/results/mtb/deletions_analysis/",
                                      "--threshold" = 15))
                                      
parser <- add_argument(parser = parser, 
                       arg = "--enzFilt",
                       help = "Filters enzymatic genes before binarization and outputs enzymatic bin file too",
                       flag = T)
                                      
parsed <- parse_args(parser)

# Load functions
source("functions.R")

# Directory stuff
#########################################################################################################
input <- parsed$input
outDir <- parsed$output

outName <- gsub(".csv", "", basename(input))



# Load data and get the dataset of enzymatic genes
#########################################################################################################
dels <- read.csv(input)
colnames(dels)[1] <- "gene"
print(sprintf("Total number of genes: %s.", as.character(nrow(dels))))

# Set a threshold of proportion of deletion gene and binarize deletion percentage based in this threshold
######################################################################################################### 
thrshld <- parsed$threshold

delsBin <- binarizeDels(DM = dels, thrshld = thrshld)
write.csv(delsBin, sprintf("%s%sBin.csv", outDir, outName))

if(parsed$enzFilt == T){
        enzDels <- filtEnzDels(dels)
        print(sprintf("Number of enzymatic genes: %s.", as.character(nrow(enzDels))))
        print(sprintf("Percentage of enzymatic genes: %s%%.", 
                      as.character(round(nrow(enzDels)/nrow(dels)*100, digits = 2))))
        enzDelsBin <- binarizeDels(DM = enzDels, thrshld = thrshld)
        enzOutName <- paste0("enz", str_match(outName, "Dels"))
        write.csv(enzDelsBin, sprintf("%s%sBin.csv", outDir, enzOutName))
        print(sprintf("%sBin.csv and %sBin.csv saved at %s.", outName, enzOutName, outDir))
}else{
        print(sprintf("%sBin.csv saved at %s.", outName, outDir))
}