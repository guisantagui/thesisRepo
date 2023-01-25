##########################################################################################################
#                                                                                                        #
# Parses stop gain SNP files and PROVEAN significant SNPs in a single matrix where columns are unique    #
# SNPs and rows are lineages. Each cell corresponds to the proportion of SNPs within the lineage with    #
# the SNP.                                                                                               #
#                                                                                                        #
##########################################################################################################

if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)

# Create parser, which accepts arguments from terminal and pass to the script
parser <- arg_parser("Parses stop gain SNP files and PROVEAN significant SNPs in a single matrix where columns are unique SNPs and rows are lineages.")

parser <- add_argument(parser = parser,
                       arg = "--output",
                       help = "Path to output directory of the whole pipeline",
                       flag = F,
                       default = "/storage/PGO/results/mtb/snpsPipe/")

parsed <- parse_args(parser)

#
# Directory stuff
##########################################################################################################
outDir <- parsed$output
provDir <- paste0(outDir, "provean/")
stopGDir <- paste0(outDir, "stopGained/")

#
# Load and format data
##########################################################################################################
snpsProv <- read.csv(paste0(provDir, "allSNPsNonSyn_prov.csv"), stringsAsFactors = F, row.names = 1)
snps <- read.csv(paste0(outDir, "allSNPsNonSyn.csv"), stringsAsFactors = F) # This line should not be necessary when running proveanRun.py again, as proportions will be kept


snpsProv <- snpsProv[, c("POS", "REF", "ALT", "LIN", "ANN", "LOCUS", "AA_one", "PROVEAN")]

snpsProv$PROVEAN <- sapply(snpsProv$PROVEAN, function(x) as.numeric(gsub("\nName:", "", strsplit(x, " ")[[1]][5])))

snpsProv$propSNP <- snps$propSNP[!is.na(snps$LOCUS)]

snpsProv <- snpsProv[!is.na(snpsProv$PROVEAN), ]


snpsProv <- snpsProv[snpsProv$PROVEAN < -2.5, ]

deletSNPs <- snpsProv[, c("LOCUS", "AA_one", "propSNP", "LIN")]

stopGFiles <- list.files(stopGDir)

stopGFiles <- stopGFiles[gsub("\\_.*", "", stopGFiles) != "NA"]


for(i in seq_along(stopGFiles)){
        stopGFile <- stopGFiles[i]
        lin <- gsub("\\_.*", "", stopGFile)
        tab <- read.table(paste0(stopGDir, stopGFile), sep = "\t", header = T)
        print(sprintf("Reading %s...", stopGFile))
        tab$LIN <- rep(lin, nrow(tab))
        print(sprintf("Adding stop gain SNPs of lineage %s to deletereous SNP dataframe...", lin))
        deletSNPs <- rbind.data.frame(deletSNPs, tab)
}

deletSNPs$snpCode <- paste(deletSNPs$LOCUS, deletSNPs$AA_one, sep = ".")

write.csv(deletSNPs, file = paste0(outDir, "deletSNPs.csv"))

print(sprintf("deletSNPs.csv saved at %s.", outDir))

#
# Create the matrix
##########################################################################################################

uniqueSNPs <- sort(unique(deletSNPs$snpCode))

deletSNPsMat <- matrix(rep(0, length(unique(deletSNPs$LIN)) * length(uniqueSNPs)),
                       nrow = length(unique(deletSNPs$LIN)),
                       ncol = length(uniqueSNPs),
                       dimnames = list(sort(unique(deletSNPs$LIN)),
                                       uniqueSNPs))

print("Generating deletereous SNPs matrix...")
for(i in seq_along(uniqueSNPs)){
        snp <- uniqueSNPs[i]
        subMat <- deletSNPs[deletSNPs$snpCode == snp, ]
        inLins <- subMat$LIN
        deletSNPsMat[match(inLins, row.names(deletSNPsMat)), i] <- subMat$propSNP#[match(inLins, row.names(deletSNPsMat))]
}

# Write the matrix to a CSV file

write.csv(deletSNPsMat, paste0(outDir, "deletSNPsMat.csv"))

print(sprintf("Deletereous SNPs matrix written in deletSNPsMat.csv and saved at %s.", outDir))