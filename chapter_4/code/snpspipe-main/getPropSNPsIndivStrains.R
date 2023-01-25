##########################################################################################################
#                                                                                                        #
# Parses the deleterious SNP file of the sampled individual strains into a CSV file containing the       # 
# proportion within the sampled members of each lineage that has each SNP                                #
#                                                                                                        #
##########################################################################################################
if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(rtracklayer)) install.packages("rtracklayer", repos='http://cran.us.r-project.org')
library(rtracklayer)
# Create parser, which accepts arguments from terminal and pass to the script
parser <- arg_parser("Parses the deleterious SNP file of the sampled individual strains into a CSV file containing the proportion within the sampled members of each lineage that has each SNP")

parser <- add_argument(parser = parser,
                       arg = c("--indivSNPDir", 
                               "--output"),
                       help = c("Path to individual strain SNP directory",
                                "Output directory"),
                       flag = c(F, F),
                       default = list("--indivSNPDir" = "/storage/PGO/results/mtb/snpsPipe_bestQual/indivSNPs/",
                                      "--output" = "/storage/PGO/results/mtb/snpsPipe_bestQual/"))
                                      
parsed <- parse_args(parser)
                                      
#
# Directory stuff
##########################################################################################################

indivSNPDir <- parsed$indivSNPDir
output <- parsed$output


#
# Read dataframes and obtain proportion
##########################################################################################################
print("Generating proportion SNPs within lineage dataframe based on the sampled individual strains...")

lins <- c("A1",
          "A2",
          "A3",
          "A4",
          "L1",
          "L2",
          "L3",
          "L4",
          "L5",
          "L6",
          "L7",
          "L8",
          "L9")

indivSNPFiles <- list.files(indivSNPDir)

allLinProps <- data.frame()
for(lin in lins){
        linFile <- indivSNPFiles[grep(lin, indivSNPFiles)]
        print(sprintf("Parsing %s...", linFile))
        linSNPs <- read.csv(paste0(indivSNPDir, linFile),
                            #row.names = 1,
                            stringsAsFactors = F,
                            check.names = F)
        colnames(linSNPs)[1] <- "strain"
        linProps <- apply(linSNPs[colnames(linSNPs) != "strain"], 2, function(x) sum(x)/length(x))
        linPropDF <- data.frame(matrix(linProps,
                                       nrow = 1,
                                       ncol = length(linProps),
                                       dimnames = list(lin, 
                                                       names(linProps))),
                                check.names = F)
        if(nrow(allLinProps) == 0){
                allLinProps <- linPropDF
        }else{
                allLinProps <- rbind.data.frame(allLinProps, linPropDF)
        }
}

write.csv(allLinProps, paste0(output, "deletSNPsMat_indiv.csv"))

print(sprintf("deletSNPsMat_indiv.csv saved at %s.", output))