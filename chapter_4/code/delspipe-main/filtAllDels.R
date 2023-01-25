#########################################################################################################
#                                                                                                       #
# Filters allDels.csv file to keep just genomes included in labkey file (high quality and with          #
# accession number).                                                                                    #
#                                                                                                       #
#########################################################################################################

if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)

# Create parser, which accepts arguments from terminal and pass to the script
parser <- arg_parser("Filters allDels.csv file to keep just genomes included in labkey file (high quality and with accession number).")

parser <- add_argument(parser = parser,
                       arg = c("--allDels", 
                               "--labkeyFile"),
                       help = c("Input deletion file",
                                "Path to labkey database file"),
                       flag = c(F, F),
                       default = list("--allDels" = "/storage/PGO/results/mtb/deletions_allGenes/allDels.csv",
                                      "--labkeyFile" = "/storage/PGO/data/mtb/genomesInfo/labkey_qual.csv"))
                                      
parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################

allDelsFile <- parsed$allDels
labkeyFile <- parsed$labkeyFile

# Create output directory to allow us to keep the original name.
outDir <- paste0(dirname(allDelsFile), "/allDelsFilt/")

if(!dir.exists(outDir)){
        dir.create(outDir)
}

# Load data 
#########################################################################################################

print("Loading allDels.csv...")
allDels <- read.csv(allDelsFile, stringsAsFactors = F)
labkey <- read.csv(labkeyFile, row.names = 1, stringsAsFactors = F)

# Filter allDels.csv
#########################################################################################################

allDelsFilt <- allDels[, gsub("\\_.*", "", colnames(allDels)) %in% labkey$G.NUMBER]

lins <- sort(unique(gsub(".*_", "", colnames(allDelsFilt))))

linListsDir <- paste0(outDir, "linLists/")
if(!dir.exists(linListsDir)){
        dir.create(linListsDir)
        print(sprintf("%s created.", linListsDir))
}else{
        print(sprintf("%s already exists.", linListsDir))
}
for(lin in lins){
        linGNums <- gsub("\\_.*", 
                         "", 
                         colnames(allDelsFilt)[grep(lin, 
                                               colnames(allDelsFilt))])
        write.table(linGNums, 
                    file = sprintf("%s%snr.txt", linListsDir, lin),
                    quote = F,
                    sep = "\t",
                    row.names = F,
                    col.names = F)
        print(sprintf("%snr.txt saved at %s.", lin, linListsDir))
}

# Add the columns of gene information

allDelsInfo <- allDels[sapply(allDels, class) != "numeric"]
allDelsFilt <- cbind.data.frame(allDelsInfo, allDelsFilt)

write.csv(allDelsFilt, file = paste0(outDir, "allDels.csv"), row.names = F)

print(sprintf("allDels.csv saved at %s.", outDir))