#########################################################################################################
#                                                                                                       #
# Plot variable importance obtained with Random Forest adding colors corresponding to RDs. Two          #
# different plots are generated, one with the RD denomination in Behr et al 1999 "Comparative genomics  # 
# of BCG vaccines by whole-genome DNA microarray" and the other with the one in Brosch et. al 2002      #
# "A new evolutionary scenario for the Mycobacterium tuberculosis complex"                              #
#                                                                                                       #
#########################################################################################################

if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(rtracklayer)) install.packages("rtracklayer", repos='http://cran.us.r-project.org')
library(rtracklayer)

# Create parser, which accepts arguments from terminal and pass to the script
parser <- arg_parser("This script plots top N important variables, one colored according with Behr et al. 1999 RD denomination and the other one with Brosch et al. 2002.")

parser <- add_argument(parser = parser,
                       arg = c("--allDels",
                               "--output",
                               "--annot",
                               "--top"),
                       help = c("Path to allDels file",
                                "Output directory where plots will be stored",
                                "Path to annotation file",
                                "Number of top important variables to plot"),
                       flag = c(F, F, F, F),
                       default = list("--allDels" = "/storage/PGO/results/mtb/deletions_allGenes/allDelsFilt/allDels.csv", 
                                      "--output" = "/storage/PGO/results/mtb/delsPipe_bestQual/RF/",
                                      "--annot" = "/storage/PGO/data/mtb/annotations/genes.gff",
                                      "--top" = "50"))
                                      
parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################
outDir <- parsed$output

# Select the important variables file:
inputFile <- list.files(outDir)
inputFile <- inputFile[grep(".csv", inputFile, fixed = T)]
inputFile <- inputFile[grep("varImp", inputFile, fixed = T)]

inputPath <- paste0(outDir, inputFile)
annoPath <- parsed$annot

print(inputPath)

# Load and parse the data
#########################################################################################################

# Variable importance data...
print(sprintf("Loading %s from %s...", inputFile, outDir))
vipRF <- read.csv(inputPath, row.names = 1)

# All deletion data (for obtaining EC numbers)
print(sprintf("Loading %s from %s...", 
              basename(parsed$allDels), 
              dirname(parsed$allDels)))

allDels <- read.csv(parsed$allDels)

# And annotation data...
annot <- readGFF(annoPath)

# Obtain dataframe of top N variables. 
top <- parsed$top

vipRFTopDF <- data.frame(locus_tag = vipRF$gene[1:top],
                         gene = annot$gene[match(vipRF$gene, annot$locus_tag)+1][1:top],
                         product = annot$product[match(vipRF$gene, annot$locus_tag)+1][1:top],
                         ECnum = allDels$EC_number[match(vipRF$gene[1:top], allDels$gene)],
                         importance = vipRF$Overall[1:top])

# Create dataframes of RDs (Brosch and Behr).
#########################################################################################################
rdTab <- data.frame(RD = paste("RD", as.character(1:16), sep = ""),
                    genes = c(paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3871"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3871")+8)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1978"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1978")+12)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1573"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1573")+13)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv0221"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv0221")+2)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3117"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3117")+4)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1506c"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1506c")+12)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2346c"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2346c")+7)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv0309"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv0309")+3)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3617"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3617")+6)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1255c"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1255c")+2)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3425"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3425")+4)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2072c"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2072c")+3)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2645"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2645")+13)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1766"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1766")+7)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1963c"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1963c")+14)], collapse = ","),
                              paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3400"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3400")+5)], collapse = ",")))


rdBrosch <- data.frame(RD = paste("RD", as.character(1:14), sep = ""),
                       genes = c(paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3871"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3871")+8)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1978"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1978")+12)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1573"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1573")+13)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1505c"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1506c")+12)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2346c"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2346c")+7)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3425"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3425")+3)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1964"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1964")+13)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3617"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3617")+6)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2072c"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2072c")+3)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv0221"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv0221")+2)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2645"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv2645")+12)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3118"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv3118")+3)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1255c"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1255c")+2)], collapse = ","),
                                 paste(annot$locus_tag[!is.na(annot$locus_tag)][which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1765c"):(which(annot$locus_tag[!is.na(annot$locus_tag)] == "Rv1765c")+8)], collapse = ",")))

# Define functions for adding the RD of each gene.
#########################################################################################################
getTopRD2PlotDF <- function(topVarImp = vipRFTopDF, RD = "behr"){
        vipRFTopRD_df <- topVarImp[nrow(topVarImp):1, ]
        rdVec <- c()
        for(i in seq_along(vipRFTopRD_df$locus_tag)){
                gene <- vipRFTopRD_df$locus_tag[i]
                trueVec <- c()
                if(RD == "behr"){
                        RDs <- rdTab
                }else if(RD == "brosch"){
                        RDs <- rdBrosch
                }
                for(j in seq_along(RDs$RD)){
                        rdGenes <- strsplit(as.character(RDs$genes[j]), ",")[[1]]
                        trueVec <- c(trueVec, gene %in% rdGenes)
                }
                if(sum(trueVec) == 0){
                        RD <- "No RD mapped."
                }else{
                        RD <- as.character(RDs$RD[trueVec])
                }
                rdVec <- c(rdVec, RD)
        }
        vipRFTopRD_df$RD <- rdVec
        return(vipRFTopRD_df)
}

vipRFTopRD_behr <- getTopRD2PlotDF(vipRFTopDF, RD = "behr")
vipRFTopRD_brosch <- getTopRD2PlotDF(vipRFTopDF, RD = "brosch")

# Prepare a dataframe with al the information for writing to CSV file.
outDF <- vipRFTopRD_brosch
colnames(outDF) <- gsub("RD", "RD_brosch", colnames(outDF))
outDF$RD_behr <- vipRFTopRD_behr$RD
outDF <- outDF[nrow(outDF):1, ]

write.csv(outDF, file = sprintf("%stop%sDels_RDs.csv", outDir, top))

print(sprintf("top%sDels_RDs.csv saved at %s.", top, outDir))

# Do the plots. 

behrPlotPath <- sprintf("%s%s_RDBehr.pdf", 
                        outDir, 
                        gsub(".csv", "", inputFile))
ggplot(data = vipRFTopRD_behr, mapping = aes(x = importance, 
                                             y = factor(as.character(locus_tag), 
                                                        levels = as.character(locus_tag)),
                                             fill = RD)) +
        geom_bar(stat = "identity") + 
        labs(x = "Variable importance", y = "Variant") + 
        theme_minimal()
ggsave(behrPlotPath, 
       width = 4, 
       height = 10)
       
print(sprintf("%s saved at %s.", basename(behrPlotPath), dirname(behrPlotPath)))

broschPlotPath <- sprintf("%s%s_RDBrosch.pdf", 
                          outDir, 
                          gsub(".csv", "", inputFile))
ggplot(data = vipRFTopRD_brosch, mapping = aes(x = importance, 
                                               y = factor(as.character(locus_tag), 
                                                          levels = as.character(locus_tag)),
                                               fill = RD)) +
        geom_bar(stat = "identity") + 
        labs(x = "Variable importance", y = "Variant") + 
        theme_minimal()
ggsave(broschPlotPath, 
       width = 4, 
       height = 10)
       
print(sprintf("%s saved at %s.", basename(broschPlotPath), dirname(broschPlotPath)))
