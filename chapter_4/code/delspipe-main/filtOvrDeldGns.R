#########################################################################################################
#                                                                                                       #
# Filter the genes that were found to be significantly differentially deleted between each lineage and  #
# reference lineage to keep just the ones that are deleted in more than the a given threshold           #
# proportion of the isolates within the lineage and intact in more than that same proportion of the     #
# isolates in reference, and that are present in iEK1011 2.0 GSMM. It also saves the RData object with  #
# the actual proportions within each lineage with the genes deleted for further exploration.            #
#                                                                                                       #
#########################################################################################################

if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(rtracklayer)) install.packages("rtracklayer", repos='http://cran.us.r-project.org')
library(rtracklayer)

parser <- arg_parser("This script takes as input the csv files created with the del proportion files (with statistically significant deletions between L4 and each lineage) and filters it to keep just the genes that are deleted in most of one of the two groups.")

parser <- add_argument(parser = parser,
                       arg = c("input", 
                               "--output", 
                               "--thrshld",
                               "--iEKGenes",
                               "--annotFile",
                               "--ref"),
                       help = c("Input directory where the files are stored",
                                "Output directory where output files will be stored", 
                                "Threshold to use in the filtering",
                                "File with genes included in iEK1011 2.0 GSMM",
                                "Annotation file",
                                "Reference lineage"),
                       default = list("input" = "/storage/PGO/results/mtb/delsPipe/signDelProp/",
                                      "--output" = "/storage/PGO/results/mtb/delsPipe/toRem/",
                                      "--thrshld" = 85,
                                      "--iEKGenes" = "/storage/PGO/data/mtb/GSMMs/iEK1011_2.0/iEK1011_2.0_genes.txt",
                                      "--annotFile" = "/storage/PGO/data/mtb/annotations/genes.gff",
                                      "--ref" = "L4"),
                       flag = c(F, F, F, F, F, F))

parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################
dataDir <- parsed$input
outDir <- parsed$output
signDelPrFiles <- list.files(dataDir)
signDelPrFiles <- signDelPrFiles[grep(".csv", signDelPrFiles)]
iEKFile <- parsed$iEKGenes
annotFile <- parsed$annotFile

# Load iEK1011 gene file and annotation.
#########################################################################################################
iEKGenes <- as.character(read.table(iEKFile, header = F, sep = "\t")$V1)
annot <- as.data.frame(readGFF(annotFile))

# Do the filtering and save files.
#########################################################################################################
lins <- substr(signDelPrFiles, 1, 2)
ref <- parsed$ref

thrshld <- parsed$thrshld

toRemDFs <- list()
genes2RemLst <- list()
for(i in seq_along(signDelPrFiles)){
        lin <- lins[i]
        path <- paste(dataDir, signDelPrFiles[i], sep = "")
        delPr <- read.csv(path)
        colnames(delPr)[1] <- "lin"
        rownames(delPr) <- delPr$lin
        toRem <- delPr[delPr$lin %in% c(lin, ref), 
                       c(T, delPr[delPr$lin == lin, 
                                  2:ncol(delPr)] > thrshld & delPr[delPr$lin == ref, 
                                                              2:ncol(delPr)] < thrshld)]
        toRemRef <- delPr[delPr$lin %in% c(lin, ref), 
                          c(T, delPr[delPr$lin == lin, 
                                     2:ncol(delPr)] < thrshld & delPr[delPr$lin == ref, 
                                                                      2:ncol(delPr)] > thrshld)]
        toRemLst <- list(toRem, toRemRef)
        names(toRemLst) <- c(paste("from", 
                                   lin, 
                                   sep = "_"), 
                             sprintf("from_%s", ref))
        toRemDFs[[i]] <- toRemLst
        genes2Rem <- sort(unique(gsub("\\_.*", "", colnames(toRem)[2:ncol(toRem)])))
        genes2Rem <- genes2Rem[genes2Rem %in% iEKGenes]
        genes2RemLst[[lin]] <- genes2Rem
        outName <- paste(outDir, lin, "toRem", ".txt", sep = "")
        write.table(genes2Rem, 
                    file = outName, 
                    sep = "\t", 
                    col.names = F,
                    row.names = F,
                    quote = F)
        print(sprintf("%s saved at %s", basename(outName), outDir))
}

save(toRemDFs, 
     file = paste(outDir, "toRemDFs.RData", sep = ""))
print(sprintf("toRemDFs.RData saved at %s.", outDir))

# Generate a data frame with the information of the genes that are gonna be removed from iEK1011 based on
# deletion data for generating lineage specific GSMMs
maxLen <- max(sapply(genes2RemLst, length))

toRem_dels_df <- data.frame(matrix(nrow = maxLen,
                           ncol = 0))
                           
for(l in lins){
         locus2Rem <- genes2RemLst[[l]]
         genes2Rem <- annot$gene[match(locus2Rem, annot$locus_tag) + 1]
         prdct2Rem <- annot$product[match(locus2Rem, annot$locus_tag) + 1]
         len2Rem <- length(locus2Rem)
         naNum <- maxLen - len2Rem
         locus2Rem <- c(locus2Rem, rep(NA, naNum))
         genes2Rem <- c(genes2Rem, rep(NA, naNum))
         prdct2Rem <- c(prdct2Rem, rep(NA, naNum))

         toRem_dels_df[[sprintf("%s_locus", l)]] <- locus2Rem
         toRem_dels_df[[sprintf("%s_gene", l)]] <- genes2Rem
         toRem_dels_df[[sprintf("%s_product", l)]] <- prdct2Rem
}
                           
write.csv(toRem_dels_df, paste0(outDir, "dels_toRemDF.csv"))
print(sprintf("dels_toRemDF.csv saved at %s.", outDir))