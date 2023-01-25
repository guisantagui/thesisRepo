#########################################################################################################
#                                                                                                       #
# Filters deletereous SNPs dataframe to keep just the ones that are present in more than a given        #
# proportion within each lineage.                                                                       #
#                                                                                                       #
#########################################################################################################

if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(rtracklayer)) install.packages("rtracklayer", repos='http://cran.us.r-project.org')
library(rtracklayer)

# Create parsed, which accepts arguments from terminal and pass to the script
parser <- arg_parser("This script plots results from Random Forest")

parser <- add_argument(parser = parser,
                       arg = c("--snpFile",
                               "--thrshld",
                               "--ref",
                               "--annoDir",
                               "--iEKGenes",
                               "--output"),
                       help = c("Directory where deletereous SNP file is located.",
                                "Threshold for considering SNP significatively present within lineage.",
                                "Reference lineage",
                                "Path to annotation directory",
                                "File with a list of genes included in iEK1011 2.0 GSMM", 
                                "Output directory where output files will be stored"),
                       flag = c(F, F, F, F, F, F),
                       default = list("--snpFile" = "/storage/PGO/results/mtb/snpsPipe/deletSNPsMat.csv",
                                      "--thrshld" = 0.85,
                                      "--ref" = "L4",
                                      "--annoDir" = "/storage/PGO/data/mtb/annotations/",
                                      "--iEKGenes" = "/storage/PGO/data/mtb/GSMMs/iEK1011_2.0/iEK1011_2.0_genes.txt",
                                      "--output" = "/storage/PGO/results/mtb/snpsPipe/toRem/"))
                       
parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################

snpFile <- parsed$snpFile
thrshld <- parsed$thrshld
annoFile <- paste0(parsed$annoDir, "genes.gff")
iEKFile <- parsed$iEKGenes
outDir <- parsed$output
outStopGain <- paste0(outDir, "stopGain/")
outMissProv <- paste0(outDir, "missProv/")

if(!dir.exists(outStopGain)){
        dir.create(outStopGain)
}

if(!dir.exists(outMissProv)){
        dir.create(outMissProv)
}

# Load data
#########################################################################################################
ref <- parsed$ref

# Load iEK1011 2.0 genes
iEKGenes <- as.character(read.table(iEKFile, header = F, sep = "\t")$V1)

# Load annotation data
annot <- as.data.frame(readGFF(annoFile))

print(sprintf("Loading deletereous SNPs dataframe %s from %s...", basename(snpFile), dirname(snpFile)))

SNPs <- read.csv(snpFile, 
                 row.names = 1, 
                 check.names = F)

# Filter SNPs dataframe to keep just genes with potentially deletereous SNPs present at more
# than a given proportion. Separate them according to if they are Stop Gained or missense PROVEAN
# significative, filter to keep just the ones included in iEK1011 2.0 model and write them in separate 
# TXT files. Write also a table with the genes and product (from annotation file) with potentially 
# deletereous SNPs in each lineage. 
#########################################################################################################

lins <- rownames(SNPs)

lins <- lins[lins != ref]

refMat <- SNPs[ref, ]
allSNPsVec <- colnames(SNPs)

toRem_SG <- list()
toRem_MS <- list()
prdcts_lst_SG <- list()
prdcts_lst_MS <- list()
genes_lst_SG <- list()
genes_lst_MS <- list()
for(i in seq_along(lins)){
        lin <- lins[i]
        subMat <- SNPs[lin, ]
        SNPs2Rem <- allSNPsVec[subMat >= thrshld & refMat <= thrshld ]
        SNPs2Rem_SG <- SNPs2Rem[grep("*", SNPs2Rem, fixed = T)]
        SNPs2Rem_MS <- SNPs2Rem[-grep("*", SNPs2Rem, fixed = T)]
        SNPs2Rem_SG <- sort(unique(gsub("\\..*", "", SNPs2Rem_SG)))
        SNPs2Rem_SG <- SNPs2Rem_SG[SNPs2Rem_SG %in% iEKGenes]
        toRem_SG[[lin]] <- SNPs2Rem_SG
        SNPs2Rem_MS <- sort(unique(gsub("\\..*", "", SNPs2Rem_MS)))
        SNPs2Rem_MS <- SNPs2Rem_MS[SNPs2Rem_MS %in% iEKGenes]
        toRem_MS[[lin]] <- SNPs2Rem_MS
        prdcts_SG <- annot$product[match(SNPs2Rem_SG, annot$locus_tag) + 1]
        prdcts_lst_SG[[lin]] <- prdcts_SG
        prdcts_MS <- annot$product[match(SNPs2Rem_MS, annot$locus_tag) + 1]
        prdcts_lst_MS[[lin]] <- prdcts_MS
        genes_SG <- annot$gene[match(SNPs2Rem_SG, annot$locus_tag) + 1]
        genes_lst_SG[[lin]] <- genes_SG
        genes_MS <- annot$gene[match(SNPs2Rem_MS, annot$locus_tag) + 1]
        genes_lst_MS[[lin]] <- genes_MS
        write.table(SNPs2Rem_SG, 
                    file = sprintf("%s%s_stopGain.txt", 
                                   outStopGain, 
                                   lin),
                    quote = F,
                    sep = "\t",
                    row.names = F,
                    col.names = F)
        print(sprintf("%s_stopGain.txt saved at %s.", lin, outStopGain))
        write.table(SNPs2Rem_MS, 
                    file = sprintf("%s%s_missProv.txt", 
                                   outMissProv, 
                                   lin),
                    quote = F,
                    sep = "\t",
                    row.names = F,
                    col.names = F)
        print(sprintf("%s_stopGain.txt saved at %s.", lin, outMissProv))
}

# Build dataframes to export...
maxLen_SG <- max(sapply(toRem_SG, length))
maxLen_MS <- max(sapply(toRem_MS, length))

SG_df <- data.frame(matrix(nrow = maxLen_SG,
                           ncol = 0))
MS_df <- data.frame(matrix(nrow = maxLen_MS,
                           ncol = 0))

for(l in lins){
         len2Rem_SG <- length(toRem_SG[[l]])
         len2Rem_MS <- length(toRem_MS[[l]])
         naNum_SG <- maxLen_SG - len2Rem_SG
         naNum_MS <- maxLen_MS - len2Rem_MS
         toRem_SG[[l]] <- c(toRem_SG[[l]], rep(NA, naNum_SG))
         genes_lst_SG[[l]] <- c(genes_lst_SG[[l]], rep(NA, naNum_SG))
         prdcts_lst_SG[[l]] <- c(prdcts_lst_SG[[l]], rep(NA, naNum_SG))
         SG_df[[sprintf("%s_locus", l)]] <- toRem_SG[[l]]
         SG_df[[sprintf("%s_gene", l)]] <- genes_lst_SG[[l]]
         SG_df[[sprintf("%s_product", l)]] <- prdcts_lst_SG[[l]]
         
         toRem_MS[[l]] <- c(toRem_MS[[l]], rep(NA, naNum_MS))
         genes_lst_MS[[l]] <- c(genes_lst_MS[[l]], rep(NA, naNum_MS))
         prdcts_lst_MS[[l]] <- c(prdcts_lst_MS[[l]], rep(NA, naNum_MS))
         MS_df[[sprintf("%s_locus", l)]] <- toRem_MS[[l]]
         MS_df[[sprintf("%s_gene", l)]] <- genes_lst_MS[[l]]
         MS_df[[sprintf("%s_product", l)]] <- prdcts_lst_MS[[l]]
}

write.csv(SG_df, paste0(outStopGain, "stopGain_toRemDF.csv"))
print(sprintf("stopGain_toRemDF.csv saved at %s.", outStopGain))

write.csv(MS_df, paste0(outMissProv, "missProv_toRemDF.csv"))
print(sprintf("missProv_toRemDF.csv saved at %s.", outMissProv))

print(sprintf("SNP based toRem files writen at %s.", outDir))