#########################################################################################################
#                                                                                                       #
# Sum the individual Gnumber SNPs affecting a single gene in a single isolate                           #
#                                                                                                       #
#########################################################################################################

if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)

# Create parsed, which would accept arguments from terminal and pass to the script
parser <- arg_parser("Sum the individual Gnumber SNPs affecting a single gene in a single isolate")

parser <- add_argument(parser = parser,
                       arg = "--indSNPDir",
                       help = "Directory where individual G-num SNP files are",
                       flag = F,
                       default = "/storage/PGO/results/mtb/snpsPipe/indivSNPs/")
                       
parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################
outDir <- parsed$indSNPDir
indSNPFile <- paste0(outDir, "deletSNPsMat_byStrain.csv")

# LoadData
#########################################################################################################

indSNPs <- read.csv(indSNPFile, stringsAsFactors = F)
colnames(indSNPs)[1] <- "gNum"

print(head(colnames(indSNPs)))

sumGeneSNPs <- function(snpMat){
        snpMatSNPs <- colnames(snpMat)
        snpMatSNPs <- snpMatSNPs[snpMatSNPs != "gNum"]
        unGns <- unique(gsub("\\..*", "", snpMatSNPs))
        unGns <- sort(unGns)
        summd <- data.frame()
        for(i in seq_along(unGns)){
                gene <- unGns[i]
                print(sprintf("Summing SNPs in %s...", gene))
                subMat <- snpMat[, grep(gene, colnames(snpMat))]
                if(is.null(dim(subMat))){
                        toCBind <- subMat
                }else{
                        toCBind <- rowSums(subMat)
                }
                if(ncol(summd) == 0){
                        summd <- data.frame(matrix(toCBind, 
                                                   ncol = 1, 
                                                   nrow = nrow(snpMat)))
                }else{
                        summd <- cbind.data.frame(summd, toCBind)
                }
        }
        print("Number of g numbers")
        print(length(make.unique(snpMat$gNum)))
        print("Nrow output")
        print(nrow(summd))
        colnames(summd) <- unGns
        print(head(snpMat$gNum))
        rownames(summd) <- make.unique(snpMat$gNum)
        summd$gNum <- snpMat$gNum
        return(summd)
}

indSNPs_sum <- sumGeneSNPs(indSNPs)

write.csv(indSNPs_sum, file = paste0(outDir, "deletSNPsMat_byStrain_sums.csv"))
print(sprintf("deletSNPsMat_byStrain_sums.csv saved at %s", outDir))