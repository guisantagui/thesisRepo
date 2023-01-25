library(factoextra)
library(FactoMineR)
library(gplots)
library(caret)

dataDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/data/"
delDir <- paste0(dataDir, "deletions/")
snpDir <- paste0(dataDir, "snps/")
whatMod <- "delsAllSNPs"


outDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/results/plots/"
resDir <- gsub("plots/", "", outDir)

delFile <- paste0(delDir, "dels_toRemDF.csv")
stopGainSNPFile <- paste0(snpDir, "stopGain_toRemDF.csv")
missProvSNPFile <- paste0(snpDir, "missProv_toRemDF.csv")


remRxnsDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/results/simulations/delsAllSNPs/"

rxnDetsFile <- paste0(dataDir, "rxnDets/iEK1011_2.0_rxnDets.csv")
remRxns <- read.csv(paste0(remRxnsDir, sprintf("rxnsRem_tab_%s.csv", whatMod)), 
                    row.names = 1, 
                    stringsAsFactors = F)

remGenes <- read.csv(paste0(gsub("plots/", "", outDir),
                            sprintf("%s_mat_cur.csv", gsub("ls", "l_", whatMod))),
                     row.names = 1)

delets <- read.csv(delFile, row.names = 1, stringsAsFactors = F)
stopGa <- read.csv(stopGainSNPFile, row.names = 1, stringsAsFactors = F)
msProv <- read.csv(missProvSNPFile, row.names = 1, stringsAsFactors = F)


rxnDets <- read.csv(rxnDetsFile, 
                    row.names = 1, 
                    stringsAsFactors = F)

rxnDets$subsystem <- gsub("Transpor", "Transport", rxnDets$subsystem)
rxnDets$subsystem <- gsub("Arabinogalactan bioynthesis", 
                          "Arabinogalactan biosynthesis", 
                          rxnDets$subsystem)


remRxns$name <- rxnDets$name[match(rownames(remRxns), rxnDets$id)]

remRxns$subsystem <- rxnDets$subsystem[match(rownames(remRxns), rxnDets$id)]

remRxns$GPR <- rxnDets$GPR[match(rownames(remRxns), rxnDets$id)]

getRemGensRxn <- function(removed){
        linRemVec <- c()
        for(i in 1:nrow(remRxns)){
                #print(i)
                genes <- unlist(strsplit(remRxns$GPR[i], " and | or "))
                lins <- unique(gsub("\\_.*","", colnames(removed)))
                linRemStr <- ""
                for(j in seq_along(lins)){
                        lin <- lins[j]
                        rem <- removed[, sprintf("%s_locus", lin)]
                        rem <- rem[!is.na(rem)]
                        remRxnLin <- rem[rem %in% genes]
                        if(length(remRxnLin) > 0){
                                linRem <- sprintf("%s: %s", lin, paste(remRxnLin, collapse = ", "))
                                if(nchar(linRemStr) > 0){
                                        linRemStr <- paste(linRemStr, linRem, sep = ". ")
                                }else{
                                        linRemStr <- linRem
                                }
                        }
                }
                if(nchar(linRemStr) == 0){
                        linRemStr = NA
                }
                linRemVec <- c(linRemVec, linRemStr)
        }
        length(linRemVec)
        names(linRemVec) <- rownames(remRxns)
        return(linRemVec)
}

remRxns$delets <- getRemGensRxn(delets)
remRxns$stopGa <- getRemGensRxn(stopGa)
remRxns$msProv <- getRemGensRxn(msProv)

View(remRxns)

write.csv(remRxns, file = paste0(resDir, "remRxnsTab.csv"))

CAMat_remRxns <- data.frame(matrix(nrow = 0, ncol = 13,
                                   dimnames = list(NULL,
                                                   colnames(remRxns)[1:13])))


uniqueSubsysts <- unique(remRxns$subsystem)
for(i in seq_along(uniqueSubsysts)){
        subsyst <- uniqueSubsysts[i]
        subRemRxns <- remRxns[remRxns$subsystem == subsyst, 1:13]
        subRemRxns <- abs(subRemRxns -1)
        subSystSums <- colSums(subRemRxns)
        subSystSums <- data.frame(matrix(subSystSums, 
                                         nrow = 1, 
                                         ncol = ncol(subRemRxns),
                                         dimnames = list(subsyst,
                                                         colnames(subRemRxns))))
        CAMat_remRxns <- rbind.data.frame(CAMat_remRxns, subSystSums)
}

#rownames(CAMat_remRxns) <- uniqueSubsysts


CA_subsystsDels <- CA(CAMat_remRxns)

fviz_ca_biplot(CA_subsystsDels, repel = T)

pdf(file = sprintf("%sremRxnsSubsystHM_%s.pdf", outDir, whatMod))
heatmap.2(as.matrix(CAMat_remRxns), 
          hclustfun = function(x) hclust(x, method = "ward.D"),
          trace = "none", 
          scale = "none", 
          col = greenred(12), 
          breaks = 13, 
          density.info = "none", 
          margins = c(4, 18),
          cexCol = 1.5,
          colCol = c("#c4bf62",
                     "#87e0a0",
                     "#6db3d6",
                     "#f279ce",
                     "#ff30c1",
                     "#001aff",
                     "#8826b5",
                     "#ff0000",
                     "#871414",
                     "#24ad37",
                     "#fbff00",
                     "#ff9d00",
                     "#37ff30"))
dev.off()


pdf(file = sprintf("%sremRxnsSubsystHM_%s.pdf", outDir, whatMod), height = 50, width = 50)
heatmap.2(as.matrix(CAMat_remRxns), 
          hclustfun = function(x) hclust(x, method = "ward.D"),
          trace = "none", 
          scale = "none", 
          col = greenred(12), 
          breaks = 13, 
          density.info = "none", 
          margins = c(0, 0),
          cexCol = 1.5,
          colCol = c("#c4bf62",
                     "#87e0a0",
                     "#6db3d6",
                     "#f279ce",
                     "#ff30c1",
                     "#001aff",
                     "#8826b5",
                     "#ff0000",
                     "#871414",
                     "#24ad37",
                     "#fbff00",
                     "#ff9d00",
                     "#37ff30"))
dev.off()


subsystRemRxnsPCA <- prcomp(as.matrix(t(CAMat_remRxns)))

fviz_pca_ind()


remGenesCaret <- as.data.frame(t(remGenes))
remGenesCaret$linGroup <- gsub("[[:digit:]]", "", rownames(remGenesCaret))

remRxnsCaret <- as.data.frame(t(remRxns[, rownames(remGenesCaret)]))
remRxnsCaret$linGroup <- gsub("[[:digit:]]", "", rownames(remRxnsCaret))

set.seed(123)
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid",
                          allowParallel = T, 
                          sampling = "up")

trainFitRemGenes <- train(linGroup ~ . ,
                          data = remGenesCaret,
                          method = "rf",
                          metric = "Accuracy",
                          trControl = trControl)

trainFitRemRxns <- train(linGroup ~ . ,
                         data = remRxnsCaret,
                         method = "rf",
                         metric = "Accuracy",
                         trControl = trControl)

predict(trainFitRemGenes, remGenesCaret[, 1:(ncol(remGenesCaret)-1)])
predict(trainFitRemRxns, remRxnsCaret[, 1:(ncol(remRxnsCaret)-1)])

varImp(trainFitRemGenes)

varImp(trainFitRemRxns)
