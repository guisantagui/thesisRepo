#########################################################################################################
#                                                                                                       #
# Run Random Forest on the input data, training to do classifications of animal or human associated     #
# lineages. Plot confusion matrix and barplot of top X important variables, if indicated. Export        #
# dataframe of variable importance. Has the option of running RF with just the genes included in        #
# iEK1011 Genome Scale Metabolic Model. Adapted to work with SNPs data.                                 #
#                                                                                                       #
#########################################################################################################
if(!require(caret)) install.packages("caret", repos='http://cran.us.r-project.org')
library(caret)
if(!require(randomForest)) install.packages("randomForest", repos='http://cran.us.r-project.org')
library(randomForest)
if(!require(dplyr)) install.packages("dplyr", repos='http://cran.us.r-project.org')
library(dplyr)
if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(ggplot2)) install.packages("ggplot2", repos='http://cran.us.r-project.org')
library(ggplot2)

# Create parsed, which would accept arguments from terminal and pass to the script
parser <- arg_parser("This script does random forest separating animal and human lineages")

parser <- add_argument(parser = parser,
                       arg = c("input", 
                               "--output", 
                               "--weights",
                               "--sample",
                               "--plots",
                               "--iEKFilt"),
                       help = c("Input file",
                                "Output directory where output files will be stored", 
                                "Do class weights for compensate class imbalance",
                                "Select sample algorythm (up, down or cust)",
                                "Output plots of confusion matrix and variable importance",
                                "Filter input dataset to genes contained in iEK1011 model"),
                       flag = c(F, F, T, F, T, T))
                       
parser <- add_argument(parser = parser, 
                       arg = "--iEKGenes",
                       help = "Path to .txt file containing a list with the genes included in iEK1011 Genome Scale Metabolic Model",
                       default = "/storage/PGO/data/mtb/GSMMs/iEK1011_2.0/iEK1011_2.0_genes.txt")

parser <- add_argument(parser = parser, 
                       arg = "--lineages",
                       help = "Select what lineages are gonna be used in RF (default: all of them)",
                       default = c("L1,L2,L3,L4,L5,L6,L7,L8,L9,A1,A2,A3,A4"))
                       
parser <- add_argument(parser = parser, 
                       arg = "--top",
                       help = "Number of top genes to be exported in text file",
                       default = 20)
                       
parser <- add_argument(parser = parser, 
                       arg = "--metricRF",
                       help = "Quality metric to be used in RF model",
                       default = "Accuracy")

parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################
linPatt <- gsub(",", "|", parsed$lineages)
resuDir <- parsed$output
metricRF <- parsed$metricRF

# Load data
#########################################################################################################
print(sprintf("Loading %s from %s...", basename(parsed$input), dirname(parsed$input)))
delsBin <- read.csv(parsed$input, row.names = 1)
genes <- rownames(delsBin)

# Load iEK1011 list of genes
if(parsed$iEKFilt == T){
        iEK1011_genes <- as.character(read.table(parsed$iEKGenes, sep = "\t")$V1)
        iEK1011_genes[iEK1011_genes == "Rv2357"] <- "Rv2357c"
        iEK1011_genes[iEK1011_genes == "Rv3556"] <- "Rv3556c"
}

# Set output name according to the chosen settings. 
#########################################################################################################

# Name for the outputs based on flags and input file's name 
outName <- basename(parsed$input)

if(outName != "allDels.csv" & outName != "deletSNPsMat_sums.csv"){
        outName <- gsub("Bin.csv", "", outName)
}
if(outName == "deletSNPsMat_byStrain_sums.csv"){             
        outName <- gsub(".csv", "", outName)
        outName <- gsub("Mat", "", outName)
        outName <- gsub("_sums", "", outName)
}else{
        outName <- gsub(".csv", "", outName)
        outName <- gsub("_sums", "", outName)
        outName <- paste(outName, "NoBin", sep = "")
}

print(outName)

if(parsed$weights == T){
        outName <- paste(outName, "Wghtd", sep = "")
}

if(is.na(parsed$sample)){
        outName <- outName
}else if(parsed$sample == "down"){
        outName <- paste(outName, "Dwn", sep = "")
}else if(parsed$sample == "up"){
        outName <- paste(outName, "Up", sep = "")
}else if(parsed$sample == "cust"){
        outName <- paste(outName, "Cust", sep = "")
}
if(parsed$lineages != "L1,L2,L3,L4,L5,L6,L7,L8,L9,A1,A2,A3,A4"){
        outName <- paste(outName, gsub(",", "", parsed$lineages), sep = "_")
}
if(parsed$iEKFilt == T){
        outName <- paste(outName, "iEKFilt", sep = "_")
}

# Parse the data to have the format needed by RF function
#########################################################################################################

# Create datasets without gene information
# delsBinRdx <- delsBin[, 9:ncol(delsBin)]

print(outName)
if(outName == "deletSNPs_byStrain"){
        gNums <- delsBin$gNum
        print(head(gNums))
        print(length(gNums))
        print(dim(delsBin))
        colnames(delsBin)
        # delsBinRdx <- delsBin[, sapply(delsBin, class) == "numeric"]
        delsBinRdx <- delsBin[, colnames(delsBin) != "gNums"]
        rownames(delsBinRdx) <- make.unique(as.character(gNums))
}else{
        delsBinRdx <- delsBin[, sapply(delsBin, class) == "numeric"]
        rownames(delsBinRdx) <- genes
}


# If --iEKFilt flag is present filter the delsBin dataset previously of running RF to include only genes in iEK1011
if(parsed$iEKFilt == T){
        delsBinRdx <- delsBinRdx[rownames(delsBinRdx) %in% iEK1011_genes, ]
        print("in")
        print(iEK1011_genes[iEK1011_genes %in% rownames(delsBinRdx)])
        print("not in")
        print(iEK1011_genes[!iEK1011_genes %in% rownames(delsBinRdx)])
}

if(outName != "deletSNPs_byStrain"){
        delsBinRdx <- as.data.frame(t(delsBinRdx))
}


# Filter the dataframe to keep just selected lineages
delsBinRdx <- delsBinRdx[grep(linPatt, rownames(delsBinRdx)), ]


# Add column indicating if is animal or human lineage
linGroup <- factor(gsub(".", 
                        "", 
                        gsub("[[:digit:]]", 
                             "", 
                             gsub(".*_", 
                                  "", 
                                  rownames(delsBinRdx))), 
                        fixed = T))
delsBinRdx$linGroup <- linGroup


# Do sampling for compensate for class imbalance, if chosen.
#########################################################################################################

# If sample algorithm flag is present set in in trainControl

# custSamp: samples in a way that all lineages have proportional representation and
# animal and human associated lineages have the same members. Performs subsampling 
# without resampling in overrepresented lineages and upsampling with replacement  
# in underrepresented lineages. 
custSamp <- list(name = "Custom sampling for having all lineages equally represented", 
                 func = function(x, y){
                   dat <- if (is.data.frame(x)) x else as.data.frame(x)
                   dat$.y <- y
                   linMembs <- table(gsub(".*_", "", rownames(dat)))
                   linNames <- names(linMembs)
                   linGroupMembs <- sort(table(gsub("[[:digit:]]", 
                                                    "", 
                                                    linNames)), 
                                         decreasing = T)
                   medLinMembs <- median(linMembs)
                   toSampFromBigGroup <- medLinMembs
                   toSampFromSmaGroup <- medLinMembs*((medLinMembs*max(linGroupMembs))/(medLinMembs*min(linGroupMembs)))
                   sampSize <- c(toSampFromBigGroup, toSampFromSmaGroup)
                   names(sampSize) <- names(linGroupMembs)
                   sampled <- c()
                   for(i in 1:length(linNames)){
                           lin <- linNames[i]
                           size <- round(sampSize[match(gsub("[[:digit:]]", "", lin), names(sampSize))])
                           linSize <- linMembs[match(lin, linNames)]
                           sampsInLin <- rownames(dat)[grep(lin, rownames(dat))]
                           if(size <= linSize){
                                   linSamp <- sample(sampsInLin, size = size, replace = F)
                           }else if(size > linSize){
                                   linSamp <- sample(sampsInLin, size = size, replace = T)
                           }
                           sampled <- c(sampled, linSamp)
                   }
                   sampDF <- dat[match(sampled, rownames(dat)), !grepl(".y", colnames(dat), fixed = TRUE)]
                   classes <- gsub("[[:digit:]]", "", gsub(".*_", "", rownames(sampDF)))
                   classes <- as.factor(gsub(".", "", classes, fixed = T))
                   output <- list(x = sampDF, 
                                  y = classes)
                   return(output)
                   },
                 first = TRUE)


if(is.na(parsed$sample)){
        trControl <- trainControl(method = "cv",
                                  number = 10,
                                  search = "grid",
                                  allowParallel = T)
}else if(parsed$sample == "down"){
        trControl <- trainControl(method = "cv",
                                  number = 10,
                                  search = "grid",
                                  allowParallel = T, 
                                  sampling = "down")
}else if(parsed$sample == "up"){
        trControl <- trainControl(method = "cv",
                                  number = 10,
                                  search = "grid",
                                  allowParallel = T, 
                                  sampling = "up")
}else if(parsed$sample == "cust"){
        trControl <- trainControl(method = "cv",
                                  number = 10,
                                  search = "grid",
                                  allowParallel = T, 
                                  sampling = custSamp)
}

if(!is.na(parsed$sample)){
        print(paste("Subsampling method selected:", parsed$sample))
}

# Fit RF model and obtain the quality values. If chosen, plot confusion matrix. 
#########################################################################################################

# print(head(delsBinRdx))

# If weights flag is present, do weights vector and add it as argument to train function
print("Fitting Random Forest model to the data...")
if(parsed$weights == T){
        weightL <- (table(delsBinRdx$linGroup)[1]/sum(table(delsBinRdx$linGroup)))/(table(delsBinRdx$linGroup)[2]/sum(table(delsBinRdx$linGroup)))
        weightA <- (table(delsBinRdx$linGroup)[1]/sum(table(delsBinRdx$linGroup)))/(table(delsBinRdx$linGroup)[1]/sum(table(delsBinRdx$linGroup)))
        modelWeights <- c(weightA, weightL)
        names(modelWeights) <- c("A", "L")
        modelWeights <- modelWeights[match(allDelsBinRdx$linGroup, names(modelWeights))]
        #set.seed(123) # Already tried, same
        set.seed(999)
        #set.seed(666) # Also not bad, but Rv0512 is not in rank, but it is pykA
        #set.seed(777) # Works well, but Rv0512 doesn't appear in rank. Try again
        trainFitDels <- train(linGroup ~ . ,
                              data = delsBinRdx,
                              method = "rf",
                              metric = metricRF,
                              trControl = trControl, 
                              weights = modelWeights)
}else{
        #set.seed(123) # Already tried, same
        set.seed(999)
        #set.seed(666) # Also not bad, but Rv0512 is not in rank, but it is pykA
        #set.seed(777) # Works well, but Rv0512 doesn't appear in rank. Try again
        trainFitDels <- train(linGroup ~ . ,
                              data = delsBinRdx,
                              method = "rf",
                              metric = metricRF,
                              trControl = trControl)
}


trainFitFile <- paste(c("trainFit",
                        outName,
                        ".RData"), 
                      collapse = "")
                         
save(trainFitDels, file = paste(resuDir, trainFitFile, sep = ""))
print(paste(trainFitFile, "saved at", resuDir))

predDels <- predict(trainFitDels, delsBinRdx)

confMatDels <- confusionMatrix(predDels, delsBinRdx$linGroup)

tableDels <- data.frame(confMatDels$table)

plotTableDels <- tableDels %>%
        mutate(goodbad = ifelse(tableDels$Prediction == tableDels$Reference, "good", "bad")) %>%
        group_by(Reference) %>%
        mutate(prop = Freq/sum(Freq))

# fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups (see dplyr code above as well as original confusion matrix for comparison)
# If plots flag is present output confusion matrix
if(parsed$plots == T){
        confMatFile <- paste("confMat", outName, ".pdf", sep = "")
        #pdf(paste(resuDir, confMatFile, sep = ""), width = 7, height = 6, onefile = T)
        p <- ggplot(data = plotTableDels, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = prop)) +
                     geom_tile() +
                     geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1, size = 20) +
                     scale_fill_manual(values = c(good = "green", bad = "red")) +
                     theme_bw() +
                     xlim(rev(levels(tableDels$Reference)))
        ggsave(paste(resuDir, confMatFile, sep = ""), 
               plot = p, 
               width = 7, 
               height = 6,
               device = "pdf")
        #dev.off()
        print(paste(confMatFile, "saved at", resuDir))
}


RF_Dels_acc <- confMatDels$overall[1]
RF_Dels_NIR <- confMatDels$overall[5]
RF_Dels_PVal <- confMatDels$overall[6]

print("Parameters of Random Forest classifying as animal or human lineages:")
print(paste("Accuracy:", RF_Dels_acc))
print(paste("NIR:", RF_Dels_NIR))
print(paste("p-value [Acc > 0.NIR]:", RF_Dels_PVal))

# Parameters from cross-validated model:
bestTry <- trainFitDels$bestTune[, 1]
print(sprintf("mtry finally used in the model: %s.", bestTry))
fitResults <- trainFitDels$results[trainFitDels$results$mtry == bestTry, ]
print(sprintf("Accuracy CV: %s.", fitResults$Accuracy))
print(sprintf("Kappa CV: %s.", fitResults$Kappa))
print(sprintf("Accuracy SD CV: %s.", fitResults$AccuracySD))
print(sprintf("Kappa SD CV: %s.", fitResults$KappaSD))
# Get variable importance. If chosen, plot barplot of top X important variables. 
#########################################################################################################
# Create varImp dataframe
varImpDels <- varImp(trainFitDels)

varImp_DF <- data.frame(gene = rownames(varImpDels$importance)[order(varImpDels$importance, 
                                                                     decreasing = T)],
                        Overall = varImpDels$importance[order(varImpDels$importance, 
                                                                 decreasing = T), ])
rownames(varImp_DF) <- rownames(varImpDels$importance)[order(varImpDels$importance, decreasing = T)]

varImpFile <- paste("varImp", outName, "DF.csv", sep = "")
write.csv(varImp_DF, file = paste(resuDir, varImpFile, sep = ""))
print(paste(varImpFile, "saved at", resuDir))

# Plot top X important variables if plot flag is present
if(parsed$plots == T){
        varImpPlotFile <- paste("varImp", outName, ".pdf", sep = "")
        #pdf(paste(resuDir, varImpPlotFile, sep = ""), height = 4, width = 2.2)
        p <- ggplot(data = head(varImp_DF, as.numeric(parsed$top)),
                     aes(x = reorder(gene, Overall), y = Overall)) +
                     labs(x = "Genes", y = "Variable Importance") +
                     geom_bar(stat = "identity") +
                     coord_flip()
        #dev.off()
        ggsave(paste(resuDir, varImpPlotFile, sep = ""), 
               plot = p,
               height = 4, 
               width = 2.2,
               device = "pdf")
        print(paste(varImpPlotFile, "saved at", resuDir))
}

topGenes <- varImp_DF$gene[1:parsed$top]

topGenesFile <- paste(outName, "TopGenes.txt", sep = "")

# Write text file with top X important genes (default = 20) in classification. 
write.table(topGenes, 
            file = paste(resuDir, 
                         topGenesFile, 
                         sep = ""), 
            sep = "\t", 
            quote = F, 
            row.names = F, 
            col.names = F)

print(paste(topGenesFile, "saved at", resuDir))

print("Done!")