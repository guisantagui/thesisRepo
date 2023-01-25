####################################################################################################
####################################################################################################
#                                                                                                  #
# Impute the missing values of the parsed metabolite data, log transform it and test the           #
# performance of different normalization methods, finally choosing CCMN and exporting it as a csv. #
#                                                                                                  #
####################################################################################################
####################################################################################################

# Download NormalizeMets from https://cran.r-project.org/src/contrib/Archive/NormalizeMets/
if(!require(NormalizeMets)) install.packages("C:/Users/Guillem/Documents/PhD/comput/NormalizeMets_0.25.tar.gz",
                                             repos = NULL, 
                                             type = "source")
library(NormalizeMets)
if(!require(cluster)) install.packages("cluster")
library(cluster)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
#detach("package:metabolomics", unload = T)
#library(metabolomics)
# Directory stuff
####################################################################################################
rootDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/PA_paper_scripts/" ## Change this 
                                                                                      ## path 
                                                                                      ## according 
                                                                                      ## user's 
                                                                                      ## direcory
                                                                                      ## location.
dataDir <- paste0(rootDir, "data/")
outDir <- paste0(rootDir, "results/")
plotDir <- paste0(rootDir, "plots/")
figDir <- paste0(plotDir, "figS5/")
preProcDir <- paste0(plotDir, "preProc/")

# Create directories if necessary
if(!dir.exists(plotDir)){
        dir.create(plotDir)
}
if(!dir.exists(outDir)){
        dir.create(outDir)
}
if(!dir.exists(figDir)){
        dir.create(figDir)
}
if(!dir.exists(preProcDir)){
        dir.create(preProcDir)
}

####################################################################################################
#                                                                                                  #
# Load the data.                                                                                   #
#                                                                                                  #
####################################################################################################


# Dictionary of metabolite names between different metabolite sets
dictionary <- read.csv(file = paste0(outDir, "dictionary_solvedPut.csv"), 
                       row.names = 1, 
                       stringsAsFactors = F)

# Parsed metabolomic dataset
mets_parsed <- read.csv(file = paste0(outDir, "mets_peakRaw_parsed.csv"), 
                        row.names = 1,
                        stringsAsFactors = F, 
                        header = T)
# When csv is imported R changes the column names. Correct it with dictionary
colnames(mets_parsed)[2:ncol(mets_parsed)] <- dictionary$Consensus[match(colnames(mets_parsed)[2:ncol(mets_parsed)],
                                                                         make.names(dictionary$Consensus))]

# Functions for assess best normalization method. 
####################################################################################################

metvar <- function(valuemat, groups) {
        ns <- nrow(valuemat)
        nmet <- ncol(valuemat)
        unig <- unique(groups)
        ng <- length(unig)
        meanvalues <- matrix(nrow=ng, ncol=nmet)
        row.names(meanvalues) <- unig
        colnames(meanvalues) <- colnames(valuemat)
        meandifvalues <- matrix(nrow=ns, ncol=nmet)
        meancvvalues <- matrix(nrow=ns,ncol=nmet)
        for (i in 1:ng) {
                meanvalues[i, ] <- colMeans(valuemat[groups == unig[i], ], na.rm = T)
                for (j in 1:nmet) {
                        meandifvalues[groups == unig[i], j] <-
                                valuemat[groups == unig[i], j] - meanvalues[i, j]
                        meancvvalues[groups == unig[i], j] <-
                                (valuemat[groups == unig[i], j] - meanvalues[i, j])/meanvalues[i, j]
                }
        }
        
        ssr <- vector()
        sst <- vector()
        ssbet <- vector()
        n1 <- vector()
        n2 <- vector()
        fstat <- vector()
        kruskals <- vector()
        kruskalp <- vector()
        kruskalpadj <- vector()
        metmean <- colMeans(valuemat, na.rm = T)
        metnames <- colnames(valuemat)
        rangefc <- vector()
        
        for (i in 1:nmet) {
                ssr[i] <- sum(meandifvalues[, i]^2, na.rm = T)
                sst[i] <- sum((valuemat[, i] - metmean[i])^2, na.rm = T)
                ssbet[i] <- sst[i] - ssr[i]
                n1[i] <- sum(!is.na(meanvalues[, i])) - 1
                n2[i] <- sum(!is.na(valuemat[, i])) - n1[i]
                fstat[i] <- (ssbet[i]/n1[i])/(ssr[i]/n2[i])
                datavec <- valuemat[, i]
                groupvec <- groups[!is.na(datavec)]
                datavec <- datavec[!is.na(datavec)]
                kresult <- kruskal.test(datavec,groupvec)
                kruskalp[i] <- kresult$p.value
                kruskals[i] <- kresult$statistic
                rangefc[i] <- 10^(max(meanvalues[, i]) - min(meanvalues[, i]))
        }
        kruskalpadj <- p.adjust(kruskalp, method = "BH")
        metvariation <- data.frame(names = metnames, ssr, ssbet, sst,
                                   n1, n2, f = fstat, k = kruskals, p = kruskalp, 
                                   padj = kruskalpadj, mean = metmean, range = rangefc)
        
        
        list(metvar=metvariation, groupmeans=meanvalues, residues=meandifvalues, cv=meancvvalues)
        
}

# Computes silhouette, intra-group average euclidean distance, inter-group average euclidean 
# distance, and other statistics that weren't finally used in the plots. 
neibcheck <- function(valuemat) {
        gclass <- c(1,1,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,
                    11,11,12,12,12,13,13,13,14,14,14,15,15,15,16,16,16,17,17,17,18,18,18,19,
                    19,19,20,20,20,21,21,21,22,22,22,23,23,23,24,24,24,25,25,25,26,26,26,27,27,27)
        gclass <- as.numeric(factor(gsub("\\_.*",
                                         "",
                                         rownames(logdataMet$featuredata)),
                                    levels = unique(gsub("\\_.*",
                                                         "",
                                                         rownames(logdataMet$featuredata)))))
        straindist <- as.matrix(dist(valuemat, method="euclidean"))
        score <- rep(0, 28)
        rscore <- rep(0, 28)
        silh <- silhouette(gclass, dmatrix = straindist)
        mgclass <- matrix(rep(gclass, 84), nrow=84, ncol=84)
        tmgclass <- t(mgclass)
        diffclass <- (mgclass != tmgclass)
        sameclassdist <- mean(straindist[!diffclass])
        difclassdist <- mean(straindist[diffclass])
        for (i in 1:28) {
                subdist <- straindist[(1+3*(i-1)):(3*i), ]
                points <- 0
                rpoints <- 0
                for (j in 1:3){
                        rankvec <- rank(subdist[j, ])
                        dvec <- subdist[j, ]
                        mvec <- dvec
                        mvec[dvec == 0] <- max(dvec)
                        dvec[dvec == 0] <- min(mvec)
                        dvec <- (dvec - min(mvec))/(max(dvec) - min(mvec))
                        points <- points + sum(rankvec[(1+3*(i-1)):(3*i)]) - 6
                        rpoints <- rpoints + sum(dvec[(1+3*(i-1)):(3*i)])
                }
                score[i] <- points
                rscore[i] <- rpoints
        }
        score <- score/3
        rscore <- rscore/6
        list(score = score, rscore = rscore, silh = silh[,3],
             sdist = sameclassdist, ddist = difclassdist)
}

neibcheck <- function(valuemat) {
        strain <- gsub("\\_.*|(PA14).*",
                       "\\1",
                       rownames(valuemat))
        gclass <- as.numeric(factor(strain,
                                    levels = unique(strain)))
        straindist <- as.matrix(dist(valuemat, method="euclidean"))
        score <- rep(0, length(unique(gsub("\\_.*", "", rownames(valuemat)))))
        rscore <- rep(0, length(unique(gsub("\\_.*", "", rownames(valuemat)))))
        silh <- silhouette(gclass, dmatrix = straindist)
        mgclass <- matrix(rep(gclass, nrow(valuemat)),
                          nrow=nrow(valuemat),
                          ncol=nrow(valuemat))
        tmgclass <- t(mgclass)
        diffclass <- (mgclass != tmgclass)
        dimnames(diffclass) <- list(rownames(valuemat),
                                    rownames(valuemat))
        sameclassdist <- mean(straindist[!diffclass])
        difclassdist <- mean(straindist[diffclass])
        for (i in 1:28) {
                subdist <- straindist[(1+3*(i-1)):(3*i), ]
                points <- 0
                rpoints <- 0
                for (j in 1:3){
                        rankvec <- rank(subdist[j, ])
                        dvec <- subdist[j, ]
                        mvec <- dvec
                        mvec[dvec == 0] <- max(dvec)
                        dvec[dvec == 0] <- min(mvec)
                        dvec <- (dvec - min(mvec))/(max(dvec) - min(mvec))
                        points <- points + sum(rankvec[(1+3*(i-1)):(3*i)]) - 6
                        rpoints <- rpoints + sum(dvec[(1+3*(i-1)):(3*i)])
                }
                score[i] <- points
                rscore[i] <- rpoints
        }
        score <- score/3
        rscore <- rscore/6
        list(score = score, rscore = rscore, silh = silh[,3],
             sdist = sameclassdist, ddist = difclassdist)
}


# Do P-value histogram of the kruskal-wallis tests computed within metvar function
doPValHist <- function(metVar, normMethod){
        plot <- qplot(metVar$metvar$padj, geom = "histogram",
              main = normMethod, 
              xlab = "p-value",
              ylab = "Frequency",
              fill = I("black"),
              alpha = I(.5), 
              xlim = c(-0.0106, 0.6),
              ylim = c(0, 85)) +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"))
        return(plot)
}

# Do within or across groups RLA plot
doRLAPlot <- function(featData, whichRLA = NA, normMeth = ""){
        featData <- data.frame(featData)
        stop_quietly <- function() {
                opt <- options(show.error.messages = FALSE)
                on.exit(options(opt))
                stop()
        }
        whichRLA <- tolower(whichRLA)
        if(!(whichRLA == "wg" | whichRLA == "ag")){
                print("Please specify if RLA plot is within-group (WG) or across-group (AG)")
                stop_quietly()
        }
        if(whichRLA == "ag"){
                if(nchar(normMeth) > 0){
                        main <- sprintf("%s, across-group RLA plot", normMeth)
                }else{
                        main <- "Across-group RLA plot"
                }
                featData <- data.frame(apply(featData, 2, function(x) (x - median(x))/sd(x)))
        }
        abuDF <- data.frame(matrix(nrow = 0, 
                                   ncol = 4, 
                                   dimnames = list(NULL,
                                                   c("sample",
                                                     "strain",
                                                     "metabolite",
                                                     "abundance"))))
        for(i in 1:nrow(featData)){
                samp <- rownames(featData)[i]
                sampDF <- data.frame(t(featData[i, ]))
                colnames(sampDF) <- "abundance"
                sampDF$metabolite <- rownames(sampDF)
                sampDF$sample <- rep(samp, nrow(sampDF))
                sampDF$strain <- gsub("\\_.*", "", sampDF$sample)
                sampDF <- sampDF[, c("sample", "strain", "metabolite", "abundance")]
                abuDF <- rbind.data.frame(abuDF, sampDF)
        }
        rownames(abuDF) <- 1:nrow(abuDF)
        if(whichRLA == "wg"){
                for(strain in unique(abuDF$strain)){
                        strainDF <- abuDF[abuDF$strain == strain, ]
                        for(met in unique(strainDF$metabolite)){
                                abStrainMet <- strainDF$abundance[strainDF$metabolite == met]
                                meanMet <- median(abStrainMet)
                                sdMet <- sd(abStrainMet)
                                abuDF$abundance[abuDF$strain == strain & abuDF$metabolite == met] <- sapply(abuDF$abundance[abuDF$strain == strain & abuDF$metabolite == met],
                                                                                                            function(x) (x - meanMet))
                        }
                }
                if(nchar(normMeth) > 1){
                        main <- sprintf("%s, within-group RLA plot", normMeth)
                }else{
                        main <- "Within-group RLA plot"
                }
        }
        strainCols <- rainbow(length(unique(abuDF$strain)))
        plot <- ggplot(data = abuDF, mapping = aes(x = sample, y = abundance, fill = strain)) +
                geom_boxplot() +
                ylim(c(-2, 2)) +
                scale_discrete_manual(aesthetics = "fill",
                                      values = strainCols) +
                labs(title = main, y = "relative log abundance") +
                theme(title = ggtext::element_markdown(),
                      axis.text.x = element_text(angle = 90),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "none")
        plot
}

# Do PCA, colored according to batch
plotBatchPCA <- function(PC, x = "PC1", y = "PC2", main = NULL, labels = F, legend = T){
        batchDF <- read.csv(file = paste0(outDir, "batchDF.csv"),
                            row.names = 1, 
                            stringsAsFactors = F)
        data <- data.frame(obsnames=row.names(PC$x), PC$x)
        data <- data[, c("obsnames", x, y)]
        data$batch <- batchDF$batch[match(gsub("\\_.*", "", data$obsnames),
                                          batchDF$strain)]
        propVar <- summary(PC)$importance[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)
        
        plot <- ggplot(data, aes(x = data[, x], 
                                 y = data[, y], 
                                 color = batch,
                                 label = obsnames)) + 
                geom_hline(yintercept = 0, alpha = 0.6) +
                geom_vline(xintercept = 0, alpha = 0.6) +
                geom_point() + 
                labs(title = main,
                     x = sprintf("%s (%s%%)", x, propX),
                     y = sprintf("%s (%s%%)", y, propY))
                
        if(labels){
                plot <- plot +
                        geom_text_repel()
        }
        if(legend){
                plot <- plot + 
                        theme(title = ggtext::element_markdown(),
                              axis.title.y = ggtext::element_markdown(),
                              panel.background = element_blank(),
                              panel.border = element_rect(colour = "black", fill=NA, size=1),
                              panel.grid.major = element_line(colour = "#d4d4d4"),
                              legend.position = "right")
        }else{
                plot <- plot +
                        theme(title = ggtext::element_markdown(),
                              axis.title.y = ggtext::element_markdown(),
                              panel.background = element_blank(),
                              panel.border = element_rect(colour = "black", fill=NA, size=1),
                              panel.grid.major = element_line(colour = "#d4d4d4"),
                              legend.position = "none")
        }
        plot
}
####################################################################################################
# Do imputation of missing values (NAs).                                                           #
####################################################################################################

# Do scatter plot of metabolite quality according to NA proportion
####################################################################################################
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

metQual <- data.frame("metabolite" = colnames(mets_parsed[2:ncol(mets_parsed)]),
                      "missingProp" = apply(mets_parsed[, 2:ncol(mets_parsed)], 
                                            2, 
                                            function(x) sum(is.na(x))/nrow(mets_parsed)),
                      "rangePeaks" = apply(mets_parsed[, 2:ncol(mets_parsed)],
                                           2,
                                           function(x) max(x[!is.na(x)]) - min(x[!is.na(x)])),
                      stringsAsFactors = F)

# Change the names of the metabolites for the figure to the ones suggested by Kyu Rhee or, 
# in the case that it's a NA, to the one indicated in the "definitiveNames" column of 
# metabolite dictionary
metQualNewNames <- dictionary$KyuRheeNames[match(metQual$metabolite, dictionary$Consensus)]
metQualNewNames[is.na(metQualNewNames)] <- dictionary$definitiveNames[match(metQual$metabolite, 
                                                                            dictionary$Consensus)][is.na(metQualNewNames)]

metQual$metabolite <- metQualNewNames


# Plot a scatterplot of max value-min value (peak range) in front of Missing value
# proportion (NA proportion).
ggplot(data = metQual, aes(x = missingProp, y = rangePeaks)) + geom_point() +
        geom_text(label=metQual$metabolite, nudge_x  = 0.04, nudge_y = 1) + 
        ylab("Peak range") + 
        xlab("Missing Value Proportion") + 
        theme_minimal()
ggsave(filename = paste0(figDir, "figS5_B.1.pdf"), 
       width = 6, 
       height = 6)


# Check what are the compounds that have NAs.
####################################################################################################
strains <- gsub("\\_.*|(PA14).*", 
                rownames(mets_parsed), 
                rep = "\\1")

hasNAs <- colnames(mets_parsed)[2:ncol(mets_parsed)][which(apply(mets_parsed[, 2:ncol(mets_parsed)], 
                                                                 2, 
                                                                 function(x) sum(is.na(x))) > 0)]

# Impute the values that are missing in all replicates by NA, and the rest
# by other replicates' average value.
####################################################################################################
mets_parsImp <- mets_parsed
for(i in seq_along(hasNAs)){
        impMet <- mets_parsed[[hasNAs[i]]]
        for(j in seq_along(strains)){
                strainSel <- gsub("\\_.*|(PA14).*",
                                  rownames(mets_parsed),
                                  rep = "\\1") %in% strains[j]
                subVec <- impMet[strainSel]
                if(sum(is.na(subVec)) == length(subVec)){
                        impMet[strainSel] <- rep(NA, length(subVec))
                }else{
                        meanImp <- mean(subVec[!is.na(subVec)])
                        impMet[strainSel] <- rep(meanImp, length(subVec))
                }
        }
        mets_parsImp[[hasNAs[i]]] <- impMet 
}

metQualPostImp <- data.frame("metabolite" = colnames(mets_parsImp[2:ncol(mets_parsImp)]),
                             "missingProp" = apply(mets_parsImp[, 2:ncol(mets_parsImp)], 
                                                   2, 
                                                   function(x) sum(is.na(x))/nrow(mets_parsImp)),
                             "rangePeaks" = apply(mets_parsImp[, 2:ncol(mets_parsImp)],
                                                  2,
                                                  function(x) max(x[!is.na(x)]) - min(x[!is.na(x)])),
                             stringsAsFactors = F)

metQualPostImp$metabolite <- metQualNewNames

ggplot(data = metQualPostImp, aes(x = missingProp, y = rangePeaks)) + geom_point() +
        geom_text(label=metQualPostImp$metabolite, nudge_x  = 0.04, nudge_y = 1) + 
        ylab("Peak range") + 
        xlab("Missing Value Proportion") + 
        theme_minimal()
ggsave(filename = paste0(figDir, "figS5_B.2.pdf"), width = 6, height = 6)

# Remove all compounds that after imputation still have any NAs
####################################################################################################

# Compute NA proportion per compound
impNAProp <- apply(mets_parsImp, 2, function(x) sum(is.na(x))/length(x))

# remove metabolites with a proportion of NAs > 0.00
toRemove <- names(impNAProp[which(impNAProp > 0.00)])

mets_parsImp <- mets_parsImp[, !colnames(mets_parsImp) %in% toRemove]

# Obtain the set of Internal Standards for normalization
####################################################################################################

# Use metvar function to do several statistical tests
metvarprenorm <- metvar(mets_parsImp[, 2:ncol(mets_parsImp)],
                        mets_parsImp[, 1])


# Use the padj output from metvar function to get the metabolites that are constant accross
# all the samples. This is the result of a Kruskal Wallis test between each strain replicates, 
# adjusted with Benjamini-Hochberg method.
IS <- c(rep(0, length(metvarprenorm$metvar$padj)))
for(i in 1:length(metvarprenorm$metvar$padj)){
        if(metvarprenorm$metvar$padj[i]>0.05){
                IS[i] <- 1
        }
}

# Prepare alldata object for NormalizeMets package
featuredataMet <- mets_parsImp[, 2:ncol(mets_parsImp)]
sampledataMet <- data.frame(Group = mets_parsImp[, 1], 
                            Species = c(rep("P. aeruginosa", length(mets_parsImp[, 1]))))


metabolitedataMet <- data.frame(names = colnames(featuredataMet), IS = IS)

alldataMet <- list(featuredata = featuredataMet, sampledata = sampledataMet, 
                   metabolitedata = metabolitedataMet)

unNormMetImpData <- featuredataMet

# save(unNormMetImpData, file = "unNormMetImpData.RData")

####################################################################################################
#                                                                                                  #
# LogTransform Data and plot log abundances (Fig S5 A).                                            #
#                                                                                                  #
####################################################################################################

# Log transform and add sample data and metabolite data to create the alldata object type for 
# the NormalizePackage input
logdataMet<- LogTransform(alldataMet$featuredata, zerotona = T)
logdataMet$sampledata <- sampledataMet
logdataMet$sampledata <- sampledataMet
logdataMet$metabolitedata <- metabolitedataMet

# Prepare log transformed dataframe for plot
metnames=colnames(logdataMet$featuredata)
samplename=rownames(logdataMet$featuredata)
logdatatable=data.frame(sample=rep('.',length(samplename)*length(metnames)),
                        strain=rep('.',length(samplename)*length(metnames)),
                        metabolite=rep('.',length(samplename)*length(metnames)),
                        peak=rep(0,length(samplename)*length(metnames)),
                        significant=rep(0,length(samplename)*length(metnames)),
                        stringsAsFactors = F)
line=0
for (i in 1:length(samplename)){
        for (j in seq_along(metnames)){
                line=line+1
                logdatatable[line,1]=samplename[i]
                logdatatable[line,2]=as.character(logdataMet$featuredata[i,1])
                logdatatable[line,3]=metnames[j]
                logdatatable[line,4]=logdataMet$featuredata[i,j]
                if (metvarprenorm$metvar$padj[metvarprenorm$metvar$names==metnames[j]]<0.05){
                        logdatatable[line,5]=1
                }
        }
}


# Change names to the ones suggested by Kyu Rhee and add the information for the plot. 
logdatatable_newNames <- dictionary$KyuRheeNames[match(logdatatable$metabolite, dictionary$Consensus)]
logdatatable$metabolite[!is.na(logdatatable_newNames)] <- logdatatable_newNames[!is.na(logdatatable_newNames)]
logdatatable$ambiguous <- rep("Non-Ambiguous", length(logdatatable$metabolite))
logdatatable$ambiguous[grep("_?", dictionary$Consensus[match(logdatatable$metabolite, dictionary$KyuRheeNames)], 
                            fixed = T)] <- "Ambiguous"
logdatatable$legend <- logdatatable$ambiguous
logdatatable$legend[logdatatable$significant == 0] <- "Internal Standards"
logdatatable$legend[logdatatable$legend == "Non-Ambiguous"] <- "Variable Metabolites"
logdatatable$legend[logdatatable$legend == "Internal Standards"] <- "Non-changing Metabolites (IS)"
logdatatable$legend[logdatatable$legend == "Ambiguous"] <- "Putative Metabolites"
logdatatable$legend <- as.factor(logdatatable$legend)

# Do plot 
metBarplot <- ggplot(data = logdatatable, mapping = aes(x = reorder(metabolite, 
                                                                    peak, 
                                                                    FUN=function(x) median(x[!is.na(x)])), 
                                                        y = peak, 
                                                        color=legend)) +
        geom_boxplot() +
        coord_flip() +
        labs(
                y= "Log Abundance",
                x= "Metabolites"
                ) +
        scale_color_manual(values = c("#000000", "red", "#33B0FF")) +
        scale_y_reverse() +
        theme(axis.text.y = element_text(angle = 180, hjust = 0, size = 11),
              axis.text.x = element_text(angle = 90, hjust = 1, size =  11),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.grid.major = element_line(colour = "#d4d4d4"))

metBarplot

ggsave(paste0(figDir, "figS5_A.pdf"), width = 6.5, height = 10.2)  

# Prepare 
#missingMet <-  MissingValues(logdataMet$featuredata, sampledataMet,metabolitedataMet,
#                             feature.cutof=0.8, sample.cutoff=0.2, method="knn")


####################################################################################################
#                                                                                                  #
# Assess what is the best berforming normalization method amomng CCMN, NOMIS, RUV Random, IS,      #
# scale to median, scale to mean and scale to sum.                                                 #
#                                                                                                  #
####################################################################################################

# Use neibcheck function to compute silhouette, intra group and inter group euclidean distances 
# to assess best normalization method
neibprenorm <- neibcheck(logdataMet$featuredata)

# CCMN
ccmnNormMets <- NormQcmets(logdataMet$featuredata, 
                           method = "ccmn", 
                           qcmets = which(logdataMet$metabolitedata$IS == 1), 
                           factors = logdataMet$sampledata$Group,
                           saveoutput = T)
metvarCcmnNorm <- metvar(ccmnNormMets$featuredata, 
                         logdataMet$sampledata$Group)
neibCcmnNorm <- neibcheck(ccmnNormMets$featuredata)

# NOMIS
nomisNormMets <- NormQcmets(logdataMet$featuredata, 
                            method = "nomis", 
                            qcmets = which(logdataMet$metabolitedata$IS == 1))
metvarNomisNorm <- metvar(nomisNormMets$featuredata, 
                          logdataMet$sampledata$Group)
neibNomisNorm <- neibcheck(nomisNormMets$featuredata)

# RUV Random
ruvRandNormMets <- NormQcmets(logdataMet$featuredata, 
                              method = "ruvrand", 
                              k = 2, 
                              qcmets = which(logdataMet$metabolitedata$IS == 1),
                              plotk = T)
metvarRuvRandNorm <- metvar(ruvRandNormMets$featuredata, 
                            logdataMet$sampledata$Group)
neibRuvRandNorm <- neibcheck(ruvRandNormMets$featuredata)

# RUV Random for clustering
ruvRandClustNormMets <- NormQcmets(logdataMet$featuredata, method = "ruvrandclust", k = 2, 
                                   qcmets = which(logdataMet$metabolitedata$IS == 1))

metvarRuvRandClustNorm <- metvar(ruvRandClustNormMets$featuredata, 
                            logdataMet$sampledata$Group)
neibRuvRandClustNorm <- neibcheck(ruvRandClustNormMets$featuredata)

# IS
isNormMets <- NormQcmets(logdataMet$featuredata, 
                         method = "is", 
                         isvec = which(logdataMet$metabolitedata$IS == 1))#,
                         #isvec = logdataMet$featuredata[, which(logdataMet$metabolitedata$IS == 1)[1]])
metvarIsNorm <- metvar(isNormMets$featuredata,
                       logdataMet$sampledata$Group)
neibIsNorm <- neibcheck(isNormMets$featuredata)

# Scale to median
medianNormMets <- NormScaling(logdataMet$featuredata, 
                              method = "median")
metvarMedianNorm <- metvar(medianNormMets$featuredata, 
                           logdataMet$sampledata$Group)
neibMedianNorm <- neibcheck(medianNormMets$featuredata)

# Scale to mean
meanNormMets <- NormScaling(logdataMet$featuredata, 
                            method = "mean")
metvarMeanNorm <- metvar(meanNormMets$featuredata, 
                         logdataMet$sampledata$Group)
neibMeanNorm <- neibcheck(meanNormMets$featuredata)

# Scale to sum
sumNormMets <- NormScaling(logdataMet$featuredata, 
                           method = "sum")
metvarSumNorm <- metvar(sumNormMets$featuredata, 
                        logdataMet$sampledata$Group)
neibSumNorm <- neibcheck(sumNormMets$featuredata)

# Do list of the metvar/neib outputs for normalize method 
# report
normRepList <- list(preNorm = list(metvar = metvarprenorm,
                                   neib = neibprenorm),
                    nomis = list(metvar = metvarNomisNorm,
                                 neib = neibNomisNorm),
                    ccmn = list(metvar = metvarCcmnNorm,
                                neib = neibCcmnNorm),
                    ruvRand = list(metvar = metvarRuvRandNorm,
                                   neib = neibRuvRandNorm),
                    ruvRandClust = list(metvar = metvarRuvRandClustNorm,
                                        neib = neibRuvRandClustNorm),
                    is = list(metvar = metvarIsNorm,
                              neib = neibIsNorm),
                    median = list(metvar = metvarMedianNorm,
                                  neib = neibMedianNorm),
                    mean = list(metvar = metvarMeanNorm,
                                neib = neibMeanNorm),
                    sum = list(metvar = metvarSumNorm,
                               neib = neibSumNorm))

normMethDF <- data.frame(id = names(normRepList),
                         name = c("Pre-normalization",
                                  "NOMIS",
                                  "CCMN",
                                  "RUV-random",
                                  "RUV-random clust",
                                  "IS",
                                  "Median",
                                  "Mean",
                                  "Sum"),
                         stringsAsFactors = F)
# Do a dataframe with previous list
normRepDF <- data.frame(method = names(normRepList), 
                        silh = rep(0, length(normRepList)),
                        score = rep(0, length(normRepList)), 
                        rscore = rep(0, length(normRepList)), 
                        difmet = rep(0, length(normRepList)), 
                        sd = rep(0, length(normRepList)), 
                        dd = rep(0, length(normRepList)), 
                        stringsAsFactors = F)

# Do a list of the outputs of normalization methods
normOutsList <- list(preNorm = logdataMet,
                     nomis = nomisNormMets,
                     ccmn = ccmnNormMets,
                     ruvRand = ruvRandNormMets,
                     ruvRandClust = ruvRandClustNormMets,
                     is = isNormMets,
                     median = medianNormMets,
                     mean = meanNormMets,
                     sum = sumNormMets)

for(i in seq_along(normRepList)){
        normRepDF$silh[i] <- mean(normRepList[[i]]$neib$silh)
        normRepDF$score[i] <- mean(normRepList[[i]]$neib$score)
        normRepDF$rscore[i] <- mean(normRepList[[i]]$neib$silh)
        normRepDF$difmet[i] <- sum(normRepList[[i]]$metvar$metvar$padj < 0.01)
        normRepDF$sd[i] <- normRepList[[i]]$neib$sdist
        normRepDF$dd[i] <- normRepList[[i]]$neib$ddist
}

# Plot mean intra-group distance vs mean inter-group distance for each normalization method
####################################################################################################

ggplot(data = normRepDF, mapping = aes(x = sd, y = dd)) +
        geom_point() +
        geom_label(aes(label=method),nudge_y=0.025) +
        labs(
                y= "Mean distance between strains",
                x= "Mean distance between replicates"
        )

# IS and sum methods produce intra group and inter group mean distances very high.
# They are not working well and make the rest of the methods appear all stacked in the plot. 
# Remove them from plot
normRepDF_4Plot <- normRepDF[!normRepDF$method %in% c("sum", "is"),]

normRepDF_4Plot$method <- normMethDF$name[match(normRepDF_4Plot$method,
                                                normMethDF$id)]

ggplot(data = normRepDF_4Plot, mapping = aes(x = sd, y = dd)) +
        geom_point() +
        geom_label(aes(label=method),nudge_y=0.025) +
        labs(
                y = "Mean distance between strains",
                x = "Mean distance between replicates"
        ) +
        theme(title = ggtext::element_markdown(),
                axis.title.y = ggtext::element_markdown(),
                panel.background = element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                panel.grid.major = element_line(colour = "#d4d4d4")) +
        xlim(c(1.6, 2.5))

ggsave(paste0(preProcDir, "normMet_avDists.pdf"), 
       height = 1000, 
       width = 1000, 
       units = "px")

# Plot silhouette value boxplots for each normalization method
####################################################################################################

silDF <- data.frame(matrix(ncol = 2, 
                           nrow = 0, 
                           dimnames = list(NULL, 
                                           c("method", "silh"))))
for(i in seq_along(normRepList)){
        normMeth <- names(normRepList)[i]
        silh <- normRepList[[i]]$neib$silh
        toBind <- data.frame(method = rep(normMeth, length(silh)),
                             silh = silh)
        silDF <- rbind.data.frame(silDF,
                                  toBind)
}

# Remove sum and is from the dataframe (because they don't work well). 
silDF <- silDF[!silDF$method %in% c("sum", "is"), ]
silDF$method <- normMethDF$name[match(silDF$method, normMethDF$id)]

ggplot(data = silDF, mapping = aes(x = reorder(method, silh, FUN=median), y = silh)) +
        geom_boxplot() +
        #geom_line(data=normreport,mapping=aes(x=method,y=silh,group=1)) +
        theme(axis.text.x = element_text(angle = 90)) +
        #coord_flip() +
        labs(
                y= "Silhouete",
                x= "Method"
        ) +
        theme(title = ggtext::element_markdown(),
              axis.title.y = ggtext::element_markdown(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.grid.major = element_line(colour = "#d4d4d4"))
ggsave(paste0(preProcDir, "normMet_silh.pdf"), 
       height = 1000, 
       width = 1000, 
       units = "px") 

# Plot adjusted p-value histograms of the Kruskal-Wallis tests that were run within metvar
####################################################################################################

doPValHist(metvarprenorm, "Pre-normalization")
ggsave(paste0(preProcDir, "pvaluehist_preNorm.pdf"), 
       height = 800, width = 800, units = "px")
doPValHist(metvarNomisNorm, "NOMIS")
ggsave(paste0(preProcDir, "pvaluehist_nomis.pdf"), 
       height = 800, width = 800, units = "px")
doPValHist(metvarCcmnNorm, "CCMN")
ggsave(paste0(preProcDir, "pvaluehist_ccmn.pdf"), 
       height = 800, width = 800, units = "px")
doPValHist(metvarRuvRandNorm, "RUV-random")
ggsave(paste0(preProcDir, "pvaluehist_ruvRand.pdf"), 
       height = 800, width = 800, units = "px")
doPValHist(metvarRuvRandClustNorm, "RUV-random clust")
ggsave(paste0(preProcDir, "pvaluehist_ruvRandClust.pdf"), 
       height = 800, width = 800, units = "px")
doPValHist(metvarMeanNorm, "Mean")
ggsave(paste0(preProcDir, "pvaluehist_mean.pdf"), 
       height = 800, width = 800, units = "px")
doPValHist(metvarMedianNorm, "Median")
ggsave(paste0(preProcDir, "pvaluehist_median.pdf"), 
       height = 800, width = 800, units = "px")

# Obtain Within-Group and across-Group Relative Log Abundance (RLA) plots for each normalization 
# method. 
# Do PCAs coloring according to batch before and after each normalization method too.
####################################################################################################
RLAPlots <- list()
batchPCAs <- list()
for(i in seq_along(normOutsList)){
        # RLA plots
        normMeth <- names(normOutsList)[i]
        normName <- normMethDF$name[match(normMeth, normMethDF$id)]
        RlaPlots(featuredata = normOutsList[[i]]$featuredata, 
                 groupdata = logdataMet$sampledata$Group, type = "wg",
                 interactiveplot = F, interactiveonly = F, 
                 cols = rainbow(28), minoutlier = 0.5, las = 3,
                 plotname = sprintf("%sRLA_WG_%s", preProcDir, normMeth),
                 saveplot = T, savetype = "pdf",
                 ylim = c(-2, 2))
        RlaPlots(featuredata = normOutsList[[i]]$featuredata, 
                 groupdata = logdataMet$sampledata$Group, type = "ag",
                 interactiveplot = F, interactiveonly = F, 
                 cols = rainbow(28), minoutlier = 0.5, las = 3,
                 plotname = sprintf("%sRLA_AG_%s", preProcDir, normMeth),
                 saveplot = T, savetype = "pdf",
                 ylim = c(-2, 2))
        # Custom function RLA plots
        wgRLA <- doRLAPlot(normOutsList[[i]]$featuredata, "WG",
                           normMeth = normName)
        wgRLA
        ggsave(filename = sprintf("%sRLA_wg_ggplt_%s.pdf",
                                  preProcDir,
                                  normMeth),
               height = 1000,
               width = 2600,
               units = "px")
        agRLA <- doRLAPlot(normOutsList[[i]]$featuredata, "AG",
                           normMeth = normName)
        agRLA
        ggsave(filename = sprintf("%sRLA_ag_ggplt_%s.pdf",
                                  preProcDir,
                                  normMeth),
               height = 1000,
               width = 2600,#2800
               units = "px")
        rlaP <- list(wg = wgRLA,
                     ag = agRLA)
        RLAPlots[[normMeth]] <- rlaP
        # Batch PCAs
        batchPCA <- plotBatchPCA(prcomp(normOutsList[[i]]$featuredata, 
                                        scale. = T, 
                                        center = T),
                                 main = normName,
                                 labels = F,
                                 legend = F)
        batchPCA
       
        ggsave(filename =  sprintf("%sbatchPCA_%s.pdf",
                                   preProcDir,
                                   normMeth),
               height = 3.2, width = 3)
        
        batchPCAs[[normMeth]] <- batchPCA
}

# CCMN method is the one that has a higher silhouette values, meaning that with this method
# the replicates are more similar to each other in relation to the other groups. It's also the 
# method that better minimizes the average distance between replicates while maximizes the 
# average distance between groups (strains). The within group and across group RLA plots also
# show that this method keeps the across group variability while making the replicates closer
# to each other, meaning that it's removing the unwanted variation. 

####################################################################################################
#                                                                                                  #
# Export CCMN-normalized metabolite dataframe.                                                     #
#                                                                                                  #
####################################################################################################


ccmn_export <- ccmnNormMets$featuredata
colnames(ccmn_export) <- colnames(mets_parsImp)[match(colnames(ccmn_export), 
                                                      make.names(colnames(mets_parsImp)))]

# write a .csv file...
write.csv(ccmn_export, 
          file = paste0(outDir, "mets_ccmn.csv"))

# And save also as .tsv file, as writing .cvs modifies slightly values. Not a big deal, 
# but makes metabolite names to appear in a slightly altered way in the further heatmap. 
write.table(format(ccmn_export, digits = 22), 
            file = paste0(outDir, "mets_ccmn.tsv"), 
            quote = F, 
            row.names = T, 
            sep = "\t")
