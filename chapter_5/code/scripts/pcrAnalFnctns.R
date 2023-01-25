# Functions needed for running plcABCD qPCR analysis
###########################################################################################################################

# Parse results files. 
parsePcrRes <- function(resFile){
        if(!require(readxl)) install.packages("readxl")
        library(readxl)
        results <- as.data.frame(read_xlsx(resFile, sheet = 3))
        colnames(results) <- make.names(colnames(results))
        results <- results[which(results$Block.Type == "Well"):nrow(results), ]
        colnames(results) <- make.names(results[1, ])
        results <- results[2:nrow(results), ]
        results <- results[!is.na(results$Sample.Name), ]
        rownames(results) <- 1:nrow(results)
        return(results)
}

# Orders oligos in increasing number, while keeping the primer pair order as
# sigA, Mb1784c, Rv2349c, Rv2350c and Rv2351c
orderOut <- function(input){
        uniqOligs <- sort(unique(input$oligo))
        primerOrd <- c("sigA", "Rv1106", "Mb1784c", "Rv2349c", "Rv2350c", "Rv2351c", "Rv2379c", "Rv2383c", "Rv3545c", "mcr7")
        resOrdered <- data.frame(matrix(nrow = 0,
                                        ncol = ncol(input),
                                        dimnames = list(NULL,
                                                        colnames(input))))
        for(olig in uniqOligs){
                ordSubMat <- input[input$oligo == olig, ]
                #print(ordSubMat)
                primOrdSub <- primerOrd[primerOrd %in% ordSubMat$primerPair]
                ordSubMat <- ordSubMat[match(primOrdSub, ordSubMat$primerPair), ]
                resOrdered <- rbind.data.frame(resOrdered, ordSubMat)
        }
        return(resOrdered)
}

parseSamps <- function(parsdRslts, filtD1 = T, filtD3 = F){
        colnames(parsdRslts) <- make.names(colnames(parsdRslts))
        samps <- parsdRslts[gsub("_", "", gsub("[[:digit:]]", "", gsub("\\ .*", "", parsdRslts$Sample.Name))) == "R", ]
        #print(samps)
        colnames(samps) <- make.names(colnames(samps))
        
        out <- data.frame(sampleName = samps$Sample.Name,
                          oligo = sapply(samps$Sample.Name, 
                                         function(x) strsplit(x, 
                                                              " ")[[1]][1]),
                          primerPair = sapply(samps$Sample.Name, 
                                              function(x) strsplit(x, 
                                                                   " ")[[1]][2]),
                          Ct = as.numeric(samps$CT))
        #print(out)
        out$oligo <- as.character(out$oligo)
        out$oligo <- gsub("\\_.*", "", out$oligo)
        out$oligo <- sapply(out$oligo, function(x) paste0("R", 
                                                          paste0(rep("0", 4 - nchar(gsub("R", "", x))), collapse = ""),
                                                          gsub("R", "", x)))
        out$sampleName <- as.character(out$sampleName)
        out$primerPair <- as.character(out$primerPair)
        out$n_num <- oligoInfo$STRAIN[match(out$oligo, oligoInfo$NAME)]
        out$lin <- oligoInfo$LINEAGE[match(out$oligo, oligoInfo$NAME)]
        out$host <- oligoInfo$HOST[match(out$oligo, oligoInfo$NAME)]
        out$daysPI <- oligoInfo$SAMPLING.DATE[match(out$oligo, oligoInfo$NAME)] - oligoInfo$INFECTION.DATE[match(out$oligo, oligoInfo$NAME)]
        if(filtD1 == T){
                out <- out[out$daysPI == 1, ]
        }else if(filtD3 == T){
                out <- out[out$daysPI == 3, ]
        }
        
        out <- orderOut(out)
        
        return(out)
}

# Calculate primer efficiency. 
getEff <- function(LM){
        slope <- LM$coefficients[2]
        ampFact <- 10^(-1/slope)
        eff <- (ampFact - 1) * 100
        names(ampFact) <- NULL
        names(eff) <- NULL
        out <- list(ampFact = ampFact, efficiency = eff)
        return(out)
}

# Clean pcr samples, choosing among duplicates the ones with best Ct (lowest) in sigA gene. 
cleanSamps <- function(resultsWDups, 
                       filtGenes = c("sigA", "Rv1106", "Mb1784c", "Rv2349c", "Rv2350c", "Rv2351c", "Rv2379c", "Rv2383c", "Rv3545c", "mcr7"), 
                       refThrshld = 34){
        unOligs <- unique(resultsWDups$oligo)
        resBestCtAll <- data.frame(matrix(nrow = 0, ncol = ncol(resultsWDups), dimnames = list(NULL, colnames(resultsWDups))))
        # Select sample with lowest Ct among duplicates
        for(olig in unOligs){
                resOlig <- resultsWDups[resultsWDups$oligo == olig, ]
                #print(resOlig)
                minIdx <- which.min(resOlig$Ct[resOlig$primerPair == "sigA"])
                minRowName <- gsub("\\ .*", "", row.names(resOlig[resOlig$primerPair == "sigA"])[minIdx])
                resBestCt <- resultsWDups[gsub("\\ .*", "", rownames(resultsWDups)) == minRowName, ]
                resBestCtAll <- rbind.data.frame(resBestCtAll, resBestCt)
        }
        resBestCtAll <- resBestCtAll[resBestCtAll$primerPair %in% filtGenes, ]
        resBestCtAll <- resBestCtAll[order(resBestCtAll$primerPair, decreasing = T), ]
        # Remove samples with Ct higher than 30 in sigA gene. 
        badQual <- resBestCtAll$oligo[resBestCtAll$primerPair == "sigA" & resBestCtAll$Ct > refThrshld]
        resBestCtAll <- resBestCtAll[!resBestCtAll$oligo %in% badQual, ]
        resBestCtAll <- orderOut(resBestCtAll)
        return(resBestCtAll)
}

remBovis <- function(df){
        noBov <- df[df$lin != "M. bovis", ]
        noBov <- noBov[tolower(noBov$lin) != "bovis", ]
        return(noBov)
}

# Obtain gene expression relative to sigA. 
getRelExp <- function(resDF, ref = "sigA"){
        genes <- unique(resDF$primerPair)
        genes <- genes[genes != ref]
        factRef <- effs[[ref]]$ampFact
        refDF <- resDF[resDF$primerPair == "sigA", ]
        outDF <- data.frame(matrix(nrow = 0, ncol = (ncol(resDF) + 1),
                                   dimnames = list(NULL,
                                                   c(colnames(refDF), "relExp"))))
        for(gene in genes){
                if(gene %in% names(effs)){
                        geneFact <- effs[[gene]]$ampFact
                }else{
                        geneFact <- mean(unlist(lapply(effs, function(x) x$ampFact)))
                }
                resDF_gene <- resDF[resDF$primerPair == gene, ]
                relExp <- geneFact^(-resDF_gene$Ct)/factRef^(-refDF$Ct[match(resDF_gene$oligo, refDF$oligo)])
                resDF_gene$relExp <- relExp
                outDF <- rbind.data.frame(outDF, resDF_gene)
        }
        refDF$relExp <- rep(1, nrow(refDF))
        outWRef <- rbind.data.frame(refDF, outDF)
        outDF <- orderOut(outDF)
        return(outDF)
}

# Filter to obtain best samples based on best signal in sigA gene. 
getBestSamps <- function(refDF, n = 4){
        refDF$case <- paste(refDF$lin, 
                            refDF$host, 
                            sep = "_")
        bestSamples <- c()
        for(case in unique(refDF$case)){
                caseSamps <- refDF[refDF$case == case, ]
                if(nrow(caseSamps) > n){
                        bestSamps <- caseSamps$oligo[order(caseSamps$Ct)][1:n]
                }else{
                        bestSamps <- caseSamps$oligo
                }
                bestSamples <- c(bestSamples, bestSamps)
        }
        return(bestSamples)
}

# Do boxplots.
doExpBoxPlot <- function(data, plateNorm = T){
        data$lin <- factor(data$lin, levels = c("Bovis", "Chimp", "L6", "L5"))
        data$host <- factor(data$host, levels = c("BOMAC", "THP-1"))
        if(plateNorm == FALSE){
                boxPlot <- ggplot(data, aes(x = lin, fill = lin, y = relExp))
        }else{
                boxPlot <- ggplot(data, aes(x = lin, fill = lin, y = relExp_plateNorm))
        }
        boxPlot <- boxPlot +
        geom_boxplot_pattern(pattern_color = "black",
                             pattern_fill = "black",
                             pattern_spacing = 0.015,
                             aes(pattern = host)) + 
                scale_pattern_manual(values = c("none", "stripe"),
                                     labels = c("BOMAC", "THP-1")) +
                labs(title = unique(data$primerPair),
                     x = "Strain",
                     y = "Expression (relative to *sigA*)") +
                #scale_fill_manual(values = c("#c4bf62", "#24ad37", "#871414"), 
                #                  labels = c("Chimp", "L6", "L5")) +
                scale_fill_manual(values = c("#f279ce", "#c4bf62", "#24ad37", "#871414"), 
                                  labels = c("Bovis", "Chimp", "L6", "L5")) +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"))
        return(boxPlot)
}

doBoxPlots <- function(parsRes, plateNorm = FALSE){
        parsRes$lin <- factor(parsRes$lin, levels = c("Bovis", "Chimp", "L6", "L5"))
        parsRes$host <- factor(parsRes$host, levels = c("BOMAC", "THP-1"))
        #plcA_res <- remBovis(parsRes[parsRes$primerPair == "Rv2351c", ])
        plcA_res <- parsRes[parsRes$primerPair == "Rv2351c", ]
        #plcB_res <- remBovis(parsRes[parsRes$primerPair == "Rv2350c", ])
        plcB_res <- parsRes[parsRes$primerPair == "Rv2350c", ]
        #plcC_res <- remBovis(parsRes[parsRes$primerPair == "Rv2349c", ])
        plcC_res <- parsRes[parsRes$primerPair == "Rv2349c", ]
        plcD_res <- parsRes[parsRes$primerPair == "Mb1784c", ]
        if(plateNorm == FALSE){
                plcAPlot <- ggplot(plcA_res, aes(x = lin, fill = lin, y = relExp))
                plcBPlot <- ggplot(plcB_res, aes(x = lin, fill = lin, y = relExp))
                plcCPlot <- ggplot(plcC_res, aes(x = lin, fill = lin, y = relExp))
                plcDPlot <- ggplot(plcD_res, aes(x = lin, fill = lin, y = relExp))
        }else{
                plcAPlot <- ggplot(plcA_res, aes(x = lin, fill = lin, y = relExp_plateNorm))
                plcBPlot <- ggplot(plcB_res, aes(x = lin, fill = lin, y = relExp_plateNorm))
                plcCPlot <- ggplot(plcC_res, aes(x = lin, fill = lin, y = relExp_plateNorm))
                plcDPlot <- ggplot(plcD_res, aes(x = lin, fill = lin, y = relExp_plateNorm))
        }
        plcAPlot <- plcAPlot +           # Create boxplot with pattern
                geom_boxplot_pattern(pattern_color = "black",
                                     pattern_fill = "black",
                                     pattern_spacing = 0.015,
                                     aes(pattern = host)) + 
                scale_pattern_manual(values = c("none", "stripe"),
                                     labels = c("BOMAC", "THP-1")) +
                labs(title = "plcA (Rv2351c)",
                     x = "Strain",
                     y = "Expression (relative to *sigA*)") +
                #scale_fill_manual(values = c("#c4bf62", "#24ad37", "#871414"), 
                #                  labels = c("Chimp", "L6", "L5")) +
                scale_fill_manual(values = c("#f279ce", "#c4bf62", "#24ad37", "#871414"), 
                                  labels = c("Bovis", "Chimp", "L6", "L5")) +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"))
        
        plcBPlot <- plcBPlot +           # Create boxplot with pattern
                geom_boxplot_pattern(pattern_color = "black",
                                     pattern_fill = "black",
                                     pattern_spacing = 0.015,
                                     aes(pattern = host)) + 
                scale_pattern_manual(values = c("none", "stripe"),
                                     labels = c("BOMAC", "THP-1")) +
                labs(title = "plcB (Rv2350c)",
                     x = "Strain",
                     y = "Expression (relative to *sigA*)") +
                #scale_fill_manual(values = c("#c4bf62", "#24ad37", "#871414"), 
                #                  labels = c("Chimp", "L6", "L5")) +
                scale_fill_manual(values = c("#f279ce", "#c4bf62", "#24ad37", "#871414"), 
                                  labels = c("Bovis", "Chimp", "L6", "L5")) +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"))
        
        plcCPlot <- plcCPlot +           # Create boxplot with pattern
                geom_boxplot_pattern(pattern_color = "black",
                                     pattern_fill = "black",
                                     pattern_spacing = 0.015,
                                     aes(pattern = host)) + 
                scale_pattern_manual(values = c("none", "stripe"),
                                     labels = c("BOMAC", "THP-1")) +
                labs(title = "plcC (Rv2349c)",
                     x = "Strain",
                     y = "Expression (relative to *sigA*)") +
                #scale_fill_manual(values = c("#c4bf62", "#24ad37", "#871414"), 
                #                  labels = c("Chimp", "L6", "L5")) +
                scale_fill_manual(values = c("#f279ce", "#c4bf62", "#24ad37", "#871414"), 
                                  labels = c("Bovis", "Chimp", "L6", "L5")) +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"))
        
        plcDPlot <- plcDPlot +           # Create boxplot with pattern
                geom_boxplot_pattern(pattern_color = "black",
                                     pattern_fill = "black",
                                     pattern_spacing = 0.015,
                                     aes(pattern = host)) + 
                scale_pattern_manual(values = c("none", "stripe"),
                                     labels = c("BOMAC", "THP-1")) +
                labs(title = "plcD (Mb1784c)",
                     x = "Strain",
                     y = "Expression (relative to *sigA*)") +
                scale_fill_manual(values = c("#f279ce", "#c4bf62", "#24ad37", "#871414"), 
                                  labels = c("Bovis", "Chimp", "L6", "L5")) +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"))
        
        allPlot <- ggarrange(plcAPlot + 
                                     rremove("xlab"), 
                             plcBPlot + 
                                     rremove("ylab") + 
                                     rremove("xlab"), 
                             plcCPlot, 
                             plcDPlot + 
                                     rremove("ylab"),
                             legend = "bottom",
                             common.legend = T)
        out <- list(sepGeneExp = list(plcA = plcA_res,
                                      plcB = plcB_res,
                                      plcC = plcC_res,
                                      plcD = plcD_res),
                    plot = allPlot)
        return(out)
}

# Do statistical tests between the two types of host cells. 
doShap <- function(gene, lin, host, which = "all", plateNorm = FALSE){
        if(which == "all"){
                allDat <- plotAll$sepGeneExp
        }else if(which == "new"){
                allDat <- plotNew$sepGeneExp
        }
        if(gene == "plcA"){
                data <- allDat$plcA
        }else if(gene == "plcB"){
                data <- allDat$plcB
        }else if(gene == "plcC"){
                data <- allDat$plcC
        }else if(gene == "plcD"){
                data <- allDat$plcD
        }
        if(plateNorm == FALSE){
                shap <- shapiro.test(data$relExp[data$lin == lin & data$host == host])
        }else{
                shap <- shapiro.test(data$relExp_plateNorm[data$lin == lin & data$host == host])
        }
        return(shap)
}

doMW <- function(gene, lin, which = "all", plateNorm = FALSE){
        if(which == "all"){
                allDat <- plotAll$sepGeneExp
        }else if(which == "new"){
                allDat <- plotNew$sepGeneExp
        }
        if(gene == "plcA"){
                data <- allDat$plcA
        }else if(gene == "plcB"){
                data <- allDat$plcB
        }else if(gene == "plcC"){
                data <- allDat$plcC
        }else if(gene == "plcD"){
                data <- allDat$plcD
        }
        if(plateNorm == FALSE){
                bomac <- data$relExp[data$lin == lin & data$host == "BOMAC"]
                thp1 <- data$relExp[data$lin == lin & data$host == "THP-1"]
        }else{
                bomac <- data$relExp_plateNorm[data$lin == lin & data$host == "BOMAC"]
                thp1 <- data$relExp_plateNorm[data$lin == lin & data$host == "THP-1"]
        }
        print(sprintf("BOMAC expression: %s", paste(bomac, collapse = ", ")))
        print(sprintf("THP-1 expression: %s", paste(thp1, collapse = ", ")))
        MW <- wilcox.test(bomac,
                          thp1)
        return(MW)
}

doTTst <- function(gene, lin, which = "all", plateNorm = FALSE){
        if(which == "all"){
                allDat <- plotAll$sepGeneExp
        }else if(which == "new"){
                allDat <- plotNew$sepGeneExp
        }
        if(gene == "plcA"){
                data <- allDat$plcA
        }else if(gene == "plcB"){
                data <- allDat$plcB
        }else if(gene == "plcC"){
                data <- allDat$plcC
        }else if(gene == "plcD"){
                data <- allDat$plcD
        }
        if(plateNorm == FALSE){
                bomac <- data$relExp[data$lin == lin & data$host == "BOMAC"]
                thp1 <- data$relExp[data$lin == lin & data$host == "THP-1"]
        }else{
                bomac <- data$relExp_plateNorm[data$lin == lin & data$host == "BOMAC"]
                thp1 <- data$relExp_plateNorm[data$lin == lin & data$host == "THP-1"]
        }
        print(sprintf("BOMAC expression: %s", paste(bomac, collapse = ", ")))
        print(sprintf("THP-1 expression: %s", paste(thp1, collapse = ", ")))
        tTst <- t.test(bomac,
                       thp1)
        return(tTst)
}

# Obtain full R number (RNA codes), with the full number of zeros in between
getFullRNums <- function(rNumInt){
        if(rNumInt != "NO_RNA"){
                nums <- sprintf("%0.3d", as.numeric(gsub("R0", 
                                                         "", 
                                                         strsplit(rNumInt, 
                                                                  " - ")[[1]]))[1]:as.numeric(gsub("R0", 
                                                                                                   "", 
                                                                                                   strsplit(rNumInt, 
                                                                                                            " - ")[[1]]))[2])
                strNums <- c()
                for(i in seq_along(nums)){
                        num <- paste0("R", 
                                      paste0(rep("0", (4 - nchar(nums[i]))), 
                                             collapse = ""), 
                                      as.character(nums[i]))
                        strNums <-  c(strNums, num)
                }
                strNums <- paste0(strNums, collapse = ", ")
        }else{
                strNums <- rNumInt
        }
        return(strNums)
}

getSigACPosCts <- function(parsedExp){
        controlsLog <- grepl("C+", parsedExp$Sample.Name, fixed = T)
        if(sum(controlsLog) == 0){
                controlsLog <- grepl("3e+6", parsedExp$Sample.Name, fixed = T)
        }
        controls <- parsedExp[controlsLog, ]
        controls <- controls[grepl("sigA", controls$Sample.Name), ]
        sigACt <- mean(as.numeric(controls$CT))
        return(sigACt)
}


# Filter citometry data to keep just samples that have data for the four types of cell death,
# and add columns stating number of replicates, cases and samples ids
keepCompSamps <- function(cito){
        compSmpIdxs <- c()
        for(i in 1:(nrow(cito) - 3)){
                idxWndw <- i:(i + 3)
                #print(idxWndw)
                #print(cito$ID_cito[idxWndw])
                cito$Muerte_cito[idxWndw]
                
                #print(idxWndw %in% compSmpIdxs)
                datesInf <- cito$Fecha_INF[idxWndw]
                uniqDates <- length(unique(datesInf))
                #print(datesInf)
                #if(cito$Cepas[i] == "L5" & cito$Célula[i] == "BoMac"){
                
                #}
                if(all(cito$Muerte_cito[idxWndw] == unique(cito$Muerte_cito)) & (sum(idxWndw %in% compSmpIdxs) == 0) & uniqDates == 1){
                        #print(idxWndw %in% compSmpIdxs)
                        #print(datesInf)
                        compSmpIdxs <- c(compSmpIdxs, idxWndw)
                }
        }
        cito <- cito[compSmpIdxs, ]
        cito <- cito[order(paste0(cito$Cepas, cito$Célula)), ]
        rownames(cito) <- 1:nrow(cito)
        cases <- unique(paste(cito$Cepas, cito$Célula, sep = "_"))
        #print(cases)
        repVec <- c()
        for(i in seq_along(cases)){
                case <- cases[i]
                caseDeaths <- cito$Muerte_cito[paste(cito$Cepas, cito$Célula, sep = "_") == case]
                nReps <- length(caseDeaths)/4
                repIDs <- rep(paste0("R", 1:nReps), each = 4)
                repVec <- c(repVec, repIDs)
        }
        #print(repVec)
        #View(cito)
        cito$Rep <- repVec
        
        cito$Case <- paste(cito$Cepas, cito$Célula, sep = "_")
        cito$Case <- gsub("BoMac", "BOMAC", cito$Case)
        cito$Case <- gsub("THP1", "THP-1", cito$Case)
        cito$SampID <- paste(cito$Case, cito$Rep, sep = "_")
        cito$Valor_cito <- as.numeric(cito$Valor_cito)
        return(cito)
}

# Create a numeric matrix where each row is a sample and each column a type of death
createNumMat <- function(cito, useFCs = F){
        colnames(cito) <- gsub("sample", "SampID", colnames(cito))
        colnames(cito) <- gsub("valor", "Valor_cito", colnames(cito))
        colnames(cito) <- gsub("tipo_muerte", "Muerte_cito", colnames(cito))
        samps <- unique(cito$SampID)
        death <- sort(unique(cito$Muerte_cito))
        outDF <- data.frame(matrix(nrow = 0,
                                   ncol = length(death),
                                   dimnames = list(NULL, 
                                                   death)))
        for(samp in samps){
                if(useFCs){
                        valVec <- cito$FC_ctrl[cito$SampID == samp]
                }else{
                        valVec <- cito$Valor_cito[cito$SampID == samp]
                }
                deathVec <- cito$Muerte_cito[cito$SampID == samp]
                valDF <- data.frame(matrix(valVec,
                                           nrow = 1, 
                                           ncol = length(death),
                                           dimnames = list(samp, deathVec)))
                valDF <- valDF[, order(colnames(valDF))]
                outDF <- rbind.data.frame(outDF, valDF)
        }
        return(outDF)
}

getExpMat <- function(resDF, plateNorm = F, genes = c("Rv2349c", "Rv2350c", "Rv2351c", "Mb1784c")){
        plcGenes <- genes
        resDF <- resDF[resDF$primerPair %in% plcGenes, ]
        resDF$case <- paste(resDF$lin, resDF$host, sep = "_")
        cases <- unique(resDF$case)
        resDFOrd <- data.frame(matrix(nrow = 0, 
                                      ncol = ncol(resDF),
                                      dimnames = list(NULL,
                                                      colnames(resDF))))
        for(case in cases){
                resSubCases <- resDF[resDF$case == case, ]
                resSubCases <- resSubCases[order(resSubCases$oligo), ]
                resDFOrd <- rbind.data.frame(resDFOrd, resSubCases)
        }
        resDFOrd$sampID <- paste(resDFOrd$case, resDFOrd$oligo, sep = "_")
        sampIDs <- unique(resDFOrd$sampID)
        # Filter sample IDs to keep just the ones that have data of all the genes (this is necessary when including the new cholesterol and mbt genes, 
        # as some samples don't have data)
        sampIDs <- table(resDFOrd$sampID) == max(table(resDFOrd$sampID)) # Edit
        sampIDs <- names(sampIDs)[sampIDs]
        expMat <- data.frame(matrix(nrow = 0, 
                                    ncol = length(plcGenes), 
                                    dimnames = list(NULL,
                                                    sort(plcGenes))))
        #print(expMat)
        for(i in seq_along(sampIDs)){
                sampID <- sampIDs[i]
                if(plateNorm == F){
                        sampRelExps <- resDFOrd$relExp[resDFOrd$sampID == sampID]
                }else if(plateNorm == T){
                        sampRelExps <- resDFOrd$relExp_plateNorm[resDFOrd$sampID == sampID]
                }
                #print(sampRelExps)
                sampGenes <- resDFOrd$primerPair[resDFOrd$sampID == sampID]
                names(sampRelExps) <- sampGenes
                #print(sampGenes)
                notInPlc <- plcGenes[!plcGenes %in% sampGenes]
                if(length(notInPlc) > 0){
                        notInPlcExps <- rep(0, length(notInPlc))
                        names(notInPlcExps) <- notInPlc
                        sampRelExps <- c(sampRelExps, notInPlcExps)
                }
                #print("JA")
                #print(notInPlc)
                sampRelExps <- data.frame(matrix(sampRelExps[order(names(sampRelExps))],
                                                 nrow = 1,
                                                 dimnames = list(sampID, sort(names(sampRelExps)))))
                sampRelExps <- sampRelExps[, colnames(sampRelExps) %in% plcGenes]
                #print(sampRelExps)
                expMat <- rbind.data.frame(expMat, 
                                           sampRelExps)
                #print(expMat)
        }
        return(expMat)
}

getExpCitMat <- function(citoMat, expMat, seed = 777){
        expCitMat <- data.frame(matrix(nrow = 0, ncol = 8, dimnames = list(NULL, 
                                                                           c(colnames(expMat),
                                                                             sort(unique(citoMat$Muerte_cito))))))
        for(i in seq_along(unique(citoMat$Case))){
                case <- unique(citoMat$Case)[i]
                citSubMat <- citoMat[citoMat$Case == case, ]
                expSubMat <- expMat[stringr::str_extract(rownames(expMat), "[^_]*_[^_]*") == case, ]
                numToSamp <- nrow(expSubMat)
                reps <- unique(citSubMat$Rep)
                if(numToSamp > length(reps)){
                        #set.seed(123)
                        #set.seed(321)
                        #set.seed(666)
                        #set.seed(777) # The used one # Default in function. Added seed argument to be able to do several rounds to identify enriched
                        #set.seed(555) # More significatives when running separate cell types
                        set.seed(seed)
                        sampledReps <- sample(reps, numToSamp, replace = T)
                }else{
                        #set.seed(123)
                        #set.seed(321)
                        #set.seed(666)
                        #set.seed(777) # The used one # Default in function. Added seed argument to be able to do several rounds to identify enriched
                        #set.seed(555) # More significatives when running separate cell types
                        set.seed(seed)
                        sampledReps <- sample(reps, numToSamp, replace = F)
                }
                citoSampMat <- data.frame(matrix(nrow = 0, 
                                                 ncol = 4, 
                                                 dimnames = list(NULL, 
                                                                 sort(unique(citoMat$Muerte_cito)))))
                for(j in seq_along(sampledReps)){
                        rep <- sampledReps[j]
                        deathType <- citSubMat$Muerte_cito[citSubMat$Rep == rep]
                        deathValue <- citSubMat$Valor_cito[citSubMat$Rep == rep]
                        deathValue <- as.numeric(as.character(deathValue))
                        names(deathValue) <- deathType
                        deathValue <- deathValue[order(names(deathValue))]
                        deathValue <- data.frame(matrix(deathValue, 
                                                        nrow = 1, 
                                                        dimnames = list(rep, 
                                                                        names(deathValue))))
                        citoSampMat <- rbind.data.frame(citoSampMat,
                                                        deathValue)
                }
                expCitSubMat <- cbind.data.frame(expSubMat, citoSampMat)
                expCitMat <- rbind.data.frame(expCitMat, expCitSubMat)
        }
        return(expCitMat)
}

# Do the biplot of the mixed expression and citometry data. 
pcBiplot <- function(PC, x="PC1", y="PC2", varPlotFilt = NULL, formatGeneNames = T, labels = F){
        if(!require(ggrepel)) install.packages(ggrepel)
        library(ggrepel)
        # Generate dataframe of gene names and locus tag equivalence
        geneNames <- data.frame(locus = c("Rv2349c",
                                          "Rv2350c",
                                          "Rv2351c",
                                          "Mb1784c"),
                                gene = c("plcC",
                                         "plcB",
                                         "plcA",
                                         "plcD"),
                                stringsAsFactors = F)
        bactCols <- data.frame(Bacteria = c("L5", "L6", "Bovis", "Chimp", "ctrl"),
                               legendName = c("L5", "L6", "*M. bovis*", "Chimpanzee Bacillus", "Uninfected"),
                               color = c("#871414", "#24ad37", "#f279ce", "#c4bf62", "#808080"),
                               stringsAsFactors = F)
        data <- data.frame(obsnames=row.names(PC$x), PC$x)
        
        data <- data[, c("obsnames", x, y)]
        data$Bacteria <- sapply(data$obsnames, function(x) strsplit(as.character(x), "_")[[1]][1])
        data$Host <- sapply(data$obsnames, function(x) strsplit(as.character(x), "_")[[1]][2])
        data$Host <- gsub("BOMAC|bomac", "BoMac", data$Host)
        data$Host <- gsub("thp1", "THP-1", data$Host)
        data$Host <- factor(data$Host, levels = c("THP-1", "BoMac"))
        data$Bacteria <- gsub("bovis", "Bovis", data$Bacteria)
        data$Bacteria <- gsub("chimp", "Chimp", data$Bacteria)
        data$Bacteria <- gsub("control", "ctrl", data$Bacteria)
        bactLevls <- c("L5", "L6", "Bovis", "Chimp", "ctrl")[c("L5", "L6", "Bovis", "Chimp", "ctrl") %in% data$Bacteria]
        data$Bacteria <- factor(data$Bacteria, levels = bactLevls)
        propVar <- summary(PC)$importance[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)
        plotCols <- bactCols$color[match(levels(data$Bacteria), bactCols$Bacteria)]
        if(all(is.na(data$Bacteria))){
                plot <- ggplot(data, aes(data[, 2], 
                                         data[, 3], 
                                         shape = Host))
        }else{
                plot <- ggplot(data, aes(data[, 2], 
                                         data[, 3], 
                                         col = Bacteria,
                                         shape = Host)) + 
                        scale_discrete_manual("Bacteria",
                                              aesthetics = "colour",
                                              values = plotCols,
                                              labels = bactCols$legendName[match(levels(data$Bacteria),
                                                                                 bactCols$Bacteria)])
        }
        plot <- plot +
                geom_hline(yintercept = 0, alpha = 0.6) +
                geom_vline(xintercept = 0, alpha = 0.6) +
                geom_point() +
                xlab(sprintf("%s (%s%%)", x, propX)) +
                ylab(sprintf("%s (%s%%)", y, propY)) + 
                theme(axis.text.y = element_text(size=10),
                      axis.text.x = element_text(size=10),
                      axis.title.x = element_text(size=10),
                      axis.title.y = element_text(size=10),
                      legend.text = ggtext::element_markdown(size = 8),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black", fill=NA, size=1))
        if(labels == T){
                plot <- plot +
                        geom_text_repel(aes(x = data[, 2], 
                                            y = data[, 3],
                                            col = Bacteria,
                                            label = obsnames),
                                        data = data)
        }
        datapc <- data.frame(varnames=rownames(PC$rotation), 
                             PC$rotation,
                             stringsAsFactors = F)
        mult <- min(
                (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
                (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
        )
        datapc <- transform(datapc,
                            v1 = .7 * mult * (get(x)),
                            v2 = .7 * mult * (get(y))
        )
        datapc$x0 <- rep(0, nrow(datapc))
        datapc$y0 <- rep(0, nrow(datapc))
        if(!is.null(varPlotFilt)){
                datapc <- datapc[datapc$varnames %in% varPlotFilt, ]
        }
        if(formatGeneNames == T){
                outGeneNames <- geneNames$gene[match(datapc$varnames[datapc$varnames %in% geneNames$locus],
                                                     geneNames$locus)]
                outLabs <- sprintf("%s (%s)", outGeneNames, datapc$varnames[datapc$varnames %in% geneNames$locus])
                datapc$varnames[datapc$varnames %in% geneNames$locus] <- outLabs
                datapc$varnames <- gsub("_", " ", datapc$varnames)
        }
        plot <- plot +
                geom_text_repel(data=datapc, 
                                aes(x=v1, y=v2, label=varnames), 
                                color = "black", 
                                size = 3,
                                inherit.aes = F) + 
                geom_segment(data = datapc, aes(x=x0, 
                                                y=y0, 
                                                xend=v1, 
                                                yend=v2),
                             arrow = arrow(length=unit(0.2,"cm"),
                                           type = "closed",
                                           angle = 20), 
                             alpha=0.75, 
                             color="black", 
                             size = 0.5,
                             inherit.aes = F)
        plot
}

pcBiplot_date <- function(PC, x = "PC1", y = "PC2", sampInfo, linNames = F){
        data <- data.frame(obsnames=row.names(PC$x), PC$x)
        
        data <- data[, c("obsnames", x, y)]
        data$date <- as.factor(sampInfo$Infeccion_fecha[match(data$obsnames, sampInfo$sample)])
        data$Host <- sapply(data$obsnames, function(x) strsplit(as.character(x), "_")[[1]][2])
        data$Host <- gsub("BOMAC|bomac", "BoMac", data$Host)
        data$Host <- gsub("thp1", "THP-1", data$Host)
        data$Host <- factor(data$Host, levels = c("THP-1", "BoMac"))
        data$Bacteria <- sapply(data$obsnames, function(x) strsplit(as.character(x), "_")[[1]][1])
        data$Bacteria <- gsub("bovis", "Bovis", data$Bacteria)
        data$Bacteria <- gsub("chimp", "Chimp", data$Bacteria)
        data$Bacteria <- gsub("control", "ctrl", data$Bacteria)
        bactLevls <- c("L5", "L6", "Bovis", "Chimp", "ctrl")[c("L5", "L6", "Bovis", "Chimp", "ctrl") %in% data$Bacteria]
        data$Bacteria <- factor(data$Bacteria, levels = bactLevls)
        propVar <- summary(PC)$importance[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)
        plot <- ggplot(data, aes(data[, 2], 
                                 data[, 3], 
                                 col = date,
                                 shape = Host)) + 
                geom_hline(yintercept = 0, alpha = 0.6) +
                geom_vline(xintercept = 0, alpha = 0.6) +
                geom_point() +
                xlab(sprintf("%s (%s%%)", x, propX)) +
                ylab(sprintf("%s (%s%%)", y, propY)) + 
                theme(axis.text.y = element_text(size=10),
                      axis.text.x = element_text(size=10),
                      axis.title.x = element_text(size=10),
                      axis.title.y = element_text(size=10),
                      legend.text = ggtext::element_markdown(size = 8),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black", fill=NA, size=1))
        if(linNames == T){
                plot <- plot +
                        geom_text_repel(aes(x = data[, 2], 
                                            y = data[, 3],
                                            col = date,
                                            label = Bacteria),
                                        data = data)
        }
        
        
        datapc <- data.frame(varnames=rownames(PC$rotation), 
                             PC$rotation,
                             stringsAsFactors = F)
        mult <- min(
                (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
                (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
        )
        datapc <- transform(datapc,
                            v1 = .7 * mult * (get(x)),
                            v2 = .7 * mult * (get(y))
        )
        datapc$x0 <- rep(0, nrow(datapc))
        datapc$y0 <- rep(0, nrow(datapc))
        plot <- plot +
                geom_text_repel(data=datapc, 
                                aes(x=v1, y=v2, label=varnames), 
                                color = "black", 
                                size = 3,
                                inherit.aes = F) + 
                geom_segment(data = datapc, aes(x=x0, 
                                                y=y0, 
                                                xend=v1, 
                                                yend=v2),
                             arrow = arrow(length=unit(0.2,"cm"),
                                           type = "closed",
                                           angle = 20), 
                             alpha=0.75, 
                             color="black", 
                             size = 0.5,
                             inherit.aes = F)
        plot
}

# Do screePlot of a PCA object
screePlot <- function(pca, nComps = 4){
        screeDF <- data.frame(PC = paste0("PC", as.character(1:length(pca$sdev))),
                              varExp = pca$sdev^2/sum(pca$sdev^2))
        screeDF <- screeDF[1:nComps, ]
        plot <-ggplot(data = screeDF, aes(x = PC, y = varExp, group = 1)) + 
                geom_bar(stat = "identity") + 
                geom_line() +
                geom_point() +
                # geom_line(data = screeDF, aes(x = PC, y = varExp)) + 
                xlab("Principal Component") + 
                ylab("Variance Explained") +
                ggtitle("Scree Plot") +
                ylim(0, 1) + 
                theme(axis.text.y = element_text(size=10),
                      axis.text.x = element_text(size=10),
                      axis.title.x = element_text(size=10),
                      axis.title.y = element_text(size=10),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black", fill=NA, size=1))
        plot
}

# Do regression plot of bidimensional model
doRegPlot <- function(model, varY, varX){
        bactCols <- data.frame(Bacteria = c("L5", "L6", "Bovis", "Chimp"),
                               color = c("#871414", "#24ad37", "#f279ce", "#c4bf62"),
                               stringsAsFactors = F)
        geneNames <- data.frame(locus = c("Rv2349c",
                                          "Rv2350c",
                                          "Rv2351c",
                                          "Mb1784c"),
                                gene = c("plcC",
                                         "plcB",
                                         "plcA",
                                         "plcD"))
        shapeVec <- c(16, 17)
        names(shapeVec) <- c("THP-1", "BOMAC")
        modelDF <- model$model
        modelDF$Bacteria <- sapply(rownames(modelDF), function(x) strsplit(x, "_")[[1]][1])
        print(modelDF$Bacteria)
        modelDF$Host <- sapply(rownames(modelDF), function(x) strsplit(x, "_")[[1]][2])
        modelDF$Bacteria <- factor(modelDF$Bacteria, 
                                   levels = unique(modelDF$Bacteria))
        modelDF$Host <- factor(modelDF$Host,
                               levels = c("THP-1", "BOMAC"))
        plotCols <- bactCols$color[match(levels(modelDF$Bacteria), bactCols$Bacteria)]
        #print(resFit$Bacteria)
        print(plotCols)
        plotShapes <- shapeVec[match(levels(modelDF$Host),
                                     names(shapeVec))]
        plot <- ggplot(data = modelDF, aes_string(x = varX, y = varY)) +
                geom_smooth(size = 1, method = "lm") +
                geom_point(aes_string(y = varY, col = "Bacteria", shape = "Host")) + 
                scale_discrete_manual("Bacteria",
                                      aesthetics = "colour",
                                      values = plotCols) +
                scale_shape_manual("Host", values = plotShapes) +
                #geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
                # theme_bw(base_size = 14) +
                theme(axis.title.x = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      plot.caption=element_text(size = 8.5,
                                                #hjust = 0,
                                                margin = margin(5,0,0,0)))# +
        #labs(x = xAxisTitle,
        #     y = yAxisTitle,
        #     caption = statStr)
        plot
}


# This function does added variable plots. Added variable plots statistical basis is partial regression, or Frisch-Waugh-Lowell 
# theorem, which establishes the individual residuals of each independent variable of a multiple regression can be obtained
# by doing a regression of the residuals of Y regressed on the independent variables with the exception of the objective 
# independent variable and the residuals of the objective independent variable on the rest of independent variables. 
# Therefore it fits a linear model on the variables in the global model without the objective variable, with the type of 
# death as response variable by one side, and by other side a model on the variables in the global model without the objective 
# variable with the objective variable as response variable. We obtain the residuals of both models. The residuals reflect the 
# amount of the response variable (either type of death or objective variable) not explained by the variables other than the 
# objective variable. The strongest the relationship between the two sets of residuals it the more improves the addition the 
# objective to the global prediction. The slope of the line corresponds to the slope of the objective variable in the global 
# model. (Gallup et al 2019, doi:10.1177/1536867X19874236)
doAddVarPlots1 <- function(model, variable, formatGeneNames = T){
        if(!require(broom)) install.packages("broom")
        library(broom)
        bactCols <- data.frame(Bacteria = c("L5", "L6", "Bovis", "Chimp"),
                               color = c("#871414", "#24ad37", "#f279ce", "#c4bf62"),
                               stringsAsFactors = F)
        geneNames <- data.frame(locus = c("Rv2349c",
                                          "Rv2350c",
                                          "Rv2351c",
                                          "Mb1784c"),
                                gene = c("plcC",
                                         "plcB",
                                         "plcA",
                                         "plcD"))
        shapeVec <- c(16, 17)
        names(shapeVec) <- c("THP-1", "BOMAC")
        modMat <- model.matrix(model)
        respVar <- all.vars(formula(model))[1]
        modVars <- all.vars(formula(model))[-1]
        varIdx <- which(colnames(modMat) == variable)
        respVec <- model$model[, colnames(model$model) == respVar]
        respFit <- lm(formula(paste0(respVar,
                                     " ~ ",
                                     paste(modVars[modVars != variable],
                                           collapse = " + "))),
                      data = model$model)
        varFit <- lm(formula(paste0(variable,
                                    " ~ ",
                                    paste(modVars[modVars != variable],
                                          collapse = " + "))),
                     data = model$model)
        res <- data.frame(variable = varFit$residuals,
                          resp = respFit$residuals)
        resFit <- lm(formula("resp ~ variable"), res)
        resFit <- augment(resFit, newdata = res, se_fit = T)
        resFit <- mutate(resFit, 
                         lower = .fitted - 2*.se.fit,
                         upper = .fitted + 2*.se.fit,
                         pred = .fitted)
        resFit <- as.data.frame(resFit)
        resFit$Bacteria <- sapply(resFit$.rownames, 
                                  function(x) strsplit(x, split = "_")[[1]][1])
        bactFactLevels <- c("L5", "L6", "Bovis", "Chimp")[c("L5", "L6", "Bovis", "Chimp") %in% resFit$Bacteria]
        #print(bactFactLevels)
        resFit$Bacteria <- factor(resFit$Bacteria, 
                                  levels = unique(resFit$Bacteria))
        resFit$Host <- factor(sapply(resFit$.rownames, 
                                     function(x) strsplit(x, split = "_")[[1]][2]),
                              levels = c("THP-1", "BOMAC"))
        yAxisTitle <- paste0(toupper(substr(gsub("_", " ", respVar), 1, 1)), 
                             substr(gsub("_", " ", respVar), 2, nchar(respVar)),
                             " | others")
        if(formatGeneNames == T){
                xAxisTitle = sprintf("*%s* (%s) | others", 
                                     geneNames$gene[match(variable, geneNames$locus)],
                                     variable)
        }else{
                xAxisTitle = paste(variable, "others", sep = " | ")
        }
        modSum <- summary(model)[[4]]
        coef <- model$coefficients[names(model$coefficients) == variable]
        coef <- as.character(round(coef, digits = 3))
        tStat <- modSum[rownames(modSum) == variable, 3]
        tStat <- as.character(round(tStat, digits = 3))
        pVal <- modSum[rownames(modSum) == variable, 4]
        pVal <- as.character(round(pVal, digits = 3))
        statStr <- sprintf("\u03B2 = %s          t-statistic = %s          p-value = %s",
                           coef,
                           tStat,
                           pVal)
        plotCols <- bactCols$color[match(levels(resFit$Bacteria), bactCols$Bacteria)]
        #print(resFit$Bacteria)
        #print(plotCols)
        plotShapes <- shapeVec[match(levels(resFit$Host),
                                     names(shapeVec))]
        plot <- ggplot(data = resFit, aes(x = variable, y = pred)) +
                geom_line(size = 1) +
                geom_point(aes(y = resp, col = Bacteria, shape = Host)) + 
                scale_discrete_manual("Bacteria",
                                      aesthetics = "colour",
                                      values = plotCols) +
                scale_shape_manual("Host", values = plotShapes) +
                geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
                # theme_bw(base_size = 14) +
                theme(axis.title.x = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      plot.caption=element_text(size = 8.5,
                                                #hjust = 0,
                                                margin = margin(5,0,0,0))) +
                labs(x = xAxisTitle,
                     y = yAxisTitle,
                     caption = statStr)
        plot
}


doMergedAddVarPlots <- function(model, formatGeneNames = T){
        variables <- all.vars(formula(model))
        plotList <- list()
        for(v in variables){
                p <- doAddVarPlots1(model, v, formatGeneNames = formatGeneNames)
                plotList[[v]] <- p
        }
        allPlots <- ggarrange(plotList[["Rv2351c"]],
                              plotList[["Rv2350c"]] + 
                                      rremove("ylab"), 
                              plotList[["Rv2349c"]], 
                              plotList[["Mb1784c"]] + 
                                      rremove("ylab"),
                              legend = "bottom",
                              common.legend = T)
        modSum <- summary(model)
        rSquared <- as.character(round(modSum[8][[1]], 
                                       digits = 3))
        adjRSqrd <- as.character(round(modSum[9][[1]],
                                       digits = 3))
        f <- modSum[10][[1]]
        fStat <- as.character(round(f[1],
                                    digits = 3))
        pVal <- as.character(round(pf(f[1], f[2], f[3], lower.tail = F),
                                   digits = 3))
        statStr <- sprintf("R^2 = %s          adjusted R^2 = %s          F-statistic = %s          p-value = %s",
                           rSquared,
                           adjRSqrd,
                           fStat,
                           pVal)
        allPlots <- annotate_figure(allPlots,
                                    bottom = text_grob(statStr,
                                                       size = 10))
        allPlots
}

# This function does added variable plots by obtaining predictions based on the global model
# about how the response (type of death) is gonna behave if the variables other than the 
# objective are kept constant to the value of the median. 
doAddVarPlots2 <- function(model, data){
        if(!require(dplyr)) install.packages("dplyr")
        library(dplyr)
        if(!require(purrr)) install.packages("purrr")
        library(purrr)
        if(!require(broom)) install.packages("broom")
        library(broom)
        bactCols <- data.frame(Bacteria = c("L5", "L6", "Bovis", "Chimp"),
                               color = c("#871414", "#24ad37", "#f279ce", "#c4bf62"),
                               stringsAsFactors = F)
        geneNames <- data.frame(locus = c("Rv2349c",
                                          "Rv2350c",
                                          "Rv2351c",
                                          "Mb1784c"),
                                gene = c("plcC",
                                         "plcB",
                                         "plcA",
                                         "plcD"))
        preddat_fun = function(data, allvars, var){
                sums = summarise_at(data, 
                                    vars( one_of(allvars), -one_of(var) ), 
                                    median) 
                cbind( select_at(data, var), sums)
        }
        respVar <- all.vars(formula(model))[1]
        mod_vars <- all.vars(formula(model))[-1]
        pred_dats = mod_vars %>%
                set_names() %>%
                map(~preddat_fun(data, mod_vars, .x))
        preds = pred_dats %>%
                map(~augment(model, newdata = .x, se_fit = TRUE) ) %>%
                map(~mutate(.x, 
                            lower = .fitted - 2*.se.fit,
                            upper = .fitted + 2*.se.fit,
                            pred = .fitted))
        yAxisTitle <- paste0(toupper(substr(gsub("_", " ", respVar), 1, 1)), 
                             substr(gsub("_", " ", respVar), 2, nchar(respVar)))
        dataList <- list()
        plotList <- list()
        for(i in seq_along(preds)){
                gene <- names(preds)[i]
                dataPlot <- cbind.data.frame(preds[[i]], data[, c("viable",
                                                                  "early_apoptotic",
                                                                  "late_apoptotic",
                                                                  "necrotic")])
                dataPlot$Host <- sapply(rownames(dataPlot),
                                        function(x) strsplit(x, "_")[[1]][2])
                dataPlot$Bacteria <- sapply(rownames(dataPlot),
                                        function(x) strsplit(x, "_")[[1]][1])
                bactFactLevels <- c("L5", "L6", "Bovis", "Chimp")[c("L5", "L6", "Bovis", "Chimp") %in% dataPlot$Bacteria]
                #print(bactFactLevels)
                dataPlot$Bacteria <- factor(dataPlot$Bacteria, 
                                          levels = unique(dataPlot$Bacteria))
                print(dataPlot$Bacteria)
                plotCols <- bactCols$color[match(levels(dataPlot$Bacteria),
                                                 bactCols$Bacteria)]
                genePlot <- ggplot(data = dataPlot, aes_string(x = gene, y = "pred") ) +
                        geom_line(size = 1) +
                        geom_point(aes_string(y = respVar, col = "Bacteria", shape = "Host")) + 
                        scale_discrete_manual("Bacteria",
                                              aesthetics = "colour",
                                              values = plotCols) +
                        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
                        geom_rug(sides = "b") +
                        theme_bw(base_size = 14) +
                        labs(x = gene,
                             y = yAxisTitle)
                plotList[[gene]] <- genePlot
                dataList[[gene]] <- dataPlot
        }
        globPlot <- ggarrange(plotlist = plotList, common.legend = T)
        return(list(data = dataList, 
                    plots = plotList, 
                    globPlot = globPlot))
}


doLogit <- function(variables, data){
        formula <- "host ~ "
        formula <- paste0(formula, 
                          paste0(variables, 
                                 collapse = " + "))
        formula <- as.formula(formula)
        logit <- glm(formula = formula,
                     data = data,
                     family = binomial)
        print(summary(logit))
        return(logit)
}

# Plot logistic regression results, only works when done in single explanatory variable
plotLogit <- function(model){
        bactCols <- data.frame(Bacteria = c("L5", "L6", "Bovis", "Chimp"),
                               color = c("#871414", "#24ad37", "#f279ce", "#c4bf62"),
                               stringsAsFactors = F)
        geneNames <- data.frame(locus = c("Rv2349c",
                                          "Rv2350c",
                                          "Rv2351c",
                                          "Mb1784c"),
                                gene = c("plcC",
                                         "plcB",
                                         "plcA",
                                         "plcD"))
        data <- model$model
        data$bacteria <- sapply(rownames(data), function(x) strsplit(x, "_")[[1]][1])
        variable = variable.names(model)[2]
        xAxisTitle = sprintf("*%s* (%s) expression relative to *sigA*",
                             geneNames$gene[match(variable, 
                                                  geneNames$locus)],
                             variable)
        yAxisTitle = "Probability of being a THP-1 infection"
        coef0 <- as.character(round(model$coefficients[2], digits = 3))
        coef1 <- as.character(round(model$coefficients[2], digits = 3))
        modSum <- summary(model)
        se <- as.character(round(modSum$coefficients[2, 2], digits = 3))
        zStat <- as.character(round(modSum$coefficients[2, 3], digits = 3))
        pVal <- as.character(round(modSum$coefficients[2, 4], digits = 3))
        statStr = sprintf(#"\u03B20 = %s    \u03B21 = %s    SE = %s    z-statistic = %s    p-value = %s",
                "\u03B2 = %s    SE = %s    z-statistic = %s    p-value = %s",
                          #coef0,
                          coef1,
                          se,
                          zStat,
                          pVal)
        plot <- ggplot(data, aes_string(x = variable,
                                 y = (as.numeric(data$host) - 1))) +
                geom_point(alpha=1,
                           aes(shape = host,
                               color = bacteria)) +
                scale_discrete_manual("Bacteria",
                                      aesthetics = "colour",
                                      values = bactCols$color[match(unique(data$bacteria),
                                                                    bactCols$Bacteria)]) +
                stat_smooth(method="glm", se=T, method.args = list(family=binomial)) +
                theme(axis.title.x = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      plot.caption=element_text(size = 8.5,
                                                hjust = 0,
                                                margin = margin(5,0,0,0))) +
                labs(x = xAxisTitle,
                     y = yAxisTitle,
                     caption = statStr)
        plot
}
