# This scripts builds dataframes with the info of the genes removed from L4 for generating each 
# lineage model

dataDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/data/"
delDir <- paste0(dataDir, "deletions/")
snpDir <- paste0(dataDir, "snps/")
whatMod <- "delsAllSNPs"
delPropDir <- paste0(delDir, "delProp/")
provDir <- paste0(dataDir, "proveanRes/")
# Deletions 

outDir <- sprintf("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbMods_paperGSMMs/results/built_models/%s_mods/remModDFs/", whatMod)

if(!dir.exists(outDir)){
        dir.create(outDir)
}

delFile <- paste0(delDir, "dels_toRemDF.csv")
stopGainSNPFile <- paste0(snpDir, "stopGain_toRemDF.csv")
deletSNPsFile <- paste0(snpDir, "deletSNPs.csv")
missProvSNPFile <- paste0(snpDir, "missProv_toRemDF.csv")
provSignFile <- paste0(provDir, "proveanSign.csv")
delMatFile <- paste0(gsub("data/", 
                          "results/", 
                          dataDir), 
                     sprintf("%s_mat_cur.csv", 
                             gsub("ls", 
                                  "l_", 
                                  whatMod)))


lins <- c("A1", "A2", "A3", "A4", "L1", "L2", "L3", "L5", "L6", "L7", "L8", "L9")


dels <- read.csv(delFile, row.names = 1, stringsAsFactors = F)
stopGainSNPs <- read.csv(stopGainSNPFile, row.names = 1, stringsAsFactors = F)
missProvSNPs <- read.csv(missProvSNPFile, row.names = 1, stringsAsFactors = F)
deletSNPs <- read.csv(deletSNPsFile, row.names = 1, stringsAsFactors = F)
provSign <- read.csv(provSignFile, row.names = 1, stringsAsFactors = F)
delMat <- read.csv(delMatFile, row.names = 1, stringsAsFactors = F)
for(i in seq_along(lins)){
        lin <- lins[i]
        remFromMod <- delMat[, lin]
        remFromMod <- rownames(delMat)[remFromMod == 0]
        linDelProp <- read.csv(sprintf("%s%sCSSignDelProp.csv", delPropDir, lin), 
                               row.names = 1, 
                               stringsAsFactors = F)
        delsLin <- dels[, grep(lin, colnames(dels))]
        colnames(delsLin) <- gsub(".*\\_", 
                                  "", 
                                  colnames(delsLin))
        delsLin <- delsLin[!is.na(delsLin$locus), ]
        
        stopGaLin <- stopGainSNPs[, grep(lin, colnames(stopGainSNPs))]
        colnames(stopGaLin) <- gsub(".*\\_", 
                                    "", 
                                    colnames(stopGaLin))
        stopGaLin <- stopGaLin[!is.na(stopGaLin$locus), ]
        
        provMiLin <- missProvSNPs[, grep(lin, colnames(missProvSNPs))]
        colnames(provMiLin) <- gsub(".*\\_", 
                                    "", 
                                    colnames(provMiLin))
        provMiLin <- provMiLin[!is.na(provMiLin$locus), ]

        linProps <- c()
        for(j in seq_along(delsLin$locus)){
                locus <- delsLin$locus[j]
                locLinProp <- linDelProp[lin, sprintf("%s_15", locus)]
                linProps <- c(linProps, locLinProp)
        }
        names(linProps) <- as.character(delsLin$locus)
        delsLin$delProp <- linProps
        delsNotRem <- delsLin[!delsLin$locus %in% remFromMod, ]
        delsLin <- delsLin[delsLin$locus %in% remFromMod, ]
        
        if(nrow(delsLin) > 0){
                rownames(delsLin) <- 1:nrow(delsLin)
                write.csv(delsLin, file = sprintf("%s%s_modDels.csv", outDir, lin))
                print(sprintf("%s_modDels.csv saved at %s", lin, outDir))
        }
        if(nrow(delsNotRem) > 0){
                rownames(delsNotRem) <- 1:nrow(delsNotRem)
                write.csv(delsNotRem, file = sprintf("%s%s_modDels_resc.csv", outDir, lin))
                print(sprintf("%s_modDels_resc.csv saved at %s", lin, outDir))
        }
        if(nrow(stopGaLin) > 0){
                linSGSNPs <- data.frame(matrix(nrow = 0, 
                                               ncol = 5, 
                                               dimnames = list(NULL, 
                                                               c("locus", 
                                                                 "gene", 
                                                                 "product", 
                                                                 "variant", 
                                                                 "propSNP"))))
                for(j in seq_along(stopGaLin$locus)){
                        locus <- stopGaLin$locus[j]
                        delSGSnps <- deletSNPs[deletSNPs$LOCUS == locus & deletSNPs$LIN == lin, ]
                        delSGSnps <- delSGSnps[delSGSnps$propSNP >= 0.85, ]
                        delSGSnps <- delSGSnps[grep("*", delSGSnps$AA_one, fixed = T), ]
                        delSGSnps$gene <- rep(stopGaLin$gene[match(locus, 
                                                                   stopGaLin$locus)], 
                                              nrow(delSGSnps))
                        delSGSnps$product <- rep(stopGaLin$product[match(locus, 
                                                                         stopGaLin$locus)], 
                                                 nrow(delSGSnps))
                        delSGSnps <- delSGSnps[, c("LOCUS", "gene", "product", "AA_one", "propSNP")]
                        colnames(delSGSnps) <- c("locus", "gene", "product", "variant", "propSNP")
                        linSGSNPs <- rbind.data.frame(linSGSNPs, delSGSnps)
                }
                linSGSNPsNotRem <- linSGSNPs[!linSGSNPs$locus %in% remFromMod, ]
                linSGSNPs <- linSGSNPs[linSGSNPs$locus %in% remFromMod, ]
                if(nrow(linSGSNPs) > 0){
                        linSGSNPs <- linSGSNPs[!duplicated(paste(linSGSNPs$locus, 
                                                                 linSGSNPs$variant, 
                                                                 sep = ".")), ]
                        rownames(linSGSNPs) <- 1:nrow(linSGSNPs)
                        write.csv(linSGSNPs, file = sprintf("%s%s_modSGSNPs.csv", outDir, lin))
                        print(sprintf("%s_modSGSNPs.csv saved at %s", lin, outDir))
                }
                if(nrow(linSGSNPsNotRem) > 0){
                        linSGSNPsNotRem <- linPMSNPsNotRem[!duplicated(paste(linSGSNPsNotRem$locus, 
                                                                             linSGSNPsNotRem$variant, 
                                                                             sep = ".")), ]
                        rownames(linSGSNPsNotRem) <- 1:nrow(linSGSNPsNotRem)
                        write.csv(linSGSNPsNotRem, file = sprintf("%s%s_modSGSNPs_resc.csv", outDir, lin))
                        print(sprintf("%s_modSGSNPs_resc.csv saved at %s", lin, outDir))
                }
        }else{
                print(sprintf("%s doesn't have any stop gained SNP at a proportion equal or higher than 85%%.", lin))
        }
        if(nrow(provMiLin) > 0){
                linPMSNPs <- data.frame(matrix(nrow = 0, 
                                               ncol = 6, 
                                               dimnames = list(NULL, 
                                                               c("locus", 
                                                                 "gene", 
                                                                 "product", 
                                                                 "variant", 
                                                                 "propSNP",
                                                                 "provean"))))
                for(j in seq_along(provMiLin$locus)){
                        locus <- provMiLin$locus[j]
                        delPMSnps <- deletSNPs[deletSNPs$LOCUS == locus & deletSNPs$LIN == lin, ]
                        delPMSnps <- delPMSnps[delPMSnps$propSNP >= 0.85, ]
                        remSG <- grep("*", delPMSnps$AA_one, fixed = T)
                        if(length(remSG) > 0){
                                delPMSnps <- delPMSnps[-remSG, ]
                        }
                        delPMSnps$gene <- rep(provMiLin$gene[match(locus, 
                                                                   provMiLin$locus)], 
                                              nrow(delPMSnps))
                        delPMSnps$product <- rep(provMiLin$product[match(locus, 
                                                                         provMiLin$locus)], 
                                                 nrow(delPMSnps))
                        delPMSnps$provean <- provSign$PROVEAN[match(delPMSnps$snpCode, paste(provSign$LOCUS,
                                                                                             provSign$AA_one,
                                                                                             sep = "."))]
                        delPMSnps <- delPMSnps[, c("LOCUS", "gene", "product", "AA_one", "propSNP", "provean")]
                        colnames(delPMSnps) <- c("locus", "gene", "product", "variant", "propSNP", "provean")
                        linPMSNPs <- rbind.data.frame(linPMSNPs, delPMSnps)
                }
                linPMSNPsNotRem <- linPMSNPs[!linPMSNPs$locus %in% remFromMod, ]
                linPMSNPs <- linPMSNPs[linPMSNPs$locus %in% remFromMod, ]
                
                if(nrow(linPMSNPs) > 0){
                        linPMSNPs <- linPMSNPs[!duplicated(paste(linPMSNPs$locus, 
                                                                 linPMSNPs$variant, 
                                                                 sep = ".")), ]
                        rownames(linPMSNPs) <- 1:nrow(linPMSNPs)
                        write.csv(linPMSNPs, file = sprintf("%s%s_modPMSNPs.csv", outDir, lin))
                        print(sprintf("%s_modPMSNPs.csv saved at %s", lin, outDir))
                }
                if(nrow(linPMSNPsNotRem) > 0){
                        linPMSNPsNotRem <- linPMSNPsNotRem[!duplicated(paste(linPMSNPsNotRem$locus, 
                                                                             linPMSNPsNotRem$variant, 
                                                                             sep = ".")), ]
                        rownames(linPMSNPsNotRem) <- 1:nrow(linPMSNPsNotRem)
                        write.csv(linPMSNPsNotRem, file = sprintf("%s%s_modPMSNPs_resc.csv", outDir, lin))
                        print(sprintf("%s_modPMSNPs_resc.csv saved at %s", lin, outDir))
                }
        }else{
                print(sprintf("%s doesn't have any PROVEAN significative SNP at a proportion equal or higher than 85%%.", lin))
        }
}
