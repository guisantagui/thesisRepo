library(readxl)

pcrResDir <- "C:/Users/Guillem/Documents/PhD/comput/data/qPCR/results/"
resDir <- "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcPCRs/results/"
dataDir <- "C:/Users/Guillem/Documents/PhD/comput/data/"

oligoInfo <- as.data.frame(readxl::read_xlsx(paste0(dataDir, "RNA_info.xlsx")))
colnames(oligoInfo) <- make.names(colnames(oligoInfo))
oligoInfo <- oligoInfo[oligoInfo$RNA.FRACTION == "Prokaryote" | oligoInfo$RNA.FRACTION == "Total", ]
oligoInfo$INFECTION.DATE <- as.Date(oligoInfo$INFECTION.DATE, origin = "1899-12-30")
oligoInfo$SAMPLING.DATE <- as.Date(oligoInfo$SAMPLING.DATE, origin = "1899-12-30")

# Load functions
source("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcPCRs/scripts/pcrAnalFnctns.R")


pcrRes1 <- parsePcrRes(paste0(pcrResDir, "2022-02-21.xlsx"))
pcrRes2 <- parsePcrRes(paste0(pcrResDir, "2022-03-02.xlsx"))
pcrRes3 <- parsePcrRes(paste0(pcrResDir, "2022-03-11.xlsx"))





parsedRes1 <- parseSamps(pcrRes1, filtD1 = F, filtD3 = T)
parsedRes2 <- parseSamps(pcrRes2, filtD1 = F, filtD3 = T)
parsedRes3 <- parseSamps(pcrRes3, filtD1 = F, filtD3 = T)


View(parsedRes1)
View(parsedRes2)
View(parsedRes3)

parsedRes1$Ct[is.na(parsedRes1$Ct)] <- 0
parsedRes2$Ct[is.na(parsedRes2$Ct)] <- 0
parsedRes3$Ct[is.na(parsedRes3$Ct)] <- 0

CtDiff <- abs(parsedRes1$Ct - parsedRes2$Ct)

allRes <- parsedRes1
allRes$Ct_new <- parsedRes2$Ct
allRes$Ct_diff <- CtDiff
allRes$Ct_3 <- parsedRes3$Ct[match(parsedRes2$sampleName, parsedRes3$sampleName)]

        


colnames(allRes)[colnames(allRes) == "Ct"] <- "Ct_old"


allRes <- allRes[, c("sampleName", "oligo", "primerPair", "Ct_old", "Ct_new", "Ct_3", "Ct_diff", "n_num", "lin", "host", "daysPI")]


res_discrep <- allRes[allRes$Ct_diff >= 1, ]

write.csv(allRes, file = paste0(resDir, "plateDiffDayComp.csv"))
write.csv(allRes, file = paste0(resDir, "plateDiffDiscrep.csv"))
