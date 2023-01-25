print("Loading SNPMat")
snpMat <- read.csv("/storage/PGO/results/mtb/snpsPipe/deletSNPsMat.csv", row.names = 1)

print(length(colnames(snpMat)))

print(length(unique(gsub("*..", "", colnames(snpMat)))))
