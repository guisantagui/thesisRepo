# This script compares the four annotation files we have, builds a version with 
# all the genes that appear in the four and outputs a tab separated txt file to 
# be used as input for the deletion detection script. 
if(!require(rtracklayer)) install.packages("rtracklayer")
library(rtracklayer)

setwd("/storage/PGO/scripts/geneTxt4Piplne_globGeneDict_create/")

annotDir <- "/storage/PGO/data/mtb/annotations/"
txtOutDir <- "/storage/PGO/data/mtb/annotations/genesTXT4Pipeline/"
annDictDir <- "/storage/PGO/data/mtb/annotations/annotDict/"

grntAnnot <- as.data.frame(readGFF(paste(annotDir, "genes.gff", sep = "")))
ncbiAnnot <- as.data.frame(readGFF(paste(annotDir, "h37RvNCBIAnnot.gff3", sep = "")))
myc1Annot <- as.data.frame(readGFF(paste(annotDir, "Mycobacterium_tuberculosis_H37Rv_gff_v1.gff", sep = "")))
myc2Annot <- as.data.frame(readGFF(paste(annotDir, "Mycobacterium_tuberculosis_H37Rv_gff_v2.gff", sep = "")))


# Create file txt gene position file in the format of the pipeline 
genes <- grntAnnot

genes4txt <- genes[!is.na(genes$locus_tag), c("start", "end", "strand", "locus_tag")]

ncbiAnnot <- ncbiAnnot[1:nrow(ncbiAnnot) %% 2 == 1 & !is.na(ncbiAnnot$locus_tag), ]

ncbiAnnot2Add <- ncbiAnnot[!ncbiAnnot$locus_tag %in% genes4txt$locus_tag, c("start", "end", "strand", "locus_tag")]

genes4txt <- rbind.data.frame(genes4txt, ncbiAnnot2Add)

genes4txt <- genes4txt[order(genes4txt$start), ]

genes4txt$strand <- as.factor(sapply(genes4txt$strand, function(x) if(x == "+") x <- "F" else x <- "R"))

write.table(genes4txt, 
            file = paste(txtOutDir, "allGenesMTBC.txt", sep = ""), 
            sep = "\t", 
            quote = F, 
            row.names = F, 
            col.names = F)
            
print(paste("Txt file of all the h37Rv genes and coordinates has been generated at", 
            paste(txtOutDir, "allGenesMTBC.txt", sep = ""), sep = " "))

# Create a global dictionary to know coordinates, product and strand of each gene

genes4dict <- genes[!is.na(genes$locus_tag), c("start", "end", "strand", "locus_tag", "Name")]

colnames(genes4dict)[ncol(genes4dict)] <- "gene"

genes4dict$gene <- sapply(genes4dict$gene, function(x) strsplit(x, "_")[[1]][1])

genes4dict$product <- genes$product[match(genes4dict$locus_tag, genes$locus_tag) + 1]

ncbiAnnot2Add2Dict <- ncbiAnnot[!ncbiAnnot$locus_tag %in% genes4dict$locus_tag, c("start", 
                                                                                  "end", 
                                                                                  "strand", 
                                                                                  "locus_tag", 
                                                                                  "gene", 
                                                                                  "product")]

genes4dict <- rbind.data.frame(genes4dict, ncbiAnnot2Add2Dict)

genes4dict <- genes4dict[order(genes4dict$start), ]

write.csv(genes4dict, file = paste(annDictDir, "allGenesMTBC_dict.csv", sep = ""), row.names = F)

print(paste("Dictionary csv file of all the h37Rv genes and information has been generated at", 
            paste(annDictDir, "allGenesMTBC_dict.csv", sep = ""), sep = " "))