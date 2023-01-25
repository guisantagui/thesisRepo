####################################################################################################
####################################################################################################
#                                                                                                  #
# Run supervised analysis (OPLS-DA), training for rhamnolipid production phenotype.                #                                                      #
#                                                                                                  #
####################################################################################################
####################################################################################################

if(!require(ropls)) BiocManager::install("ropls")
library(ropls)
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(ggpubr)) install.packages("ggpubr")
library(ggpubr)
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(ggrepel)) install.packages("ggrepel")
library(ggrepel)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(randomcoloR)) install.packages("randomcoloR")
library(randomcoloR)

# Directory paths
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
fig5Dir <- paste0(plotDir, "fig5/")
figS7Dir <- paste0(plotDir, "figS7/")
oplsdaDir <- paste0(plotDir, "OPLS-DA/")
# Create directories if necessary
if(!dir.exists(plotDir)){
        dir.create(plotDir)
}
if(!dir.exists(outDir)){
        dir.create(outDir)
}
if(!dir.exists(fig5Dir)){
        dir.create(fig5Dir)
}
if(!dir.exists(figS7Dir)){
        dir.create(figS7Dir)
}
if(!dir.exists(oplsdaDir)){
        dir.create(oplsdaDir)
}

# Functions
####################################################################################################

# KEGGREST is giving problems, so here's a set of functions that make the same
keggLink_cust <- function(target, source){
        if(!require(RCurl)) install.packages("RCurl")
        library(RCurl)
        URL <- sprintf("https://rest.kegg.jp/link/%s/%s", target, source)
        tab <- getURL(URL)
        if(nchar(tab) > 1){
                outDF <- read.table(URL, quote="", sep="\t")
                out <- as.character(outDF[, 2])
                names(out) <- outDF[, 1]
        }else{
                out <- NULL
        }
        
        return(out)
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

# CCMN-normalized metabolomic data (.tsv file), without putative metabolites.
mets_ccmn <- read.table(file = paste0(outDir, "mets_ccmn_noPut.tsv"), sep = "\t")
# When tsv is imported R changes the column names. Correct it with dictionary
colnames(mets_ccmn) <- dictionary$KyuRheeNames[match(colnames(mets_ccmn),
                                                     make.names(dictionary$KyuRheeNames))]

# Rhamnolipid production data
rhamnMat <- read.csv(file = paste0(dataDir, "rhamnMat.csv"),
                     row.names = 1, 
                     stringsAsFactors = F)

# FELLA results table
fella_res <- read.csv(file = paste0(outDir, "fella_tab_res.csv"), 
                      row.names = 1,
                      stringsAsFactors = F)

# Differential metabolites between rhamnolipid producers and non-producers, determined
# by univariate test (Mann Whitney U test)
mwRhamn_sign <- read.csv(file = paste0(outDir, "mwRhamn_result.csv"),
                         stringsAsFactors = F)

# Pathway colors for plotting
pathCols <- read.csv(paste0(dataDir, "pathway_colors.csv"),
                     row.names = 1,
                     stringsAsFactors = F)

####################################################################################################
#                                                                                                  #
# Run supervised analysis (OPLS-DA).                                                               #
#                                                                                                  #
####################################################################################################

# Add a column to the metabolomics data indicating the rhamnolipid production phenotype 
# (2 categories, 0 = non-producer, 1 = producer).
mets_ccmn$rhamn <- as.factor(rhamnMat$rhamn2cats[match(gsub("\\_.*", 
                                                            "", 
                                                            rownames(mets_ccmn)), 
                                                 rhamnMat$strains)])

# Run OPLS-DA with 2000 label permutations for obtaining the p-value (significance of model).
# Save the output image that ropls produces
pdf(file = paste0(oplsdaDir, "OPLSDA_rhamn.pdf"), width = 10, height = 10)
OPLSDA_rhamn <- opls(mets_ccmn[, colnames(mets_ccmn) != "rhamn"],
                     mets_ccmn$rhamn, 
                     predI = 1, 
                     orthoI = 3,
                     permI = 2000)
dev.off()

# Plot OPLS-DA scores for each strain (Figure 5A)
####################################################################################################

# Create vector of strains
strains <- unique(gsub("\\_.*|(PA14).*", 
                       "\\1", 
                       oplsdaScores$sample))

# Obtain scores from OPLSDA object
oplsdaScores <- data.frame(sample = rownames(mets_ccmn),
                           t1 = OPLSDA_rhamn@scoreMN[, 1],
                           o1 = OPLSDA_rhamn@orthoScoreMN[, 1],
                           stringsAsFactors = F)

# obtain average scores dataframe for each strain
oplsdaScores_mean <- data.frame(matrix(ncol = 3, 
                                       nrow = 0, 
                                       dimnames = list(NULL, 
                                                       c("strain", 
                                                         "t1", 
                                                         "o1"))))

for(i in seq_along(strains)){
        strain <- strains[i]
        subScores <- oplsdaScores[grep(strain, oplsdaScores$sample), ]
        toBind <- data.frame(strain = strain, 
                             t1 = mean(subScores$t1),
                             o1 = mean(subScores$o1),
                             stringsAsFactors = F)
        oplsdaScores_mean <- rbind.data.frame(oplsdaScores_mean, toBind)
}

oplsdaScores_mean$rhamn <- as.factor(rhamnMat$rhamn2cats[match(oplsdaScores_mean$strain, 
                                                               gsub("(PA14).*", 
                                                                    "\\1", 
                                                                    rhamnMat$strains))])

# Do score plot of OPLS-DA (Figure 5A)
fig5A_plot <- ggplot(oplsdaScores_mean, aes(x=t1, y=o1, color=rhamn, label = strain)) + 
        geom_point(size = 5) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) + 
        scale_discrete_manual("Rhamnolipid production",
                              aesthetics = "colour",
                              values = unique(rhamnMat$rhamn2catsCols[order(rhamnMat$rhamn2cats)]),
                              labels = c("Rhamnolipid -", "Rhamnolipid +/- and +")) +
        geom_text_repel(size = 5,
                        point.padding = 0.5) + 
        #labs(col="Rhamnolipid Production") + 
        theme(axis.text.y   = element_text(size=14),
              axis.text.x   = element_text(size=14),
              axis.title.y  = element_text(size=14),
              axis.title.x  = element_text(size=14),
              panel.background = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA, size=1)
        ) +
        xlab(sprintf("Predictive component (%s%%)", round(OPLSDA_rhamn@modelDF[1, 1]*100))) +
        ylab("First orthogonal component")

fig5A_plot

ggsave(filename = paste0(fig5Dir, "fig5_A.pdf"),
       height = 5, width = 7)

# Do barplots of metabolite loadings for predictive component, with the pathways where they
# are involved indicated (Figures 5B and S7).
####################################################################################################

# Create dataframe with OPLS-DA metabolite loadings of predictor component and obtain metabolite 
# KEGG IDs
oplsdaLoads <- data.frame(metabolite = rownames(OPLSDA_rhamn@loadingMN),
                          keggID = dictionary$KEGG.IDs[match(rownames(OPLSDA_rhamn@loadingMN), 
                                                                      dictionary$KyuRheeNames)],
                          loading = OPLSDA_rhamn@loadingMN[, 1],
                          stringsAsFactors = F)

# Obtain the pathway each metabolite is involved in, among the ones enriched in FELLA
pathsPerMet <- list()
for(i in seq_along(oplsdaLoads$keggID)){
        keggID <- oplsdaLoads$keggID[i]
        if(!is.na(keggID)){
                # Obtain all the pathways each metabolite is involved in
                #paths <- keggLink("pathway", 
                #                  keggID)
                paths <- keggLink_cust("pathway", 
                                       keggID)
                # Filter to keep just the pathways we found enriched with FELLA.
                paths <- paths[gsub("path:map", "pau", paths) %in% fella_res$KEGG.id]
        }else{
                paths <- NA
        }
        pathsPerMet[[i]] <- paths
}
names(pathsPerMet) <- oplsdaLoads$metabolite
pathsPerMet

allPaths <- unique(unlist(pathsPerMet))
allPaths <- allPaths[!is.na(allPaths)]


# Get a list of which of our metabolites is involved in each one of the FELLA pathways
metsInPaths <- list()
for(i in seq_along(allPaths)){
        pathFella <- allPaths[i]
        #metsInPath <- keggLink("compound", pathFella)
        metsInPath <- keggLink_cust("compound", pathFella)
        metsInPath <- gsub("cpd:", "", metsInPath)
        mets <- metsInPath[metsInPath %in% oplsdaLoads$keggID]
        metsInPaths[[i]] <- mets
}
names(metsInPaths) <- allPaths

# Add the metabolites that don't have a role in any of the FELLA enriched pathways
metsInPaths$anyRelPath <- oplsdaLoads$keggID[!oplsdaLoads$keggID %in% unique(unlist(metsInPaths))]

# Create a loading dataframe where metabolites are stacked according to what pathway they are 
# involved in 
loadsPaths <- data.frame()
for(i in seq_along(names(metsInPaths))){
        pathID <- gsub("path:map", "pau", names(metsInPaths)[i])
        if(pathID == "anyRelPath"){
                pathName <- "Not Mapped"
        }else{
                pathName <- fella_res$KEGG.name[fella_res$KEGG.id == pathID]
                pathName <- gsub(" -....*| - Pse.*", "", pathName)
        }
        mets <- metsInPaths[[i]]
        toBind <- oplsdaLoads[match(mets, oplsdaLoads$keggID), ]
        toBind$pathID <- rep(pathID, nrow(toBind))
        toBind$pathName <- rep(pathName, nrow(toBind))
        toBind <- toBind[order(toBind$loading), ]
        rownames(toBind) <- 1:nrow(toBind)
        loadsPaths <- rbind.data.frame(loadsPaths, toBind)
}

# Make metabolites unique, as many of them are multiple times because they participate in
# different metabolic pathways. Factor too, to keep order in barplot.
loadsPaths$metabolite <- factor(make.unique(loadsPaths$metabolite),
                                levels = make.unique(loadsPaths$metabolite))

loadsPaths$pathName <- factor(loadsPaths$pathName,
                              levels = unique(loadsPaths$pathName))

# Plot barplot of all loadings stacked according to the pathway they are involved in
# (Figure S7).
####################################################################################################

figS7_plot <- ggplot(data = loadsPaths,
                     aes(x = metabolite, y = loading, fill = factor(pathName, 
                                                                    levels = unique(pathName)))) +
        labs(x = "Metabolite", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() + 
        scale_fill_manual(values = pathCols$color[match(levels(loadsPaths$pathName), 
                                                        pathCols$pathway)],
                          name = "Pathway") +
        scale_x_discrete(position = "top") +
        theme(axis.text = element_text(size = 14))

figS7_plot

ggsave(filename = paste0(figS7Dir, "figS7.pdf"), height = 20, width = 12)


# Create filtered dataset to focus just on carbon (TCA cycle and Pentose Phosphate 
# Pathway) and amino acid metabolism, for avoiding having repeated compounds and 
# for the sake of simplicity in figure 5 
paths2Keep <- c("Valine, leucine and isoleucine biosynthesis", 
                "Alanine, aspartate and glutamate metabolism", 
                "Cysteine and methionine metabolism",
                "Pentose phosphate pathway",
                "Citrate cycle (TCA cycle)")

AA_paths <- c("Valine, leucine and isoleucine biosynthesis",
              "Alanine, aspartate and glutamate metabolism", 
              "Cysteine and methionine metabolism")

loadsPathsFilt <- loadsPaths[loadsPaths$pathName %in% paths2Keep, ]

# Order for obtaining a dataframe where all the pathways related either with 
# carbon or amino-acid metabolism are contiguous. 
ordVec <- unlist(sapply(paths2Keep, function(x) which(loadsPathsFilt$pathName %in% x)))
loadsPathsFilt <- loadsPathsFilt[ordVec, ]

# Aggrupate all the pathways related to amino acids into "Amino acid metabolism"
# category. Drop pathway KEGG ID, as it no more is necessary. 
loadsPathsFilt$pathName <- gsub(paste(AA_paths, collapse = "|"), 
                                "Amino acid metabolism", 
                                loadsPathsFilt$pathName)

# Change the name of the TCA cycle 
loadsPathsFilt$pathName <- gsub("Citrate cycle (TCA cycle)", 
                                "Tricarboxylic acid (TCA) cycle", 
                                loadsPathsFilt$pathName, 
                                fixed = T)

loadsPathsFilt <- loadsPathsFilt[, c("metabolite",
                                     "keggID",
                                     "loading",
                                     "pathName")]

# Remove .N from metabolite names to see what are the duplicates.
loadsPathsFilt$metabolite <- gsub("\\..*", "", loadsPathsFilt$metabolite)

dupMets <- loadsPathsFilt$metabolite[duplicated(loadsPathsFilt$metabolite)]

dupMetsAAMet <- dupMets[dupMets %in% loadsPathsFilt$metabolite[loadsPathsFilt$pathName == "Amino acid metabolism"]]
dupMetsTCACyc <- dupMets[dupMets %in% loadsPathsFilt$metabolite[loadsPathsFilt$pathName == "Tricarboxylic acid (TCA) cycle"]]
dupMetsPPP <- dupMets[dupMets %in% loadsPathsFilt$metabolite[loadsPathsFilt$pathName == "Pentose phosphate pathway"]]




# Aspartate and Alanine are duplicated because they appear in multiple amino acid 
# related pathways. Now that we merged them in a single category we can remove 
# them. 

loadsPathsFilt <- loadsPathsFilt[(!duplicated(loadsPathsFilt$metabolite) & 
                                          loadsPathsFilt$pathName == "Amino acid metabolism") | 
                                         (loadsPathsFilt$pathName == "Tricarboxylic acid (TCA) cycle" | 
                                                  loadsPathsFilt$pathName == "Pentose phosphate pathway"), ]


loadsPaths$pathName[grep("Citrate", loadsPaths$metabolite)]
loadsPaths$pathName[grep("Succinate", loadsPaths$metabolite)]
loadsPaths$pathName[grep("Oxoglutarate", loadsPaths$metabolite)]
loadsPaths$pathName[grep("Fumarate", loadsPaths$metabolite)]

# Citrate, Succinate, Oxoglutarate and Fumarate map also to Amino acid metabolism, in 
# particular to "Alanine, aspartate and glutamate metabolism" in KEGG. But are 
# peripherial in this pathway, while in TCA-cycle are central. So for avoiding 
# duplicate metabolites in Figure 5B we are removing them. 

loadsPathsFilt <- loadsPathsFilt[(!loadsPathsFilt$metabolite %in% dupMetsTCACyc & 
                                          loadsPathsFilt$pathName == "Amino acid metabolism") | 
                                         (loadsPathsFilt$pathName == "Tricarboxylic acid (TCA) cycle" | 
                                                  loadsPathsFilt$pathName == "Pentose phosphate pathway"), ]

# Order the rows of Amino acid metabolism according to loading value
loadPathsFilt_ordVec <- c(order(loadsPathsFilt$loading[loadsPathsFilt$pathName == "Amino acid metabolism"]),
                          which(loadsPathsFilt$pathName %in% c("Tricarboxylic acid (TCA) cycle",
                                                               "Pentose phosphate pathway")))

loadsPathsFilt <- loadsPathsFilt[loadPathsFilt_ordVec, ]

# Initialize row names, add asterisks to significant metabolites in the Mann Whitney U test and 
# assign new factors to metabolite names to keep order for the barplot. 
rownames(loadsPathsFilt) <- 1:nrow(loadsPathsFilt)

signMets <- paste0(loadsPathsFilt$metabolite[loadsPathsFilt$metabolite %in% mwRhamn_sign$metabolite], 
                   " (*)")

loadsPathsFilt$metabolite[loadsPathsFilt$metabolite %in% mwRhamn_sign$metabolite] <- signMets

loadsPathsFilt$metabolite <- factor(loadsPathsFilt$metabolite, 
                                    levels = loadsPathsFilt$metabolite)

# Plot Figure 5B
####################################################################################################

# Create color vector
fig5B_cols <- c("#ff732b", "#2b7fff", "#00b819")

fig5B_plot <- ggplot(data = loadsPathsFilt,
                     aes(x = metabolite, y = loading, fill = pathName)) +
        labs(x = "Metabolite", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity", position = "dodge") +
        geom_vline(xintercept = 0,
                   #linetype = "dashed", 
                   alpha = 0.8) +
        coord_flip() + 
        scale_fill_manual(values = fig5B_cols,
                          name = "Pathway") +
        scale_x_discrete(position = "top") +
        theme(axis.text.y = element_text(size=20),
              axis.text.x = element_text(size=15),
              panel.background = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.line.y = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA, size=1))

fig5B_plot

ggsave(filename = paste0(fig5Dir, "fig5_B.pdf"), height = 10, width = 9)
