####################################################################################################
####################################################################################################
#                                                                                                  #
# Do Mann Whitney U-test and, with the significantly differential metabolites, do metabolic        #
# pathway enrichment with FELLA.                                                                   #
#                                                                                                  #
####################################################################################################
####################################################################################################

if(!require(FELLA)) BiocManager::install("FELLA")
library(FELLA)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(igraph)) install.packages("igraph")
library(igraph)
if(!require(magrittr)) install.packages("magrittr")
library(magrittr)

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
fellaDir <- paste0(plotDir, "FELLA/")

# Create directories if necessary
if(!dir.exists(plotDir)){
        dir.create(plotDir)
}
if(!dir.exists(outDir)){
        dir.create(outDir)
}
if(!dir.exists(fellaDir)){
        dir.create(fellaDir)
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


####################################################################################################
#                                                                                                  #
# Run univariate analysis (Mann Whitney U test).                                                   #
#                                                                                                  #
####################################################################################################

# Add a column to the metabolomics data indicating the rhamnolipid production phenotype 
# (2 categories, 0 = non-producer, 1 = producer).
mets_ccmn$rhamn <- rhamnMat$rhamn2cats[match(gsub("\\_.*", 
                                                  "", 
                                                  rownames(mets_ccmn)), 
                                             rhamnMat$strains)]

# Do Mann Whitney U test
mwRhamn_pVal <- apply(mets_ccmn[, colnames(mets_ccmn) != "rhamn"],
                      2,
                      function(x) wilcox.test(x=x[mets_ccmn$rhamn == 0], 
                                              y=x[mets_ccmn$rhamn == 1])$p.value)

# Adjust p-values with Benjamini-Hochberg method
mwRhamn_pAdj <- p.adjust(mwRhamn_pVal, method = "BH")

# Filter significative metabolites, with a significance value of 0.05 and write results
mwRhamn_pAdj_sign <- mwRhamn_pAdj[mwRhamn_pAdj < 0.05]

mwRhamn_result <- data.frame(metabolite = names(mwRhamn_pAdj_sign),
                             pVal_adj =mwRhamn_pAdj_sign)

write.csv(mwRhamn_result, 
          file = paste0(outDir, "mwRhamn_result.csv"),
          row.names = F)

####################################################################################################
#                                                                                                  #
# Run FELLA. Note: results might change slightly, as KEGG gets updated frequently.                 #
#                                                                                                  #
####################################################################################################

# Get KEGG IDs of the differential metabolites between rhamnolipid producers and non-producers
signMetKEGGIDs <- dictionary$KEGG.IDs[match(names(mwRhamn_pAdj_sign), 
                                            dictionary$KyuRheeNames)]

# Build graph of KEGG entries for Pseudomonas aeruginosa strain PA14 ("pau" code in KEGG). Exclude 
# From graph generic pathway entries ("Metabolic pathways", "Carbon metabolism", 
# "2-Oxocarboxylic acid metabolism", "Biosynthesis of amino-acids").
graph <- buildGraphFromKEGGREST(
        organism = "pau",
        filter.path = c("01100", 
                        "01200", 
                        "01210", 
                        #"01212", 
                        "01230",
                        "02010")
)

# Write KEGG knowledge model based on the graph
buildDataFromGraph(
        keggdata.graph = graph,
        databaseDir = NULL,
        internalDir = TRUE,
        matrices = "none",
        normality = "diffusion",
        niter = 100)


# Loads the last KEGG model stored in the computer into fella.data object.
fella.data <- loadKEGGdata(
        internalDir = T,
        loadMatrix = "none")

fella.data

# Get ids of compounds, reactions and enzymes in the fella.data
id.cpd <- getCom(fella.data, level = 5, format = "id") %>% names
id.rx <- getCom(fella.data, level = 4, format = "id") %>% names
id.ec <- getCom(fella.data, level = 3, format = "id") %>% names

# Run FELLA
fellaRes<- enrich(
        compounds = signMetKEGGIDs,
        data = fella.data,
        method = "diffusion",
        approx = "normality")

# Test if, among the compounds in FELLA input some were not mapped to graph.
# All the compounds are mapped.
excludedComps <- getExcluded(fellaRes)
if(length(excludedComps) == 0){
        print("Any compunds in FELLA input were excpluded from the graph.")
}else{
        print(sprintf("Compounds excluded from FELLA graph: %s.",
                      paste(excludedComps, collapse = ", ")))
}

# Plot FELLA results graph
graph_res <- generateResultsGraph(
        object = fellaRes,
        method = "diffusion",
        nlimit = 350,
        data = fella.data)

pdf(paste0(fellaDir, "FELLA_graph_res.pdf"), width = 15, height = 15)
plotGraph(
        graph_res
)
dev.off()

tab_res <- generateResultsTable(
        object = fellaRes,
        data = fella.data,
        method = "diffusion",
        nlimit = 500)

write.csv(tab_res, file = paste0(outDir, "FELLA_tab_res.csv"))