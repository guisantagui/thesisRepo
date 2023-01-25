#########################################################################################################
#                                                                                                       #
# Plots a dendrogram that represents simplified phylogenic relationships between Mtb lineages. (The     #
# tree for constructing the dendrogram was built by hand).                                              #
#                                                                                                       #
#########################################################################################################
if(!require(argparser)) install.packages("argparser", repos='http://cran.us.r-project.org')
library(argparser)
if(!require(ggplot2)) install.packages("ggplot2", repos='http://cran.us.r-project.org')
library(ggplot2)
if(!require(ggtree)) BiocManager::install("ggtree", update = F)
library(ggtree)
if(!require(phylobase)) install.packages("phylobase")
library(phylobase)
if(!require(phytools)) install.packages("phytools")
library(phytools)

# Create parser, which accepts arguments from terminal and pass to the script
parser <- arg_parser("This script plots a simplified phylogeny of Mtb.")

parser <- add_argument(parser = parser,
                       arg = c("--treeFile", 
                               "--output"),
                       help = c("Path to the tree object for building the phylogeny.",
                                "Output directory where the plot will be stored"),
                       flag = c(F, F),
                       default = list("--treeFile" = "/home/guisana/scripts/delsPipe/",
                                      "--output" = "/storage/PGO/results/mtb/delsPipe/"))
                                      
parsed <- parse_args(parser)

# Directory stuff
#########################################################################################################
treeDir <- parsed$treeFile
outDir <- parsed$output

# create color dataframe and load the tree 
#########################################################################################################

colors <- data.frame(name = c("rosa",
                              "azul",
                              "morado",
                              "rojo",
                              "rojo_oscuro",
                              "verde_normal",
                              "amarillo",
                              "naranja",
                              "verde_radioactivo",
                              "amarillo_raro",
                              "verde_pastel",
                              "azul_pastel",
                              "rosa_pastel"),
                     hex = c("#ff30c1",
                             "#001aff",
                             "#8826b5",
                             "#ff0000",
                             "#871414",
                             "#24ad37",
                             "#fbff00",
                             "#ff9d00",
                             "#37ff30",
                             "#c4bf62",
                             "#87e0a0",
                             "#6db3d6",
                             "#f279ce"), 
                     lin = c("L1",
                             "L2",
                             "L3",
                             "L4",
                             "L5",
                             "L6",
                             "L7",
                             "L8",
                             "L9",
                             "A1",
                             "A2",
                             "A3",
                             "A4"),
                     RGB = c("255,48,193",
                             "0,26,255",
                             "136,38,181",
                             "255,0,0",
                             "135,20,20",
                             "36,173,55",
                             "251,255,0",
                             "255,157,0",
                             "55,255,48",
                             "196,191,98",
                             "135,224,160",
                             "109,179,214",
                             "242,121,206"),
                     stringsAsFactors = F)

linDend <- ape::read.tree(paste(treeDir, "linTree", sep = ""))

# Plot the tree
#########################################################################################################

linPhylo <- as(linDend, "phylo4")
linPhyloDF <- data.frame(lin = linDend$tip.label)

linTreeCols <- data.frame(color = colors$hex[match(linPhyloDF$lin, colors$lin)])
linPhylo_1 <- phylo4d(linPhylo, linTreeCols)

linTreeNodeData <- data.frame(lin = rep(NA, nNodes(linPhylo_1)),
                              color = rep("#000000", nNodes(linPhylo_1)),
                              row.names = nodeId(linPhylo_1, "internal"))

nodeData(linPhylo_1) <- linTreeNodeData

linPhyloPlot <- ggtree(linPhylo_1, aes(color=I(color)), ladderize = F) + 
        geom_tiplab(size = 8)

linPhyloPlot

ggtree::rotate(tree_view = linPhyloPlot, node = MRCA(linDend, c("L2", "L3")))
ggtree::rotate(tree_view = linPhyloPlot, node = MRCA(linDend, c("A1", "A2", "A3", "A4", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L9")))

ggsave(paste(outDir, "linPhylo.pdf"), height = 10, width = 1.5)