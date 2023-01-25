#!/bin/bash
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=snpsPipe
#SBATCH --mem=500G
#SBATCH --time=15-00:00:00
#Define sdout path
#SBATCH --output=/home/guisana/scripts/snpsPipe_filtGenomes/output_snpsPipe.txt
#Define sderr path
#SBATCH --error=/home/guisana/scripts/snpsPipe_filtGenomes/error_snpsPipe.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH --partition=long


module load R

# General variables for whole pipeline
#######################################################################################################

btwSNPFile="/storage/PGO/data/mtb/wrk_dataset/mutations/positions_btw_DR/btw.snps.snpeff.annot" # ---># Path to btwSNPs
intraSNPDir="/storage/PGO/data/mtb/wrk_dataset/mutations/positions_intra_DR/" # ---------------------># Directory to intra SNP files.
outDir="/storage/PGO/results/mtb/snpsPipe_bestQual/" # ----------------------------------------------># Output directory
annDir="/storage/PGO/data/mtb/annotations/" #--------------------------------------------------------># Directory where annotation files are allocated. genes.gff, 
                                                                                                      # h37RvNCBIAnnot.gff3, Mycobacterium_tuberculosis_H37Rv_gff_v1.gff 
                                                                                                      # and Mycobacterium_tuberculosis_H37Rv_gff_v2.gff need to be in there. 
thrshldStpGa=0.00 # ---------------------------------------------------------------------------------># Threshold of proportion of whole lineages having a stop 
                                                                                                      # gained SNP for including it in the .txt files with the 
                                                                                                      # lists of genes with stop gain SNPs per lineage.
excludePos='/storage/PGO/scripts/pipeline/mtb/exclude_without_PPE.txt' # ----------------------------># Path to file with positions to exclude from SNP analysis because they
                                                                                                      # are problematic
top=40 # --------------------------------------------------------------------------------------------># Number of top variables to show in RF plots
thrshldToRem=0.9 # ---------------------------------------------------------------------------------># Threshold proportion of a SNP within a lineage to appear in toRem files
refStrain="L4" # ------------------------------------------------------------------------------------># Lineage to what the reference strain belongs to. Here L4, as H37Rv belongs
                                                                                                      # to it. 
iEKGenes="/storage/PGO/data/mtb/GSMMs/iEK1011_2.0/iEK1011_2.0_genes.txt" # --------------------------># Path to iEK1011 genes file. 
gNumDir="/storage/PGO/results/mtb/deletions_allGenes/allDelsFilt/linLists/" # -----------------------># Path to directory where files with G numbers per lineage are stored
mapsDir="/storage/PGO/data/mtb/mappings/v1/" # ------------------------------------------------------># Path to the directory of the output of the SNP calling/deletion pipeline

# Create output directory if it doesn't exist
if [ ! -d $outDir ]; then mkdir $outDir; fi

# Parsing and PROVEAN
#######################################################################################################

# Parse SNP files 
# Rscript SNPs.R --btwSNPFile $btwSNPFile --intraSNPDir $intraSNPDir --output $outDir --annotDir $annDir --thrshldStpGa $thrshldStpGa --filtNonSyn --filtStopGa

# Run provean
module load python/3.8

inpProv=$outDir"allSNPsNonSyn.csv"
anno=$outDir"allGenesMTBC.txt"
outProv=$outDir"provean/"

if [ ! -d $outProv ]; then mkdir $outProv; fi

/usr/bin/time python3 proveanRun.py $inpProv --output $outProv --refSeq $refSeq --anno $anno

# Create a matrix with potentially deletereous SNPs, by mixing missense SNPs significative in provean and stop gained SNPs

Rscript parseDeletSNPs.R --output $outDir

# Unsupervised analysis
#######################################################################################################

# Do Correspondence Analysis of the potentially deletereous SNPs data. 
outCA=$outDir"CA/"
inpCA=$outDir"deletSNPsMat.csv"

if [ ! -d $outCA ]; then mkdir $outCA; fi

Rscript CA_snps.R $inpCA --output $outCA

# Supervised analysis
#######################################################################################################

# Run supervised analysis (Random Forest). 
inpRF=$outCA"deletSNPsMat_sums.csv"
outRFDir=$outDir"RF/"
if [ ! -d $outRFDir ]; then mkdir $outRFDir; fi
/usr/bin/time Rscript RF.R $inpRF -o $outRFDir -s up --top $top #-l A1,A2,A3,A4,L5,L6,L9, 

# Plot results of Random Forest. 
snpFile=$outDir"deletSNPsMat.csv"

/usr/bin/time Rscript RF_plots.R --rfDir $outRFDir --provean $outProv --annot $annDir --snpPath $snpFile --top $top --output $outRFDir


#######################################################################################################
# Get SNPs per sampled isolates                                                                       #
#######################################################################################################

# Parse VCF files to get a matrix of potentially deletereous SNPs by individual strains. 

outIndivSNPs=$outDir"indivSNPs/"
if [ ! -d $outIndivSNPs ]; then mkdir $outIndivSNPs; fi
/usr/bin/time python3 parseIndivSNPs.py --gNumDir $gNumDir --mapsDir $mapsDir --btwSNPFile $btwSNPFile --intraSNPDir $intraSNPDir --deletSNPFile $snpFile --output $outIndivSNPs

# From individual G-number SNP file, sum the number of times a gene is affected by a potentially deletereous SNP in each isolate
/usr/bin/time Rscript sumIndivSNPs.R --indSNPDir $outIndivSNPs

# Unsupervised analysis
#######################################################################################################

# Do Correspondence Analysis of the potentially deletereous SNPs data (with proportion file generated
# with sampled genomes). 
outCA_indiv=$outDir"CA_indiv/"
inpCA_indiv=$outDir"deletSNPsMat_indiv.csv"

if [ ! -d $outCA_indiv ]; then mkdir $outCA_indiv; fi

Rscript CA_snps.R $inpCA_indiv --output $outCA_indiv # --multSampSize

# Supervised analysis
#######################################################################################################

# Run random forest on summed individual SNPs
indRFInp=$outIndivSNPs"deletSNPsMat_byStrain_sums.csv"
rfByStrainOut=$outRFDir"byStrain/"
if [ ! -d $rfByStrainOut ]; then mkdir $rfByStrainOut; fi

# /usr/bin/time Rscript RF.R $indRFInp --output $rfByStrainOut --plots --top $top --metricRF Kappa


# Plot results of individual strain SNPs random forest:
varImpFile=$rfByStrainOut"varImpdeletSNPs_byStrainDF.csv"
provFile=$outProv"proveanSign.csv"

/usr/bin/time Rscript RF_indivSNPsHM.R --deletSNPsMat $snpFile --varImpFile $varImpFile --top $top --provFile $provFile --output $rfByStrainOut

# Create a dataframe of the proportion of individual sampled isolates with each SNP 
/usr/bin/time Rscript getPropSNPsIndivStrains.R --indivSNPDir $outIndivSNPs --output $outDir

# Create lists of genes to remove from models
#######################################################################################################

# Filter deletereous SNPs dataset to obtain lists of genes to remove from iEK1011 2.0 GSMMs (stop gain SNPs and PROVEAN significative, 
# present in more than a given proportion within each lineage. 

toRemDir=$outDir"toRem/"

if [ ! -d $toRemDir ]; then mkdir $toRemDir; fi


/usr/bin/time Rscript filtDeletSNPs.R --snpFile $snpFile --thrshld $thrshldToRem --ref $refStrain --annoDir $annDir --output $toRemDir --iEKGenes $iEKGenes


# Create lists of genes to remove from models, based on proportions obtained from the sampled strains.
#######################################################################################################

# Filter deletereous SNPs dataset to obtain lists of genes to remove from iEK1011 2.0 GSMMs (stop gain SNPs and PROVEAN significative, 
# present in more than a given proportion within each lineage. 

snpFileIndiv=$outDir"deletSNPsMat_indiv.csv"
toRemDirIndiv=$outDir"toRemIndiv/"

if [ ! -d $toRemDirIndiv ]; then mkdir $toRemDirIndiv; fi


/usr/bin/time Rscript filtDeletSNPs.R --snpFile $snpFileIndiv --thrshld $thrshldToRem --ref $refStrain --annoDir $annDir --output $toRemDirIndiv --iEKGenes $iEKGenes