#!/bin/bash
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=delsPipe
#SBATCH --mem=500G
#SBATCH --time=15-00:00:00
#Define sdout path
#SBATCH --output=/home/guisana/scripts/delsPipe/output_delsPipe.txt
#Define sderr path
#SBATCH --error=/home/guisana/scripts/delsPipe/error_delsPipe.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH --partition=long


module load R

# General variables for whole pipeline
#######################################################################################################

input="/storage/PGO/results/mtb/deletions_allGenes/allDels.csv" # -----------------------------------># Input file for whole pipeline.
outDir="/storage/PGO/results/mtb/delsPipe_bestQual/" # -------------------------------------------------------># Output directory.
binThrshld=15 # -------------------------------------------------------------------------------------># Threshold for binarizing the percentage of ORF lost.
labkeyFile="/storage/PGO/data/mtb/genomesInfo/labkey_qual.csv" # ------------------------------------># File with genomes of highest quality that have Accession numbers
refLin="L4" # ---------------------------------------------------------------------------------------># Lineage used as reference for univariate test.
linDelThrshld=90 # ----------------------------------------------------------------------------------># Threshold of proportion of isolates within a lineage with a 
                                                                                                      # statistically significative ORF lost to keeping the ORF in the 
                                                                                                      # toRem files.
statTestSign="CS" # ---------------------------------------------------------------------------------># Statistical test for obtaining differentially deleted genes.
                                                                                                      # Possible values: -CS -->Chi Squared
                                                                                                      #                  -Fi -->Fisher test
                                                                                                      #                  -MW -->Mann Whitney U test
top=50 # --------------------------------------------------------------------------------------------># Number of top variables to show in outputs of Random Forest.
btwSNPFile="/storage/PGO/data/mtb/wrk_dataset/mutations/positions_btw_DR/btw.snps.snpeff.annot" # ---># Path to btw SNP file to add stop-gain SNPs to toRem files.
iEKGenes="/storage/PGO/data/mtb/GSMMs/iEK1011_2.0/iEK1011_2.0_genes.txt" # --------------------------># Path to iEK1011 genes file.
annot="/storage/PGO/data/mtb/annotations/genes.gff" # -----------------------------------------------># Path to annotation GFF file.
treeFile="/home/guisana/scripts/delsPipe_filtGenomes/" # --------------------------------------------># Path to tree file for doing the simplified phylogeny to show in plots.

# Create output directory if it doesn't exist
if [ ! -d $outDir ]; then mkdir $outDir; fi

# Filter allDels.csv to include just genomes in labkey file (with high quality and accession number available), and reset input variable to use the filtered file.
/usr/bin/time Rscript filtAllDels.R --allDels $input --labkeyFile $labkeyFile
input=$(dirname "$input")"/allDelsFilt/allDels.csv"

# Binarize deletion data, both full set and enzymatic set. 
/usr/bin/time Rscript binarizeDels.R $input --output $outDir --threshold $binThrshld --enzFilt

# Run PCA and HCA
/usr/bin/time Rscript delsUnsup.R $input --output $outDir --enzs --plots

# Run PCA of iEK1011 Genes
/usr/bin/time Rscript iEK1011GenesPCA.R $input --iEKGenes $iEKGenes --output $outDir --plots

# Run CA
/usr/bin/time Rscript CA.R $input --output $outDir --iEKGenes $iEKGenes --enzs --iEK --plots

# Run univariate tests
/usr/bin/time Rscript univAnalAllLinsVsRef.R $input --output $outDir --ref $refLin

# Obtain dataframes of proportion of deleted isolates within each lineage with each ORF deleted with a given percentage or more. 
univSignDir=$outDir"linVs"$refLin"DiffGenes/*"$statTestSign"Sign.txt"
        # Create outPropDir if it doesn't exist. 
outPropDir=$outDir"signDelProp/"
if [ ! -d $outPropDir ]; then mkdir $outPropDir; fi
for file in $univSignDir;
 do
  fileName=$(echo "$file");
  echo "Getting within lineage deletion proportion of genes in "$fileName"..."
  /usr/bin/time Rscript getDelPropLins.R "$fileName" -o $outPropDir -l $binThrshld --high 90 -p --delFile $input
done

# Filter the sets of statistically significative genes from past step to keep just the ones that are lost in the
# target lineage above a given proportion of the total number of isolates in it while present in more than that
# same threshold of the reference lineage. 
        # Create outToRemDir if it doesn't exist. 
outToRemDir=$outDir"toRem/"
if [ ! -d $outToRemDir ]; then mkdir $outToRemDir; fi
/usr/bin/time Rscript filtOvrDeldGns.R $outPropDir --output $outToRemDir --thrshld $linDelThrshld --ref $refLin --iEKGenes $iEKGenes --annotFile $annot

# Run supervised analysis (Random Forest). 
outRFDir=$outDir"RF/"
if [ ! -d $outRFDir ]; then mkdir $outRFDir; fi
/usr/bin/time Rscript RF.R $input -o $outRFDir -s cust -p --top $top #-l A1,A2,A3,A4,L5,L6,L9, 

# Get deletion proportion of the top X important genes in the RF classification. 
RF_filt_input=$(find $outRFDir -name '*TopGenes.txt')
/usr/bin/time Rscript getDelPropLins.R $RF_filt_input -o $outRFDir -l $binThrshld --high 90 -p --rowReord L8,L3,L2,L4,L7,L1,A3,A4,A2,L6,L9,A1,L5

# Do barplots of variable importance of the Random Forest analysis with variables colored according 
# to the RDs.  
/usr/bin/time Rscript varImpWRDs.R --allDels $input --output $outDir --annot $annot --top $top

# Plot simplified phylogeny.  
/usr/bin/time Rscript doRedPhylo.R --treeFile $treeFile --output $outDir