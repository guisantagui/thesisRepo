#!/bin/bash
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=modsPipe
#SBATCH --mem=500G
#SBATCH --time=15-00:00:00
#Define sdout path
#SBATCH --output=/home/guisana/scripts/modsPipe/output_modsPipe.txt
#Define sderr path
#SBATCH --error=/home/guisana/scripts/modsPipe/error_modsPipe.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH --partition=long

# General variables for whole pipeline
#######################################################################################################
                                                                                                      #
modDir="/home/guisana/scripts/modsPipe/data/curMods/" # ---------------------------------------------># Path to curated model directory
delMatDir="/home/guisana/scripts/modsPipe/data/delMats/" # ------------------------------------------># Path to toRem matrixes directory
baseModDir="/home/guisana/scripts/modsPipe/data/models/" # ------------------------------------------># Path to directory where base reconstruction is stored
fluxSampSize=1000 # ---------------------------------------------------------------------------------># Sample size for flux sampling
output="/storage/PGO/results/mtb/gsmmSims/" # -------------------------------------------------------># General output of the pipeline
rxnDetsFile="/home/guisana/scripts/modsPipe/data/rxnDets/iEK1011_2.0_rxnDets.csv" # -----------------># iEK1011 2.0 reaction details file
mwAlpha=0.01 # --------------------------------------------------------------------------------------># Significance threshold for Mann Whitney U test

# Run simulations with deletion/all potentially deletereous SNPs models:
#######################################################################################################
module load python/3.8

# /usr/bin/time python3 gsmmSims.py --modDir $modDir --whatMods delsAllSNPs --delMatDir $delMatDir --baseModDir $baseModDir --fluxSampSize $fluxSampSize --output $output

# Run simulations with deletion/stop gained SNPs models:
#######################################################################################################
# /usr/bin/time python3 gsmmSims.py --modDir $modDir --whatMods delsSGSNPs --delMatDir $delMatDir --baseModDir $baseModDir --fluxSampSize $fluxSampSize --output $output

# Run statistical analyses on the sampled fluxes distributions:
#######################################################################################################

module load R

/usr/bin/time Rscript sampAnal.R --resDir $output --whatMod delsAllSNPs --rxnDetsFile $rxnDetsFile --medium mi7H9OADCCholMed --mwAlpha $mwAlpha --whatSampSize highGrwth
#/usr/bin/time Rscript sampAnal.R --resDir $output --whatMod delsSGSNPs --rxnDetsFile $rxnDetsFile --medium mi7H9OADCCholMed --mwAlpha $mwAlpha