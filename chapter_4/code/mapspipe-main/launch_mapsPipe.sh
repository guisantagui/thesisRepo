#!/bin/bash
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=mapsPipe
#SBATCH --mem=50G
#SBATCH --time=15-00:00:00
#Define sdout path
#SBATCH --output=/home/guisana/scripts/mapsPipe/output_mapsPipe.txt
#Define sderr path
#SBATCH --error=/home/guisana/scripts/mapsPipe/error_mapsPipe.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH --partition=long


module load python/3.8

# Change these directories according to the path where the download repo is and 
# where mappings will be stored.
rootDir="/home/guisana/scripts/mapsPipe/"
mapDir="/storage/PGO/data/mtb/mappings/v1/"

gNumDir=$rootDir"data/"
gNumFiles=$gNumDir"*.txt"

# Obtain deletions

for file in $gNumFiles;
 do
  fileName=$(echo "$file");
  echo "Obtaining deletions of gNumbers included in "$fileName"..."
  /usr/bin/time python3 /storage/PGO/scripts/deletions_allGenes/deleted_epitopes.py -f $fileName -d $mapDir
done

