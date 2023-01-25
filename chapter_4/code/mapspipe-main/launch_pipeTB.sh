#!/bin/bash
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=PipeTB_Garnatxa
#SBATCH --ntasks=10
#SBATCH --partition=short
#SBATCH --cpus-per-task=10
#SBATCH --mem=100gb
#SBATCH --time=1-00:00:0 
#Define sdout path
#SBATCH --output=/storage/PGO/PGOgit/prueba/output_sample_prueba_2.txt
#Define sderr path
#SBATCH --error=/storage/PGO/PGOgit/prueba/error_sample_prueba_2.txt


  

module load R

#comando para que podamos hacer el array con todos los bnumbers

FASTQ=$1


cat $FASTQ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 > /storage/PGO/PGOgit/prueba/input.txt 

echo "$SLURM_ARRAY_TASK_ID"


#LINEA PARA QUITAR ERROR DE JAVA QUALIMAPBAMQC http://qualimap.bioinfo.cipf.es/doc_html/faq.html#x11problem
unset DISPLAY 


#command line WITH array
python /storage/PGO/PGOgit/prueba/Pipeline_TB.py --git /storage/PGO/PGOgit/prueba/ --min-var-freq 0.10 --min-var-freq-for-hom 0.90 --read-length 20 --min-length 20 -p 6 -f /storage/PGO/PGOgit/prueba/input.txt -o /storage/PGO/PGOgit/prueba/map_ch/ --excluded-loci Locus_to_exclude_Mtb.txt 

#command line for console
#python /storage/PGO/PGOgit/prueba/Pipeline_TB.py --git /storage/PGO/PGOgit/prueba/ --min-var-freq 0.10 --min-var-freq-for-hom 0.90 --read-length 20 --min-length 20 -p 3 -f /storage/PGO/PGOgit/prueba/same_sample_prueba.txt -o /storage/PGO/PGOgit/prueba/ --excluded-loci Locus_to_exclude_Mtb.txt 


#ojo, dejar la ultima linea en blanco porque la borra
#NFIL=`wc -l same_sample_solo1.txt | cut -d ' ' -f1`

#sbatch -a 1-$NFIL  /storage/PGO/PGOgit/prueba/launch_pipeTB.sh /storage/PGO/PGOgit/prueba/same_sample_solo1.txt



