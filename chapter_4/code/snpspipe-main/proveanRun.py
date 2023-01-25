import argparse
import subprocess as sp
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import os
import copy as cp
import re

parser = argparse.ArgumentParser(description='Run PROVEAN using as input a dataframe of missense SNPs.')
parser.add_argument('input', metavar='N', type=str,
                    help='Path to input file')

parser.add_argument('--anno', 
	                  type=str,
                    help='Path to annotation file', 
                    default='/storage/PGO/data/mtb/annotations/genesTXT4Pipeline/allGenesMTBC.txt')
parser.add_argument('--refSeq', 
	                  type=str,
                    help='Path to reference sequence fasta file', 
                    default='/storage/PGO/data/mtb/wrk_dataset/fastas/MTB_ancestor_reference.fasta')
parser.add_argument('--output', 
	                  type=str,
                    help='Path to output directory', 
                    default='/storage/PGO/results/mtb/snpsPipe/provean/')


args = parser.parse_args()

# Directory stuff
##########################################################################################################
missSnpPath = args.input
annotPath = args.anno
refPath = args.refSeq
outDir = args.output

# Load and parse data
##########################################################################################################

# Parse the missense_variant SNP csv file. Obtain dataframe of unique genes and variants. 
SNPs = pd.read_csv(missSnpPath, index_col = 0, low_memory=False)
# Remove possible NaNs
SNPs = SNPs.loc[[pd.isna(x) != True for x in SNPs['LOCUS']]]
SNPs = SNPs.reset_index()

SNPs_rdx = SNPs[['POS', 'REF', 'ALT', 'LIN', 'ANN', 'LOCUS', 'AA_one', 'propSNP']]


missUniqSNPLst = list(np.unique((SNPs_rdx['LOCUS'] + '.' + SNPs_rdx['AA_one']).tolist()))
missUniqSNPLst.sort()
missUniqSNPDF = pd.DataFrame({'GENE':[x.split('.')[0] for x in missUniqSNPLst], 'AA':[x.split('.')[1] for x in missUniqSNPLst]})

# Load the position files generated from ancestral GFF file. 
posNames = ['start', 'end', 'strand', 'gene']
genePos = pd.read_csv(annotPath, sep = '\t', names = posNames, index_col = None)

# Load the ancestral FASTA as a string. 
refGenome = SeqIO.parse(refPath, 'fasta')
for record in refGenome:
    refSeq = str(record.seq)

# For each gene with missense SNPs, get AA sequence from annotation file and reference fasta file, 
# as well as .var file with the AA substitutions in the gene, and run PROVEAN. Parse scores in single 
# dataframe. 
##########################################################################################################

# Iterate over missUniqSNPDF df, write a fasta for each protein with a SNP and compute PROVEAN score 
# with the AA change.  

# This function parses output generated from PROVEAN into a pd.DataFrame
def parsePROVEAN(gene):
    provRaw = pd.read_csv('%s%s.txt'%(outDir, gene), sep = '\t', index_col = None) 
    provMut = list(provRaw.index)
    if '# VARIATION' in provMut:
        resIdx = provMut.index('# VARIATION') + 1
        var = provMut[resIdx:]
        provScore = provRaw['## PROVEAN v1.1 output ##'].tolist()
        score = provScore[resIdx:]
        outDF = pd.DataFrame({'gene': list(np.repeat(gene, len(score))), 'var': var, 'score': score})
    # This else statement is for evaluating what genes have annotation errors. 
    else:
        annotErrors.append(gene)
        outDF = pd.DataFrame()
    return outDF

      
notInPosFile = []
annotErrors = []

genesWSNPs = np.unique(missUniqSNPDF['GENE'])
provFiles = os.listdir(outDir)

print("jajajajaj")

allGenesProvDF = pd.DataFrame()
for i in range(0, len(genesWSNPs)):
    gene = genesWSNPs[i]
    if gene not in genePos['gene'].tolist() and 'c' not in gene:
        gene = gene + 'c'
    if gene in genePos['gene'].tolist():
        var = missUniqSNPDF.loc[missUniqSNPDF['GENE'] == gene, 'AA']
        varFName = outDir + gene + '.var'
        start = genePos.loc[genePos['gene'] == gene, 'start'].tolist()[0] -1
        end = genePos.loc[genePos['gene'] == gene, 'end'].tolist()[0]
        strand = genePos.loc[genePos['gene'] == gene, 'strand'].tolist()[0]
        geneSeq = refSeq[start:end]
        if strand == 'F':
            seqTr = Seq(geneSeq).translate()
        elif strand == 'R':
            seqTr = Seq(geneSeq).reverse_complement().translate()
        #print(gene)
        #print(geneSeq)
        #print(seqTr)
        seqTrRec = SeqRecord(seqTr, id = '', name = '', description = gene)
        seqTrOutPath = outDir + gene + '.fa'
        if '%s.txt'%gene not in provFiles:
            var.to_csv(varFName, header = False, index = False)
            SeqIO.write(seqTrRec, seqTrOutPath, 'fasta')
            cmd = 'provean.sh -q %s%s.fa -v %s%s.var --num_threads 12 > %s%s.txt'%(outDir, gene, outDir, gene, outDir, gene)
            print('Running PROVEAN with command line %s.'%cmd)
            sp.call(cmd, shell = True)
            # Remove .fa and .var files once PROVEAN has run. 
            sp.call('rm '+ seqTrOutPath, shell = True)
            sp.call('rm '+ varFName, shell = True)
        else:
            print('%s.txt already exists.'%gene)
        # Parse PROVEAN output:
        provDF = parsePROVEAN(gene)
        if provDF.shape[0] > 0:
            allGenesProvDF = pd.concat([allGenesProvDF, provDF])
    else:
        notInPosFile.append(gene)


print('Genes that have SNP annotation errors: ')
print(annotErrors)

print('Genes that are not in annotation file: ')
print(notInPosFile)

# Create a list of lineages affected by each SNP, removing the genes that have annotation errors, if there are, 
# and add it as a column to the output dataframe


# function for doing the string matching
match = lambda a, b: [ b.index(x) if x in b else None for x in a ]


allGenesProvDF = allGenesProvDF.reset_index(drop = True)

# Create lists of concatenated gene and its variable to make the match of the lineages.
#geneVarProvDF = (allGenesProvDF['gene'] + '_' + allGenesProvDF['var']).tolist()
#geneVarSNPsDF = (SNPs['GENE'] + '_' + SNPs['AA']).tolist()

#outLins = SNPs.iloc[match(geneVarProvDF, geneVarSNPsDF)]['LIN'].tolist()

#if len(annotErrors) > 0:
#    outLins = SNPs.loc[[x not in annotErrors for x in SNPs['GENE'].tolist()]]['LIN'].tolist()
#else:
#    outLins = SNPs['LIN'].tolist()

#print(len(outLins))
print(allGenesProvDF.shape)
print(allGenesProvDF.head)

#allGenesProvDF.to_csv('%sproveanScores.txt'%outDir, sep = '\t')

#allGenesProvDF['lin'] = outLins

#allGenesProvDF.to_csv('%sproveanScores.txt'%outDir, sep = '\t')

print('Saving PROVEAN parsed results in %sproveanScores.txt'%outDir)
allGenesProvDF.to_csv('%sproveanScores.txt'%outDir, sep = '\t')
print('proveanScores.txt saved at %s!'%outDir)
 
# Add PROVEAN scores to the non synonymous SNP dataset and write the new dataset to 
# an new CSV file
print('Adding PROVEAN scores to non synonymous SNPs dataset...')
provMutLst = list(allGenesProvDF['gene'] + '.' + allGenesProvDF['var'])
provScores = []
for i in range(0, SNPs_rdx.shape[0]):
    mut = SNPs_rdx['LOCUS'][i] + '.' + SNPs_rdx['AA_one'][i]
    if(mut in provMutLst):
        score = allGenesProvDF['score'].loc[[x == mut for x in provMutLst]]
    else:
        score = None
    provScores.append(score)
    
SNPs_rdxProv = cp.copy(SNPs_rdx)  

SNPs_rdxProv['PROVEAN'] = provScores

SNPs_W_provName = re.sub('.csv', '_prov.csv', os.path.basename(missSnpPath))

SNPs_rdxProv.to_csv('%s%s'%(outDir, SNPs_W_provName))
print('%s saved at %s!'%(SNPs_W_provName, outDir))

print('Done!')