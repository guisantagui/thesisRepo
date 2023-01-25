##########################################################################################################
#                                                                                                        #
# Builds, from the potentially deletereous SNP matrix, the same matrix but with the presence/absence     #
# of the SNPs in each one of the clinical isolates, by looking in the VCF files generated during SNP     #
# calling.                                                                                               #
#                                                                                                        #
##########################################################################################################

import argparse
import os
# import subprocess as sp
import sys
import re
import pandas as pd
import io
import copy as cp
import numpy as np
import statistics as st
import random as rand

# reqs = sp.check_output([sys.executable, '-m', 'pip', 'freeze'])
# installed_packages = [r.decode().split('==')[0] for r in reqs.split()]


parser = argparse.ArgumentParser(description='Parse SNP proportion file into a dataframe of presence/absence in individual isolates.')

parser.add_argument('--gNumDir', 
	                  type=str,
                    help='Path to directory where .txt files of the G-NUMBER lists of isolates within each lineage that were used to get SNP proportions.', 
                    default='/storage/PGO/data/mtb/wrk_dataset/gnumbers/')
                    
parser.add_argument('--mapsDir', 
	                  type=str,
                    help='Path to mappings directory, where file system of all vcf files for all G-numbers is stored', 
                    default='/storage/PGO/data/mtb/mappings/v1/')
                    
parser.add_argument('--btwSNPFile', 
	                  type=str,
                    help='Path to btw SNP file', 
                    default='/storage/PGO/data/mtb/wrk_dataset/mutations/positions_btw_DR/btw.snps.snpeff.annot')
                    
parser.add_argument('--intraSNPDir', 
	                  type=str,
                    help='Path to btw SNP file', 
                    default='/storage/PGO/data/mtb/wrk_dataset/mutations/positions_intra_DR/')
                    
parser.add_argument('--deletSNPFile', 
	                  type=str,
                    help='Deletereous SNP matrix CSV file', 
                    default='/storage/PGO/results/mtb/snpsPipe/deletSNPsMat.csv')
                    
parser.add_argument('--output', 
	                  type=str,
                    help='Output directory', 
                    default='/storage/PGO/results/mtb/snpsPipe/indivSNPs/')


args = parser.parse_args()


#
# Directory stuff
##########################################################################################################
gNumDir = args.gNumDir
mapsDir = args.mapsDir
intraSNPDir = args.intraSNPDir
btwSNPFile = args.btwSNPFile
deletSNPsFile = args.deletSNPFile
outDir = args.output

#
# Functions needed in script: 
##########################################################################################################

# Parse vcf file into pandas DF
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


# Match function
match = lambda a, b: [ b.index(x) if x in b else None for x in a ]

# Translate three letter amino acid code to one letter
def transAA(three_change):
    import Bio
    from Bio import Data
    from Bio.Data import IUPACData
    if '*' in three_change:
        fst = three_change[0:3]
        num = three_change[3:(len(three_change)-1)]
        snd = '*'
        sndTr = '*'
    else:
        fst = three_change[0:3]
        snd = three_change[-3:]
        num = three_change[3:(len(three_change)-3)]
        sndTr = Bio.Data.IUPACData.protein_letters_3to1[snd]
    fstTr = Bio.Data.IUPACData.protein_letters_3to1[fst]
    tr = fstTr + num + sndTr
    return tr

# Parse intra SNPs lineage file to have only missense_variant and stop_gained SNPs, 
# ensuring that different variants that in a lineage appear in the same position
# occupy different rows in the output dataframe
def parseIntraSNPs(filePath):
    # This function is for obtaining the changes that happen in the same position in separate 
    # rows in a dataframe,
    def getCommaSplitDF(commaDF):
        commaSplitDF = pd.DataFrame(columns=commaDF.columns)
        for i in range(0, commaDF.shape[0]):
            splitDF = pd.DataFrame(columns=commaDF.columns)
            splitAnn = commaDF['ANN'][i].split(',')
            splitAlt = commaDF['ALT'][i].split(',')
            splitAA3 = commaDF['AA_three'][i].split(',')
            for j in range(0, len(splitAnn)):
                changeSeries = commaDF.loc[i, :]
                changeSeries['ANN'] = splitAnn[j]
                changeSeries['ALT'] = splitAlt[j]
                changeSeries['AA_three'] = splitAA3[j]
                changeSeries = changeSeries.to_frame().T
                splitDF = splitDF.append(changeSeries, ignore_index = True)
            splitDF = splitDF[[x in ['missense_variant', 'stop_gained'] for x in splitDF['ANN']]]
            splitDF.reset_index(drop = True)
            commaSplitDF = commaSplitDF.append(splitDF, ignore_index = True)
        return commaSplitDF
    # Load data, tidy and change names
    intraSNPs = pd.read_csv(filePath, sep = '\t', index_col = False)
    intraSNPs = intraSNPs[['Position', 'Gene', 'REF', 'HOMO_SNP', 'annotation', 'aa change']]
    intraSNPs = intraSNPs.rename({'Position':'POS', 'Gene':'GENE', 'HOMO_SNP':'ALT', 'annotation':'ANN', 'aa change':'AA_three'}, axis = 1)
    intraSNPs = intraSNPs [[not pd.isna(x) for x in intraSNPs['ANN']]]
    intraSNPs['ANN'] = [re.sub(' ', '', x) for x in intraSNPs['ANN']]
    intraSNPs['ANN'] = [re.sub('stopgain', 'stop_gained', x) for x in intraSNPs['ANN']]
    # Create dataframes of missense_variants and stop_gained variants
    miss = intraSNPs[['missense_variant' in x for x in intraSNPs['ANN']]]
    stop = intraSNPs[['stop_gained' in x for x in intraSNPs['ANN']]]
    # Obtain subDataframes of the SNPs that are in the same position. ALT, AA_three and 
    # ANN attributes are all in same row separated by commas. We want each one in a separate
    # row.
    missComma = miss[[',' in x for x in miss['ANN']]]
    missComma = missComma.reset_index(drop = True)
    stopComma = stop[[',' in x for x in stop['ANN']]]
    stopComma = stopComma.reset_index(drop = True)
    # Remove the rows with commas from missense and stop dataframes
    miss = miss[[',' not in x for x in miss['ANN']]].reset_index(drop = True)
    stop = stop[[',' not in x for x in stop['ANN']]].reset_index(drop = True)
    # Obtain the dataframes of the SNPs that are in the same position, 
    # but this time with each one in a separated row.
    missSplit = getCommaSplitDF(missComma)
    stopSplit = getCommaSplitDF(stopComma)
    # Merge the 4 dataframes
    parsedSNPs = pd.concat([miss, stop, missSplit, stopSplit]).reset_index(drop = True)
    # Remove rows with NaNs from AA_three column, as it might be some because of the input file.
    # Remove also rows with ',', as there are some.
    parsedSNPs = parsedSNPs[[not pd.isna(x) for x in parsedSNPs['AA_three']]].reset_index(drop = True)
    parsedSNPs = parsedSNPs[[',' not in x for x in parsedSNPs['AA_three']]].reset_index(drop = True)
    # Reorder dataframe according to 'POS' and reset the index
    parsedSNPs = parsedSNPs.iloc[list(np.argsort(parsedSNPs['POS'].tolist()))].reset_index(drop = True)
    # Remove p. from the AA change
    parsedSNPs['AA_three'] = [re.sub('p\.', '', x) for x in parsedSNPs['AA_three'].tolist()]
    # Translate AA_three to AA_one and create the code for locating the SNP in VCF file
    parsedSNPs['AA_one'] = [transAA(x) for x in parsedSNPs['AA_three']]
    parsedSNPs['CODE'] = parsedSNPs['REF'].tolist() + pd.Series([str(x) for x in parsedSNPs['POS']]) + parsedSNPs['ALT'].tolist()
    return parsedSNPs

#
# Load data and parse the data: 
##########################################################################################################
print('Loading deletereous SNPs matrix %s from %s...'%(os.path.basename(deletSNPsFile), os.path.dirname(deletSNPsFile)))
deletSNPs = pd.read_csv(deletSNPsFile, index_col = 0)
# btw SNPs
btwSNPs = pd.read_csv(btwSNPFile, sep = '\t')
btwSNPs = btwSNPs[['POS', 'REF', 'ALT', 'LIN', 'Unnamed: 8', 'Unnamed: 9', 'Unnamed: 10']]
btwSNPs = btwSNPs.rename({'Unnamed: 8':'ANN', 'Unnamed: 9':'GENE', 'Unnamed: 10':'AA_three'}, axis = 1)
btwSNPs = btwSNPs[['POS', 'GENE', 'REF', 'ALT', 'LIN', 'ANN', 'AA_three']]
# Filter to keep just stop_gained and missense_variant, as potentially deletereous snps only account for 
# these classes, via the boolean lists and zip() function.
missBool = (btwSNPs['ANN'] == 'missense_variant').tolist()
stopGaBool = (btwSNPs['ANN'] == 'stop_gained').tolist()
btwSNPs = btwSNPs.loc[[a or b for a, b in zip(missBool, stopGaBool)]]
btwSNPs = btwSNPs.reset_index(drop = True)
btwSNPs['AA_three'] = [re.sub('p\.', '', x) for x in btwSNPs['AA_three'].tolist()]

#
# Build the matrix:
##########################################################################################################

# Iterate from vcf files of each strain and build the deletereous SNP dataframe by strain.

gNumFiles = [x for x in os.listdir(gNumDir) if '.txt' in x]
gNumFiles.sort()

print('Parsing through individual isolates VCF files for building potentially deletereous SNP matrix...')

outFiles = [x for x in os.listdir(outDir) if '.csv' in x]
outName = outDir + 'deletSNPsMat_byStrain.csv'


# Get sample size of number of strains for having the same number of animal and human associated lineages
    
lins = [re.sub('nr.txt', '', f) for f in gNumFiles]
linLens = [len(pd.read_csv('%s%s'%(gNumDir, f), header = None)[0].tolist()) for f in gNumFiles]
medLen = st.median(linLens)

toSampFromBig = round(cp.copy(medLen)/2)
toSampFromSma = round(toSampFromBig*9/4)
print('Number of samples from each animal associated lineages: %s'%str(toSampFromSma))
print('Number of samples from each human associated lineages: %s'%str(toSampFromBig))

# Look if there is a file already generated that might have some lineages already parsed, to not start 
# from zero. 

if len(outFiles) > 0:
    outFiles = [x for x in outFiles if x != 'deletSNPsMat_byStrain.csv']
    outFiles = [x for x in outFiles if x != 'deletSNPsMat_byStrain_sums.csv']
    print(outFiles)
    alreadyDone = [x.split('_')[len(x.split('_')) - 1] for x in outFiles]
    alreadyDone = [re.sub('.csv', '', x) for x in alreadyDone]
    for lin in alreadyDone:
        print('%s lineage G-numbers already parsed into deletSNPsMat_byStrain_%s.csv'%(lin, lin))
    gNumFiles = [x for x in gNumFiles if re.sub('nr.txt', '', x) not in alreadyDone]
    print('Reading %s from %s...'%(os.path.basename(outName), os.path.dirname(outName)))
    delSNPsByStrain = pd.read_csv(outName, index_col = 0)
    print(gNumFiles)
else:
    delSNPsByStrain = pd.DataFrame(columns = list(deletSNPs.columns))

# Shut down warning of chained assignment to avoid see the warning when "missense_variant,synonymous_variant"
# annotated SNPs ALT is splitted: ALT is of the type A,C, so we want to change it to A to avoid errors when 
# mapping to VCF files. 
  
pd.options.mode.chained_assignment = None

# Iterate over the each G number list file, generate a SNPs in lin file by merging btw SNPs and intra 
# SNPs, sample the G numbers of each lineage, access the VCF files and check what of the potentially
# deleterious SNPs appear in the VCFs of each strain. Use the SNP_lin file to determine the AA_code 
# correspondent to the NT code to add it to the VCF, and then check what of the potentially deleterious
# SNPs appear in the VCF file (the potentially deleterious SNPs are in AA_code, as they were determined
# with PROVEAN. Finally write a CSV file with the SNPs in the individual strains. 
for file in gNumFiles:
    delSNPsByStrain_byLin = pd.DataFrame(columns = list(deletSNPs.columns))
    lin = re.sub('nr.txt', '',file)
    print('Adding SNPs of lineage %s clinical isolates...'%lin)
    linGNums = pd.read_csv('%s%s'%(gNumDir, file), header = None)[0].tolist()
    # Sample the G-Nums
    if 'A' in lin:
        if toSampFromSma <= len(linGNums):
            rand.seed(10)
            linGNums = rand.sample(linGNums, toSampFromSma)
        elif toSampFromSma > len(linGNums):
            linGNums = rand.choices(linGNums, k=toSampFromSma)
    elif 'L' in lin:
        if toSampFromBig <= len(linGNums):
            rand.seed(9)
            linGNums = rand.sample(linGNums, toSampFromBig)
        elif toSampFromBig > len(linGNums):
            linGNums = rand.choices(linGNums, k=toSampFromBig)
    # Filt btw SNPs to keep the ones of the target lineage
    btwSNPs_lin = cp.copy(btwSNPs)
    btwSNPs_lin = btwSNPs_lin.loc[[lin in x for x in btwSNPs_lin['LIN']]]
    btwSNPs_lin = btwSNPs_lin.reset_index(drop = True)
    btwSNPs_lin = btwSNPs_lin[['POS', 'GENE', 'REF', 'ALT', 'ANN', 'AA_three']]
    btwSNPs_lin['AA_one'] = [transAA(x) for x in btwSNPs_lin['AA_three']]
    btwSNPs_lin['CODE'] = btwSNPs_lin['REF'].tolist() + pd.Series([str(x) for x in btwSNPs_lin['POS']]) + btwSNPs_lin['ALT'].tolist()
    # Load and parse intra SNPs of the target lineage
    intraSNPs_linPath = '%s%s_intralineage_DR.annot'%(intraSNPDir, lin)
    print('Parsing %s_intralineage_DR.annot...'%lin)
    intraSNPs_lin = parseIntraSNPs(intraSNPs_linPath)
    
    # intraSNPs_lin = pd.read_csv('%s%s_intralineage_DR.annot'%(intraSNPDir, lin), sep = '\t', index_col = False)
    # intraSNPs_lin = intraSNPs_lin[['Position', 'Gene', 'REF', 'HOMO_SNP', 'annotation', 'aa change']]
    # intraSNPs_lin = intraSNPs_lin.rename({'Position':'POS', 'Gene':'GENE', 'HOMO_SNP':'ALT', 'annotation':'ANN', 'aa change':'AA_three'}, axis = 1)
    # Create boolean lists for doing the zip comprehension list that allows to do the boolean operator in order to keep just stop_gain and missense 
    # SNPs
    # missBoolIntra = (intraSNPs_lin['ANN'] == 'missense_variant').tolist()
    # stopGaBoolIntra = (intraSNPs_lin['ANN'] == 'stop_gained').tolist()
    # missSynBoolIntra = (intraSNPs_lin['ANN'] == 'missense_variant,synonymous_variant').tolist()
    # intraSNPs_lin = intraSNPs_lin.loc[[a or b or c for a, b, c in zip(missBoolIntra, stopGaBoolIntra, missSynBoolIntra)]]
    # intraSNPs_lin = intraSNPs_lin.reset_index(drop = True)
    # intraSNPs_lin['AA_three'] = [re.sub('p\.', '', x) for x in intraSNPs_lin['AA_three'].tolist()]
    # Modify ALT and AA_three of the SNPs that are annotated as 'missense_variant,synonymous_variant' 
    # to remove the synonymous ALT and the synonymous AA (the second one).
    # for i in range(0, intraSNPs_lin.shape[0]):
    #     snpAnn = intraSNPs_lin['ANN'][i]
    #     if snpAnn == 'missense_variant,synonymous_variant':
    #         intraSNPs_lin['ALT'][i] = intraSNPs_lin['ALT'][i].split(',')[0]
    #         intraSNPs_lin['AA_three'][i] = intraSNPs_lin['AA_three'][i].split(',')[0]
    # Translate AA code to one letter and generate the nt change code we will use for mapping to VCF
    # intraSNPs_lin['AA_one'] = [transAA(x) for x in intraSNPs_lin['AA_three']]
    # intraSNPs_lin['CODE'] = intraSNPs_lin['REF'].tolist() + pd.Series([str(x) for x in intraSNPs_lin['POS']]) + intraSNPs_lin['ALT'].tolist()
    
    # Merge btw SNPs and intra SNPs
    snps_lin = pd.concat([btwSNPs_lin, intraSNPs_lin])
    snps_lin = snps_lin.reset_index(drop = True)
    snps_lin['GENE'] = [x.split('_')[len(x.split('_'))-1] for x in snps_lin['GENE'].tolist()]
    snps_lin['AA_CODE'] = snps_lin['GENE'] + '.' + snps_lin['AA_one']
    for gnum in linGNums:
        fst = gnum[0:3]
        snd = gnum[3:5]
        trd = gnum[5]
        pathToVCF =  '%s%s/%s/%s/%s.var.snp.vcf'%(mapsDir, fst, snd, trd, gnum)
        #print('Parsing %s file from lineage %s, in directory %s...'%(gnum, lin, os.path.dirname(pathToVCF)))
        vcf = read_vcf(pathToVCF)
        vcf_parsed = cp.copy(vcf)
        vcf_parsed = vcf_parsed[['POS', 'REF', 'ALT']]
        vcf_parsed['CODE'] = vcf_parsed['REF'].tolist() + pd.Series([str(x) for x in vcf_parsed['POS']]) + vcf_parsed['ALT'].tolist()
        # Filter to keep just the ones in the intra and btw files
        vcf_parsed = vcf_parsed.loc[[x in snps_lin['CODE'].tolist() for x in vcf_parsed['CODE']]]
        vcf_parsed = vcf_parsed.reset_index(drop = True)
        # Add the gene.AA_one snp code
        vcf_parsed['AA_CODE'] = snps_lin['AA_CODE'].loc[match(vcf_parsed['CODE'].tolist(), snps_lin['CODE'].tolist())].tolist()
        # Check what SNPs in delet SNPs DF are present in the vcf file 
        pattLst = [int(x in vcf_parsed['AA_CODE'].tolist()) for x in list(deletSNPs.columns)]
        #print(sum(pattLst))
        snpDict = dict(zip(list(deletSNPs.columns), pattLst))
        gNum_DF = pd.DataFrame(snpDict, index = ['%s_%s'%(gnum, lin)])
        delSNPsByStrain_byLin = pd.concat([delSNPsByStrain_byLin, gNum_DF])
        delSNPsByStrain = pd.concat([delSNPsByStrain, gNum_DF])
    delSNPsByStrain_byLin.to_csv('%sdeletSNPsMat_byStrain_%s.csv'%(outDir, lin))
    delSNPsByStrain.to_csv('%sdeletSNPsMat_byStrain.csv'%outDir)
    print('Potentially deletereous SNPs of clinical isolates of %s written to deletSNPsMat_byStrain_%s.csv in %s'%(lin, lin, outDir))
    print('Potentially deletereous SNPs of clinical isolates of %s added to deletSNPsMat_byStrain.csv in %s'%(lin, outDir))

# delSNPsByStrain.to_csv('%sdeletSNPsMat_byStrain.csv'%outDir)
# print('deletSNPsMat_byStrain.csv saved at %s.'%outDir)