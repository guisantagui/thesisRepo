import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser(description='Filter provean dataframe to keep just significative SNPs.')

parser.add_argument('--provDFFile', 
	                  type=str,
                    help='Path to provean DF file, with provean scores and each one of the lineages.', 
                    default='/storage/PGO/results/mtb/snpsPipe/provean/allSNPsNonSyn_prov.csv')
                    
parser.add_argument('--output', 
	                  type=str,
                    help='Path to mappings directory, where file system of all vcf files for all G-numbers is stored', 
                    default='/storage/PGO/results/mtb/snpsPipe/provean/')

args = parser.parse_args()

provDFFile = args.provDFFile
outDir = args.output


provDF = pd.read_csv(provDFFile)


print('Number of rows in provean DF: %s.'%str(provDF.shape[0]))
print('Obtaining provean scores from provean strings located in scores column...')
scores = []
for i in range(0, len(provDF['PROVEAN'].tolist())):
    print(i)
    scoreStr = provDF['PROVEAN'].tolist()[i]
    if not pd.isna(scoreStr):
        score = float(re.sub('\nName:', '', scoreStr.split(' ')[4]))
    else:
        score = float('NaN')
    scores.append(score)
    
provDF['PROVEAN'] = pd.Series(scores)


print('Removing NaNs')
provDF = provDF.loc[[not pd.isna(x) for x in provDF['PROVEAN']]]

provDF = provDF.reset_index(drop = True)

print('Filtering dataset to keep provean significative SNPs...')
provDF = provDF.loc[[x <= -2.5 for x in provDF['PROVEAN']]]


provDF = provDF.reset_index(drop = True)
provDF.to_csv('%sproveanSign.csv'%outDir)

print('proveanSign.csv saved at %s.'%outDir)