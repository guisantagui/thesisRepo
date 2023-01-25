"""This script takes tables (.csv) generated in the deletion analysis and creates 
a .txt file with the deletion percentage of each gnumber and places it in it's 
corresponding folder"""

import pandas as pd
import os
from os import listdir
from os.path import join

delsPath = '/storage/PGO/results/mtb/deletions_allGenes/'
mapsPath = '/storage/PGO/data/mtb/mappings/v1/'

delFiles = sorted(listdir('/storage/PGO/results/mtb/deletions_allGenes/'))
delFiles = delFiles[1:]


for f in delFiles:
    delsDF = pd.read_csv(os.path.join(delsPath, f))
    print("Starting to do files of separated gnumbers of " + f + ".")
    for g in delsDF.columns[1:]:
        pathOut = os.path.join(mapsPath, g[0:3], g[3:5], g[5])
        pathOut = pathOut + '/' + g + '_dels.txt'
        gNumDels = delsDF[g]
        gNumDels = pd.concat([delsDF[delsDF.columns[0]], gNumDels], axis = 1)
        gNumDels.to_csv(pathOut, header=None, index=None, sep='\t', mode='a')
        print("Saved " + g + " deletions at " + pathOut + ".")
        
print("Done!")