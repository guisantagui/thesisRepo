#########################################################################################################
#                                                                                                       #
# Run simulations with the GSMMs generated in jupyter notebooks.                                        #
#                                                                                                       #
#########################################################################################################
import os
import sys
sys.path.append(os.getcwd()) 

import argparse

import cobra
from cobra.flux_analysis import flux_variability_analysis
from cobra.sampling import sample
import pandas as pd
from os.path import join
import copy as cp
import numpy as np

parser = argparse.ArgumentParser(description='Run simulations with the GSMMs generated in jupyter notebooks.')

parser.add_argument('--modDir', 
	                  type=str,
                    help='Path to directory where lineage specific GSMMs are stored.', 
                    default='/home/guisana/scripts/modsPipe/data/curMods/')
                    
parser.add_argument('--whatMods', 
	                  type=str,
                    help='Choose if we are using the models accounting to deletions and stop gain SNPs or deletions and all potentially deletereous SNPs.', 
                    default='delsAllSNPs')
                    
parser.add_argument('--delMatDir', 
	                  type=str,
                    help='Directory where delMats are stored.',                    
                    default='/home/guisana/scripts/modsPipe/data/delMats/')
                    
parser.add_argument('--baseModDir', 
	                  type=str,
                    help='Directory where base reconstruction (iEK1011 2.0 and iEK1011) is stored.',                    
                    default='/home/guisana/scripts/modsPipe/data/models/')
                    
parser.add_argument('--fluxSampSize', 
	                  type=str,
                    help='Number of samples to get in flux sampling.',                    
                    default=5000)
                    
parser.add_argument('--output', 
	                  type=str,
                    help='General output directory.', 
                    default='/storage/PGO/results/mtb/gsmmSims/')
                    
args = parser.parse_args()

#
# Directory stuff
##########################################################################################################

modDir = args.modDir
delMatDir = args.delMatDir
output = args.output
whatMods = args.whatMods
baseModDir = args.baseModDir
fluxSampSize = args.fluxSampSize

if not os.path.exists(output):
    print('%s created'%output)
    os.makedirs(output)

if whatMods == 'delsAllSNPs':
    mod_dir = modDir + 'delsAllSNPs/'
    res_dir = output + 'delsAllSNPs/'
    delMatFile = delMatDir + 'del_AllSNPs_mat_cur.csv'
elif whatMods == 'delsSGSNPs':
    mod_dir = modDir + 'delsSGSNPs/'
    res_dir = output + 'delsSGSNPs/'
    delMatFile = delMatDir + 'del_sgSNPs_mat_cur.csv'

FBA_dir = res_dir + 'FBA/'

# Create directory if it doesn't exist
sampling_dir = res_dir + 'sampling/'

if not os.path.exists(sampling_dir):
    print('%s created'%sampling_dir)
    os.makedirs(sampling_dir)

    
#
# Load data and models for obtaining different media: 
##########################################################################################################   

# Load base reconstruction
base_mod = cobra.io.load_json_model(join(baseModDir, "iEK1011_2.0.json"))

# Create 7H9 OADC medium, obtained from S8 file matlab scripts from iEK1011 2.0 paper
mi7H9_OADC_Med = {'EX_glu__L_e':1,'EX_cu2_e':1000,'EX_btn_e':1,'EX_pydxn_e':1,'EX_ca2_e':1000, 'EX_mg2_e':1000, 'EX_h_e':1000, 'EX_k_e':1000, 'EX_nh4_e':10, 'EX_h2o_e':1000,'EX_pi_e':1, 'EX_cl_e':1000, 'EX_o2_e':20, 'EX_na1_e':1000, 'EX_so4_e':1000,'EX_cit_e':1,'EX_fe3_e':5, 'EX_glyc_e':1, 'EX_glc__D_e':1,'EX_ocdca_e':1}

mi7H9_OADCchol_Med = cp.copy(mi7H9_OADC_Med)
mi7H9_OADCchol_Med['EX_chsterol_e'] = 1

# Load different iEK1011 models (in different mediums) to retrieve mediums 
inVivoMed = cobra.io.load_json_model(baseModDir + 'iEK1011_inVivo_media.json').medium
deJesuMed = cobra.io.load_json_model(baseModDir + 'iEK1011_deJesusEssen_media.json').medium
drugTsMed = cobra.io.load_json_model(baseModDir + 'iEK1011_drugTesting_media.json').medium
grifEsMed = cobra.io.load_json_model(baseModDir + 'iEK1011_griffinEssen_media.json').medium
mi7H10Med = cobra.io.load_json_model(baseModDir + 'iEK1011_m7H10_media.json').medium

# Exchange reactions have different names in old model than in new, so change the keys to 
# make mediums compatible with models. 

for m in list(inVivoMed.keys()):
    if 'L' in m or 'D' in m:
        k = m[0:(len(m)-1)] + '_' + m[len(m)-1] + '_e'
    else:
        k = m + '_e'
    inVivoMed[k] = inVivoMed.pop(m)

for m in list(deJesuMed.keys()):
    if 'L' in m or 'D' in m:
        k = m[0:(len(m)-1)] + '_' + m[len(m)-1] + '_e'
    else:
        k = m + '_e'
    deJesuMed[k] = deJesuMed.pop(m)

for m in list(drugTsMed.keys()):
    if 'L' in m or 'D' in m:
        k = m[0:(len(m)-1)] + '_' + m[len(m)-1] + '_e'
    else:
        k = m + '_e'
    drugTsMed[k] = drugTsMed.pop(m)

for m in list(grifEsMed.keys()):
    if 'L' in m or 'D' in m:
        k = m[0:(len(m)-1)] + '_' + m[len(m)-1] + '_e'
    else:
        k = m + '_e'
    grifEsMed[k] = grifEsMed.pop(m)

for m in list(mi7H10Med.keys()):
    if 'L' in m or 'D' in m:
        k = m[0:(len(m)-1)] + '_' + m[len(m)-1] + '_e'
    else:
        k = m + '_e'
    mi7H10Med[k] = mi7H10Med.pop(m)
    
# Load delMat: 
delMat = pd.read_csv(delMatFile, index_col = 0)

#
# Functions needed in script: 
##########################################################################################################
def sampleFluxes(modLst, outDir, size = 100, outCsv = False, medium = mi7H10Med):
    lins = [x.id for x in defModLst]
    # sampFluxDict = {}
    sampFiles = os.listdir(outDir)
    for i in range(0, len(lins)):
        lin = lins[i]
        fileName = lin + '_samp.csv'
        if fileName not in sampFiles:
            print('Sampling %s model...'%lin)
            mod = modLst[i].copy()
            mod.medium = medium
            sampFlux = sample(mod, size)
        else: 
            # print('%s already sampled, loading %s'%(lin, fileName))
            print('%s already sampled, in %s'%(lin, fileName))
            # sampFlux = pd.read_csv(outDir + fileName, index_col = 0)
        # sampFluxDict[lin] = sampFlux
        if outCsv == True and fileName not in sampFiles:
            sampFlux.to_csv(outDir + lin + '_samp.csv')
            print('%s_samp.csv saved at %s.'%(lin, outDir))
    # return sampFluxDict

#
# ANALYSIS: 
##########################################################################################################  

# Read models
print('Loading %s models...'%whatMods)
L1 = cobra.io.load_json_model(join(mod_dir, "L1_cur.json"))
L2 = cobra.io.load_json_model(join(mod_dir, "L2_cur.json"))
L3 = cobra.io.load_json_model(join(mod_dir, "L3_cur.json"))
L4 = cobra.io.load_json_model(join(mod_dir, "L4_cur.json"))
L5 = cobra.io.load_json_model(join(mod_dir, "L5_cur.json"))
L6 = cobra.io.load_json_model(join(mod_dir, "L6_cur.json"))
L7 = cobra.io.load_json_model(join(mod_dir, "L7_cur.json"))
L8 = cobra.io.load_json_model(join(mod_dir, "L8_cur.json"))
L9 = cobra.io.load_json_model(join(mod_dir, "L9_cur.json"))
A1 = cobra.io.load_json_model(join(mod_dir, "A1_cur.json"))
A2 = cobra.io.load_json_model(join(mod_dir, "A2_cur.json"))
A3 = cobra.io.load_json_model(join(mod_dir, "A3_cur.json"))
A4 = cobra.io.load_json_model(join(mod_dir, "A4_cur.json"))

# Create a list of models

defModLst = [L1, L2, L3, L4, L5, L6, L7, L8, L9, A1, A2, A3, A4]

# Generate rxns removed files:
iEK1011_rxns = [x.id for x in base_mod.reactions]
rmRxnsTab = pd.DataFrame(index = iEK1011_rxns, columns = list(delMat.columns))
allRemRxns = []
for mod in defModLst:
    lin = mod.id
    modRxns = [x.id for x in mod.reactions]
    remRxnsLin = [x for x in iEK1011_rxns if x not in modRxns]
    allRemRxns = [*allRemRxns, *remRxnsLin]
    remRxnList = []
    for idx in list(rmRxnsTab.index):
        if idx in remRxnsLin:
            remRxnList.append(0.0)
        else:
            remRxnList.append(1.0)
    rmRxnsTab[lin] = remRxnList
allRemRxns = list(np.unique(allRemRxns))

allRemRxns = pd.Series(allRemRxns)
allRemRxns.to_csv('%srxnsRem_%s.txt'%(output, whatMods), sep = '\t', header = False, index = False)
print('rxnsRem_%s.txt saved at %s'%(whatMods, output))

rmRxnsTab = rmRxnsTab.loc[[sum(rmRxnsTab.loc[x] == 1) != rmRxnsTab.shape[1] for x in rmRxnsTab.index]]
rmRxnsTab.to_csv('%srxnsRem_tab_%s.csv'%(output, whatMods))
print('rxnsRem_tab_%s.csv saved at %s'%(whatMods, output))

#
# Flux sampling: 
########################################################################################################## 

# Do sampling in Middlebrook 7H9 OADC chol medium
print('Doing flux sampling with %s samples in Middlebrock 7H9 OADC chol medium...'%str(fluxSampSize))
samp_mi7H9OADCCholDir = sampling_dir + 'mi7H9OADCCholMed/n%s/'%str(fluxSampSize)
if not os.path.exists(samp_mi7H9OADCCholDir):
    print('%s created'%samp_mi7H9OADCCholDir)
    os.mkdir(samp_mi7H9OADCCholDir)

sampFluxDict_mi7H10Med = sampleFluxes(defModLst, size = fluxSampSize, outCsv = True, outDir = samp_mi7H9OADCCholDir, medium = mi7H9_OADCchol_Med)

# Do sampling in Middlebrook 7H10 medium
print('Doing flux sampling with %s samples in Middlebrock 7H10 medium...'%str(fluxSampSize))
samp_mi7H10Dir = sampling_dir + 'mi7H10Med/n%s/'%str(fluxSampSize)
if not os.path.exists(samp_mi7H10Dir):
    print('%s created'%samp_mi7H10Dir)
    os.mkdir(samp_mi7H10Dir)

sampFluxDict_mi7H10Med = sampleFluxes(defModLst, size = fluxSampSize, outCsv = True, outDir = samp_mi7H10Dir, medium = mi7H10Med)

# Do sampling in in vivo medium
print('Doing flux sampling with %s samples in in vivo medium...'%str(fluxSampSize))
samp_inVivoDir = sampling_dir + 'inVivoMed/n%s/'%str(fluxSampSize)
if not os.path.exists(samp_inVivoDir):
    print('%s created'%samp_inVivoDir)
    os.mkdir(samp_inVivoDir)

sampFluxDict_inVivoMed = sampleFluxes(defModLst, size = fluxSampSize, outCsv = True, outDir = samp_inVivoDir, medium = inVivoMed)

# Do sampling in Middlebrock 7H9 OADC medium
print('Doing flux sampling with %s samples in Middlebrock 7H9 OADC medium...'%str(fluxSampSize))
samp_7H9OADCDir = sampling_dir + 'mi7H9OADCMed/n%s/'%str(fluxSampSize)
if not os.path.exists(samp_7H9OADCDir):
    print('%s created'%samp_7H9OADCDir)
    os.mkdir(samp_7H9OADCDir)

sampFluxDict_mi7H9OADCMed = sampleFluxes(defModLst, size = fluxSampSize, outCsv = True, outDir = samp_7H9OADCDir, medium = mi7H9_OADC_Med)