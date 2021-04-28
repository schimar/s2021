#! /usr/bin/python

# This script needs some info HERE!

# Usage: script/allel.py vars/<zarrPath>

import os
import random
random.seed(42)
import sys
import allel as al
from itertools import combinations
import numpy as np
np.random.seed(42)
import pandas as pd
import zarr
#import scipy.parse
import scipy.spatial
#%matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
sns.set_style('whitegrid')

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

#sns.set_style('ticks')
#sns.set_context('notebook')
# --------------

# Plotting pretty figures and avoid blurry images
#%config InlineBackend.figure_format = 'retina'
# No need to include %matplotlib inline magic command. These things come built-in now.

# Ignore warnings
import warnings
warnings.filterwarnings('ignore')

# --------------------------------------------



# --------------------------------------------

if __name__ == "__main__":



    zarrPath = sys.argv[1]


    zarrname = zarrPath.strip('.zarr/')
    # create folders

    varname = zarrname.split('/')[1]


#    statsP = os.path.join(zarrname, 'stats/al/')
#    figsP = os.path.join(zarrname, 'figs/al/')
#    pcafP = os.path.join(figsP, 'pca/')      # pca figs
#    pcasP = os.path.join(statsP, 'pca/')     # pca stats
#    varpcafP = os.path.join(pcafP, 'varPca/')
#    varpcasP = os.path.join(pcasP, 'varPca/')
#    varDenseP = os.path.join(figsP, 'varDense/')
#    hetfP = os.path.join(figsP, 'hets/')
#    sfsP = os.path.join(figsP, 'sfs/')
#    fstsP = os.path.join(zarrname, 'stats/al/fst/')
#    fstfP = os.path.join(zarrname, 'figs/al/fst/')
#    fstsPnests = os.path.join(zarrname, 'stats/al/fst/nests/')
#    fstsPpops = os.path.join(zarrname, 'stats/al/fst/pops/')
#    gemmasP = os.path.join(zarrname, 'stats/gemma/')
#    selsP = os.path.join(zarrname, 'stats/al/sel/')
#    selfP = os.path.join(zarrname, 'figs/al/sel/')
#    vcftlsP = os.path.join(zarrname, 'stats/vcftools')


#    folderList = [statsP, figsP, pcasP, pcafP, varpcafP, varpcasP, varDenseP, hetfP, sfsP, fstsPnests, fstsPpops, fstfP, selsP, selfP, gemmasP, vcftlsP]
#    for folder in folderList:
#        if not os.path.exists(folder):
#            os.makedirs(folder)



    # load the data
variants = zarr.open_group(zarrPath, mode='r')

#NOTE: ids will already be read in common.smk, so no need to do this here (for now, we need to, though)
ids = pd.read_table('samples109.txt', sep='\t', index_col=False)
ids['id_nest'] = ids['id'] + '_' + ids['nest']
ids = ids.sort_values(by='nest')


nInds = ids.groupby(by= ['nest', 'pop']).count()['sample']
#np.all(list(variants['samples']) == ids['id_nest'].values)

samples = list(variants['samples'])
subsIndex = [samples.index(s) for s in ids['id_nest']]
ids['subsIndex'] = subsIndex
ids.sort_values(by=['subsIndex'], inplace= True)
##
# check if the ids are in the same order
# np.all(list(variants['samples'][:]) == ids['id_nest'])

# add variable 'agg' (with 0 = started aggression, 1 = reacted peacefully, 2 = reacted aggressively)
ids['agg'] = 0
ids['agg'][ids['started_aggression'] == 1] = 0
ids['agg'][ids['reacted_peacefully'] == 1] = 1
ids['agg'][ids['reacted_aggressively'] == 1] = 2

# create genotype array
gtvars = al.GenotypeArray(variants['calldata/GT'])

# get locus ids (for getting scafbp in resp. sets of segregating sites) NOTE: see getScafBp()
scaf = variants['variants/CHROM'][:]
bp = variants['variants/POS'][:]
scafbp = scaf + ':' + list(map(str, bp))

idx = al.UniqueIndex(np.arange(0, len(scafbp)))




# define the groups (nests & pops)
pops = {
    'all': list(range(len(ids))),
    'A': ids[ids['pop'] == 'A'].index.tolist(),
    'N': ids[ids['pop'] == 'N'].index.tolist(),
    'S': ids[ids['pop'] == 'S'].index.tolist(),
}

nests = {
    'all': list(range(len(ids))),
    'A2': ids[ids['nest'] == 'A2'].index.tolist(),
    'A5': ids[ids['nest'] == 'A5'].index.tolist(),
    'A6': ids[ids['nest'] == 'A6'].index.tolist(),
    'N1': ids[ids['nest'] == 'N1'].index.tolist(),
    'N4': ids[ids['nest'] == 'N4'].index.tolist(),
    'N6': ids[ids['nest'] == 'N6'].index.tolist(),
    'S1': ids[ids['nest'] == 'S1'].index.tolist(),
    'S2': ids[ids['nest'] == 'S2'].index.tolist(),
    'S5': ids[ids['nest'] == 'S5'].index.tolist()
}



rel = pd.read_table('vars/taSubInDel/stats/ngsRelate/stats.txt', sep= '\t')


def matrixify(df, ids, param):
    m = np.zeros(shape=(109,109))
    m = pd.DataFrame(m, index= ids, columns= ids)
    for i, ind1 in enumerate(m.index):
        for j, ind2 in enumerate(m.columns):
            if i == j:
                m.iloc[i,j] = np.NaN   #continue
            else:
                pair = df[(df['a'] == i) & (df['b'] == j)]
                try:
                    np.float(pair[param])
                    m.iloc[i,j] = np.float(pair[param])
                except:
                    continue
    return(m.transpose())

m = matrixify(rel, ids['id_nest'], 'rab')



sns.heatmap(m, cmap= 'coolwarm', annot= False)


sns.scatterplot(rel['theta'], rel['rab'])










