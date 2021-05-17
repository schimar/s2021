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
import math
import zarr
#import scipy.parse
import scipy.spatial
#%matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
sns.set_style('whitegrid')

#from sklearn.decomposition import PCA
#from sklearn.preprocessing import StandardScaler

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



    #zarrPath = sys.argv[1]
    zarrPath = 'vars/taSubInDel.zarr/'


    zarrname = zarrPath.strip('.zarr/')
    # create folders

    varname = zarrname.split('/')[1]

    statsP = os.path.join(zarrname, 'stats/ngsRelate/')
    figsP = os.path.join(zarrname, 'figs/ngsRelate/')
    gemmasP = os.path.join(zarrname, 'stats/gemma')

    folderList = [statsP, figsP]
    for folder in folderList:
        if not os.path.exists(folder):
            os.makedirs(folder)


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

    tnests = { key: nests[key] for key in nests if key not in ['all'] }


    rel = pd.read_table('vars/taSubInDel/stats/ngsRelate/stats.txt', sep= '\t')
    rel['na'] = rel['ida'].str.split('_').str[1]
    rel['nb'] = rel['idb'].str.split('_').str[1]

    nestDict = {}
    for nest in tnests.keys():
        neststats = rel[(rel['na'] == nest) & (rel['nb'] == nest)]
        nestDict[nest] = neststats[['rab', 'KING', 'R0', 'R1']]

    #print(nest, nestrab.mean(), nestrab.std(), effQnumPamilo(nestrab.mean()))

    nIndNest = [len(vals) for vals in tnests.values()]


    #Effective queen number after Pamilo (1991)
    #f = (3-r)/(3r)
    #https://www.jstor.org/stable/pdf/2462158.pdf
    #
    #Effective queen number after Pedersen & Boomsma (2001)
    #neg = 3/4*r
    #https://onlinelibrary.wiley.com/doi/full/10.1046/j.1420-9101.1999.00109.x

    def effQnumPamilo(r):
        f = (3-r)/(3*r)
        return f

    def effQnumPedBoo(r):
        f = 3/(4*r)
        return f


#pd.DataFrame([(nest, df['rab'].mean(), df['rab'].std(), effQnumPamilo(df['rab'].mean() - (1.96*df['rab'].std())), effQnumPamilo(df['rab'].mean()), effQnumPamilo(df['rab'].mean() + (1.96*df['rab'].std()))) for nest, df in nestDict.items()])

    nq = pd.DataFrame([(nest, df['rab'].mean(), df['rab'].std(), effQnumPamilo(df['rab'].mean()), effQnumPedBoo(df['rab'].mean())) for nest, df in nestDict.items()], columns= ['nest', 'meanRel', 'sdRel', 'nQueenPamilo', 'nQueenPedBoo'])
    nq.to_csv(os.path.join(statsP, 'relNqueen.txt'), header=True, index=False, sep= '\t')

    nqid = pd.Series(nq.nQueenPamilo.values, index= nq.nest).to_dict()
    ids['enq'] = ids['nest'].map(nqdi)
    ids['enq'].to_csv(os.path.join(gemmasP, 'enq.pheno'), header= False, index= False, sep= ' ')



    nest_cols = {
        'A2': sns.color_palette()[0],
        'A5': sns.color_palette()[1],
        'A6': sns.color_palette()[2],
        'N1': sns.color_palette()[3],
        'N4': sns.color_palette()[4],
        'N6': sns.color_palette()[5],
        'S1': sns.color_palette()[6],
        'S2': sns.color_palette()[7],
        'S5': sns.color_palette()[8]
    }




    ncols = 3  # how many plots per row
    nrows = math.ceil(len(nestDict) / ncols)  # how many rows of plots
    plt.figure(figsize=(10, 4))  # change the figure size as needed
    for i, (k, v) in enumerate(nestDict.items(), 1):
        plt.subplot(nrows, ncols, i)
        p = sns.distplot(v['rab'])#scatterplot(data=v, x='rab', y='KING', palette=nest_cols)
        #p.legend_.remove()
        plt.title(k)
    plt.tight_layout()
    plt.savefig(os.path.join(figsP, 'distplot.png'), bbox_inches='tight')



    # ---------------------------------------------------------

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
                        m.iloc[j,i] = np.float(pair[param])
                    except:
                        continue
        return(m)   #.transpose())

    m = matrixify(rel, ids['id_nest'], 'rab')


    fig = plt.figure(figsize=(16, 14))
    ax = sns.heatmap(m, cmap= 'coolwarm', annot= False)
    ax.set_title("Relatedness matrix")
    fig.tight_layout()
    fig.savefig(os.path.join(figsP, 'relmat.png'), bbox_inches='tight')


    #sns.scatterplot(rel['R0'], rel['rab'])








#[j for i in range(100) if i > 10 for j in range(i) if j < 20]


