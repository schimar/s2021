#! /usr/bin/python

# This script needs some info HERE!

# Usage: script/allel.py vars/<zarrPath>

import os
import random
random.seed(42)
import sys
import allel as al
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
#sns.set_style('ticks')
#sns.set_context('notebook')
# --------------------------------------------


zarrPath = sys.argv[1]
statsPath = os.path.join(zarrPath.strip('.zarr'), 'stats/al/')
figsPath = os.path.join(zarrPath.strip('.zarr'), 'figs/al/')

# load the data
variants = zarr.open_group(zarrPath, mode='r')

#NOTE: ids will already be read in common.smk, so no need to do this here (for now, we need to, though)
ids = pd.read_table('samples109.txt', sep='\t', index_col=False)
ids['id_nest'] = ids['id'] + '_' + ids['nest']
ids = ids.sort_values(by='nest')


nInds = ids.groupby(by= ['nest', 'pop']).count()['sample']
np.all(list(variants['samples']) == ids['id_nest'].values)

samples = list(variants['samples'])
subsIndex = [samples.index(s) for s in ids['id_nest']]
ids['subsIndex'] = subsIndex

##
gtvars = al.GenotypeArray(variants['calldata/GT'])

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





### pops stats

ac_pops_vars = gtvars.count_alleles_subpops(pops, max_allele=1)
segAll_vars = ac_pops_vars['all'].is_segregating()[:]
gtseg_vars = gtvars.compress(segAll_vars, axis=0)
nAltVars = gtseg_vars.to_n_alt()

gtseg_is_het = pd.DataFrame(gtseg_vars.is_het()).describe().transpose()
gtseg_is_het.index = ids['id_nest']
gtseg_is_het['het'] = gtseg_is_het['freq']/gtseg_is_het['count']
gtseg_is_het['hom'] = 1-gtseg_is_het['het']
gtseg_is_het.to_csv(os.path.join(statsPath, 'nSegHets.txt'), header=True, index=True, sep= '\t')



ac_nests_vars = gtvars.count_alleles_subpops(nests, max_allele=1)

# nSeg per pop/nest as DataFrame:
popseg = dict()
for pop in ac_pops_vars.keys():
    popseg[pop] = ac_pops_vars[pop].count_segregating()
pd.Series(popseg, index=popseg.keys()).to_csv(os.path.join(statsPath, 'nSegAlleles.pop.txt'), header= False, index=True, sep= '\t')


nestseg = dict()
for nest in ac_nests_vars.keys():
    nestseg[nest] = ac_nests_vars[nest].count_segregating()
pd.Series(nestseg, index=nestseg.keys()).to_csv(os.path.join(statsPath, 'nSegAlleles.nest.txt'), header= False, index=True, sep= '\t')




# pairwise distance matrix

dvar = al.pairwise_distance(gtvars.to_n_alt(), metric= 'cityblock')


# heatmap with dendrogram

condensedDvar = scipy.spatial.distance.squareform(dvar)

n2col = dict(zip(ids['nest'].unique(), sns.color_palette()))
rowCols = np.array(ids['nest'].map(n2col))

cDdf = pd.DataFrame(condensedDvar, index= ids['nest'], columns=ids['id'])
g = sns.clustermap(cDdf, row_colors= rowCols, cmap= 'jet')
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 12)     #ha= 'right', rotation= 40
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 12)





#f = plt.figure()
#plt.plot(range(10), range(10), "o")
#plt.show()
#
#f.savefig("foo.pdf", bbox_inches='tight')
