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

nInds = ids.groupby(by= ['nest', 'pop']).count()['sample']
np.all(list(variants['samples']) == ids['id'].values)

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










### nests stats
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




#f = plt.figure()
#plt.plot(range(10), range(10), "o")
#plt.show()
#
#f.savefig("foo.pdf", bbox_inches='tight')
