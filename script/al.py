#! /usr/bin/python

# This script needs some info HERE!

# Usage: script/allel.py vars/<zarrPath>


import random
#random.seed(42)
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

variants = zarr.open_group(zarrPath, mode='r')

#ids are already read in common.smk, so no need to do this here
# ids['nest']

gtvars = al.GenotypeArray(variants['calldata/GT'])


subpops = {
    'all': list(range(len(ids))),
    'A': ids[ids['pop'] == 'A'].index.tolist(),
    'N': ids[ids['pop'] == 'N'].index.tolist(),
    'S': ids[ids['pop'] == 'S'].index.tolist(),
}
###

ac_subpops_vars = gtvars.count_alleles_subpops(subpops, max_allele=1)

segAll_vars = ac_subpops_vars['all'].is_segregating()[:]

gtseg_vars = gtvars.compress(segAll_vars, axis=0)

nAltVars = gtseg_vars.to_n_alt()

gtseg_is_het = pd.DataFrame(gtseg_vars.is_het()).describe().transpose()

gtseg_is_het.index = ids['id']

gtseg_is_het['het'] = gtseg_is_het['freq']/gtseg_is_het['count']
gtseg_is_het['hom'] = 1-gtseg_is_het['het']


print(gtseq_is_het.head())






