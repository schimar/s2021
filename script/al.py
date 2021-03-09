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


zarrname = zarrPath.strip('Bypop.zarr')

statsP = os.path.join(zarrname, 'stats/al/')
figsP = os.path.join(zarrname, 'figs/al/')
pcafP = os.path.join(figsP, 'pca/')      # pca figs
pcasP = os.path.join(statsP, 'pca/')     # pca stats

# create folders
folderList = [statsP, figsP, pcasP, pcafP]
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
gtseg_is_het.to_csv(os.path.join(statsP, 'nSegHets.txt'), header=True, index=True, sep= '\t')



ac_nests_vars = gtvars.count_alleles_subpops(nests, max_allele=1)

# nSeg per pop/nest as DataFrame:
popseg = dict()
for pop in ac_pops_vars.keys():
    popseg[pop] = ac_pops_vars[pop].count_segregating()
pd.Series(popseg, index=popseg.keys()).to_csv(os.path.join(statsP, 'nSegAlleles.pop.txt'), header= False, index=True, sep= '\t')


nestseg = dict()
for nest in ac_nests_vars.keys():
    nestseg[nest] = ac_nests_vars[nest].count_segregating()
pd.Series(nestseg, index=nestseg.keys()).to_csv(os.path.join(statsP, 'nSegAlleles.nest.txt'), header= False, index=True, sep= '\t')




#############   pairwise distance matrix   #############

dvar = al.pairwise_distance(gtvars.to_n_alt(), metric= 'cityblock')


# heatmap with dendrogram

condensedDvar = scipy.spatial.distance.squareform(dvar)

n2col = dict(zip(ids['nest'].unique(), sns.color_palette()))
rowCols = np.array(ids['nest'].map(n2col))

cDdf = pd.DataFrame(condensedDvar, index= ids['nest'], columns=ids['id'])
g = sns.clustermap(cDdf, row_colors= rowCols, cmap= 'jet')
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 12)     #ha= 'right', rotation= 40
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 12)

g.savefig(os.path.join(figsP, 'pwDist.pdf'), bbox_inches='tight')

# ---------------------------------------------------------------------

#############   PCA   #############

def plot_ld(gn, title, filetitle):
    m = al.rogers_huff_r(gn) ** 2
    ax = al.plot_pairwise_ld(m)
    ax.set_title(title)
    ax.figure.savefig(os.path.join(pcafP, filetitle), bbox_inches='tight')



plot_ld(nAltVars[:1000], 'Pairwise LD after random downsampling', 'ld_1000.pdf')



# random subsampling of loci
n = 100000  # number of SNPs to choose randomly
vidxVars = np.random.choice(nAltVars.shape[0], n, replace=False)
vidxVars.sort()

gnrVars = nAltVars.take(vidxVars, axis=0)

plot_ld(gnrVars[:1000], 'Pairwise LD after random downsampling.', 'ld_100k_rand.pdf')


segScafs = variants['variants/CHROM'][:][segAll_vars]
segBP = variants['variants/POS'][:][segAll_vars]
segVars = pd.DataFrame({'bp': segScafs, 'scaf': segBP})


# LD pruning
def ld_prune(gn, varInf, size, step, threshold=.1, n_iter=1):
    for i in range(n_iter):
        loc_unlinked = al.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
        gn = gn.compress(loc_unlinked, axis=0)
        varInf = varInf[loc_unlinked]
    return (gn, varInf)



gnuVars, vars_ldPrd = ld_prune(nAltVars, segVars, size=50, step=20, threshold=.1, n_iter=4)


plot_ld(gnuVars[:1000], 'Pairwise LD after LD pruning.', filetitle= 'ld_prune.pdf')


#populations = ids['pops'].unique()
pop_cols = {
    'A': sns.color_palette()[0],
    'N': sns.color_palette()[3],
    'S': sns.color_palette()[6],
}

#nests = ids['nest'].unique()
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


############

def plot_pca_coords(coords, model, pc1, pc2, ax, pops, pcols):
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for pop in pops.unique():
        flt = (pops.values == pop)
        ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color= pcols[pop],
                label=pop, markersize=6, mec='k', mew=.5)
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))


def fig_pca(coords, model, title, pops, pcols, filetitle):
    # plot coords for PCs 1 vs 2, 3 vs 4
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, 0, 1, ax, pops= pops, pcols= pcols)
    ax = fig.add_subplot(1, 2, 2)
    plot_pca_coords(coords, model, 2, 3, ax, pops= pops, pcols= pcols)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(pcafP, filetitle), bbox_inches='tight')


# nests

# PCA using SVD - LD-pruned data (168k loci)
coords1var, model1var = al.pca(gnuVars, n_components=10, scaler='patterson')


fig_pca(coords1var, model1var, 'LD-pruned PCA', pops = ids['nest'], pcols= nest_cols, filetitle= 'pca_LDprune.pdf')


######
# pca without LD pruning (random subset of 100000 loci)
coords2var, model2var = al.pca(gnrVars, n_components=10, scaler='patterson')

# pops
#fig_pca(coords2var, model2var, 'Figure 4.2. Conventional PCA', pops = ids['pops'], pcols= pop_cols)

# nests
fig_pca(coords2var, model2var, 'Conventional PCA', pops = ids['nest'], pcols= nest_cols, filetitle= 'pca_100k_rand.pdf')


# now for the full set (gtseg_vars)
coords2allVars, model2allVars = al.pca(nAltVars, n_components=10, scaler='patterson')


# pops
#fig_pca(coords2allVars, model2allVars, 'Conventional PCA without LD pruning all variants', pops = ids['pops'], pcols= pop_cols)

# nests
fig_pca(coords2allVars, model2allVars, 'Conventional PCA without LD pruning', pops= ids['nest'], pcols= nest_cols, filetitle= 'pca_all.pdf')


# pca with LD pruning, without Patterson's scaling
coords3vars, model3vars = al.pca(gnuVars, n_components=10, scaler=None)


# pops
#fig_pca(coords3vars, model3vars, 'Conventional PCA LD-pruned variants without variance scaling', pops = ids['pops'], pcols= pop_cols)

# nests
fig_pca(coords3vars, model3vars, 'Conventional PCA LD-pruned variants without variance scaling.', pops = ids['nest'], pcols= nest_cols, filetitle= 'pca_LDprune_noPatterson.pdf')


# randomized PCA with LD pruning
coords5vars, model5vars = al.randomized_pca(gnuVars, n_components=10, scaler='patterson')


# pops
#fig_pca(coords5vars, model5vars, 'Randomized PCA', pops= ids['pops'], pcols= pop_cols)

# nests
fig_pca(coords5vars, model5vars, 'Randomized PCA LD-pruned variants', pops= ids['nest'], pcols= nest_cols, filetitle= 'pca_LDprune_rand.pdf')


def plotHeatPCs(coords, ids, filetitle, PCs=4):
    df = pd.DataFrame(coords[:,0:PCs].T, columns=ids, index= range(1,PCs+1))
    plt.subplots(figsize= (20,5))
    ax = sns.heatmap(df, cmap= 'plasma')
    ax.set_ylabel('PC')
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
    ax.set_xticklabels(ids, rotation= 40, ha= 'right', fontsize= 8)
    plt.tight_layout()
    plt.savefig(os.path.join(pcafP, filetitle), bbox_inches='tight')


plotHeatPCs(coords1var, ids['nest'], PCs=5, filetitle= 'pca_LDprune_Heat.pdf')
plotHeatPCs(coords2allVars, ids['nest'], PCs=5, filetitle= 'pca_all_Heat.pdf')

## get the Eigen values for PCAs

# for all (segreg.) vars (n = 1,319,775)
segScafs = bbvars['variants/CHROM'][:][segAll_vars]
segBP = bbvars['variants/POS'][:][segAll_vars]
segVars = pd.DataFrame({'bp': segScafs, 'scaf': segBP})
## for gemma: segVars.to_csv('bbvars_f1byPopSegregate.scafbp', sep= ' ', index= False, header= False)



# write first 4 eigen values to file:
pd.DataFrame([ x[:4] for x in coords1var ]).to_csv("bbvars_f1byPopLDpruned.eigen", sep= ' ', index= False, header= False)

pd.DataFrame([ x[:4] for x in coords2allVars ]).to_csv("bbvars_f1byPop.eigen", sep= ' ', index= False, header= False)


#f = plt.figure()
#plt.plot(range(10), range(10), "o")
#plt.show()
#
#f.savefig("foo.pdf", bbox_inches='tight')
