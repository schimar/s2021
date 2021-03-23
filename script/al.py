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

def plot_ld(gn, title, filename):
    m = al.rogers_huff_r(gn) ** 2
    ax = al.plot_pairwise_ld(m)
    ax.set_title(title)
    ax.figure.savefig(os.path.join(pcafP, filename), bbox_inches='tight')

def ld_prune(gn, varInf, size, step, threshold=.1, n_iter=1):
    for i in range(n_iter):
        loc_unlinked = al.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('LD pruning: iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
        gn = gn.compress(loc_unlinked, axis=0)
        varInf = varInf[loc_unlinked]
    return (gn, varInf)

def plot_pca_coords(coords, model, pc1, pc2, ax, pops, pcols):
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for pop in pops.unique():
        flt = (pops.values == pop)
        ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color= pcols[pop], label=pop, markersize=6, mec='k', mew=.5)
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))


def fig_pca(coords, model, title, pops, pcols, filename):
    # plot coords for PCs 1 vs 2, 3 vs 4
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, 0, 1, ax, pops= pops, pcols= pcols)
    ax = fig.add_subplot(1, 2, 2)
    plot_pca_coords(coords, model, 2, 3, ax, pops= pops, pcols= pcols)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(pcafP, filename), bbox_inches='tight')


def plotHeatPCs(coords, ids, filename, PCs=4):
    df = pd.DataFrame(coords[:,0:PCs].T, columns=ids, index= range(1,PCs+1))
    plt.subplots(figsize= (20,5))
    ax = sns.heatmap(df, cmap= 'plasma')
    ax.set_ylabel('PC')
    #ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
    #ax.set_xticklabels(ids, rotation= 40, ha= 'right', fontsize= 8)
    plt.tight_layout()
    plt.savefig(os.path.join(pcafP, filename), bbox_inches='tight')


def plot_windowed_variant_density(pos, window_size, filename, title=None):
    # setup windows
    bins = np.arange(0, pos.max(), window_size)
    # use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1])/2
    # compute variant density in each window
    h, _ = np.histogram(pos, bins=bins)
    y = h / window_size
    # plot
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    if title:
        ax.set_title(title)
    fig.savefig(os.path.join(varDenseP, filename), bbox_inches='tight')

def plot_variant_hist(f, filename, bins=30):
    x = variants[os.path.join('variants',f)][:]
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.despine(ax=ax, offset=10)
    ax.hist(x, bins=bins)
    ax.set_xlabel(f)
    ax.set_ylabel('No. variants')
    ax.set_title('Variant %s distribution' % f)
    fig.savefig(os.path.join(varDenseP, filename), bbox_inches='tight')


def plot_variant_hist_2d(f1, f2, downsample, filename):
    x = variants[os.path.join('variants',f1)][:][::downsample]
    y = variants[os.path.join('variants',f2)][:][::downsample]
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.despine(ax=ax, offset=10)
    ax.hexbin(x, y, gridsize=40)
    ax.set_xlabel(f1)
    ax.set_ylabel(f2)
    ax.set_title('Variant %s versus %s joint distribution' % (f1, f2))
    fig.savefig(os.path.join(varDenseP, filename), bbox_inches='tight')

def plotPropHets(propHets, ids, filename):
    plt.subplots(figsize= (20,5))
    ax = sns.barplot(np.arange(len(propHets)), propHets, hue= ids['nest'].values, dodge= False)
    ax.set_ylim([0,1])
    ax.set_xlabel('samples')
    ax.set_ylabel('proportion heterozygous')
    ax.set_title('proportion of heterozygous genotype calls')
    #ax.set_xticks(np.arange(len(propHets)))
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=True)
    ax.set_xticklabels(ids['id_nest'].values, rotation= 40, ha= 'right', fontsize= 8)
    ax.figure.savefig(os.path.join(hetfP, filename), bbox_inches='tight')

def plot_jsfs(ac, fname):
    tmpDict = { key: ac[key] for key in ac if key not in ['all'] }
    combs = list(combinations(tmpDict, 2))
    for popPair in combs:
        jsfs = al.stats.sf.joint_sfs(ac[popPair[0]][:,1], ac[popPair[1]][:,1])
        fig, ax = plt.subplots(figsize=(6,6))
        al.stats.sf.plot_joint_sfs(jsfs, ax=ax)
        ax.set_ylabel(' '.join(['Alternate allele count,', popPair[0] ]))
        ax.set_xlabel(' '.join(['Alternate allele count,', popPair[1] ]))
        fig.savefig(os.path.join(sfsP, '.'.join([fname, popPair[0], popPair[1], 'png' ])), bbox_inches='tight')

def plot_jsfs_nests(ac, fname):
    tmpDict = { key: ac[key] for key in ac if key not in ['all'] }
    for pop in ['A', 'N', 'S']:
        popDict = { k: v for k, v in tmpDict.items() if pop in k }
        nestCombs = list(combinations(popDict, 2))
        for nestPair in nestCombs:
            jsfs = al.stats.sf.joint_sfs(ac[nestPair[0]][:,1], ac[nestPair[1]][:,1])
            fig, ax = plt.subplots(figsize=(6,6))
            al.stats.sf.plot_joint_sfs(jsfs, ax=ax)
            ax.set_ylabel(' '.join(['Alternate allele count,', nestPair[0] ]))
            ax.set_xlabel(' '.join(['Alternate allele count,', nestPair[1] ]))
            fig.savefig(os.path.join(sfsP, '.'.join([fname, nestPair[0], nestPair[1], 'png' ])), bbox_inches='tight')


# -------------------------------------------------------- #
# W & C's fst with block-jackknife standard errors
def calcWCfst_bj_knife(ac, pairs, gtvars, idx1, idx2, blen):
    acu = al.AlleleCountsArray(ac[pairs[0]][:] + ac[pairs[1]][:])
    is_seg = acu.is_segregating() & (acu.max_allele() == 1)
    gtmp = gtvars.compress(is_seg, axis=0)
    segSitesPos = scafbp[getScafBp(idx, is_seg)]
    # Weir & Cockerham's
    fst, se, vb, vj = al.stats.fst.average_weir_cockerham_fst(gtseg_vars, subpops=[idx1, idx2], blen=blen, max_allele=1)
    return fst, se, vb, pairs, np.count_nonzero(is_seg), segSitesPos, is_seg

#calcWCfst_bj_knife(ac_nests_vars, ('A2', 'S1'), gtvars, nests['A2'], nests['S1'], blen=10000)
#calcWCfst_bj_knife(ac_pops_vars, ('A', 'S'), gtvars, pops['A'], pops['S'], blen=10000)

def fst_bj_knife_pair(ac, gtvars, poplvl, blen=10000, subsample= True):
    tmpDict = { key: ac[key] for key in ac if key not in ['all'] }
    combs = list(combinations(tmpDict, 2))
    res = list()
    resPerPair = dict()
    is_seg_dict = dict()
    for pair in combs:
        if subsample:
            minlen = min(len(poplvl[pair[0]]), len(poplvl[pair[1]]))
            idx1 = random.sample(poplvl[pair[0]], minlen)
            idx2 = random.sample(poplvl[pair[1]], minlen)
        else:
            idx1 = poplvl[pair[0]]
            idx2 = poplvl[pair[1]]
        wc, se, vb, pairs, nSegVars, segSitesPos, is_seg = calcWCfst_bj_knife(ac, pair, gtvars, idx1, idx2, blen)    # NOTE: add is_seg to output below???!!!
        res.append((pairs[0], len(idx1), pairs[1], len(idx2), nSegVars, wc, se))
        pwScafBp = pd.DataFrame({'scafbp': segSitesPos})
        resPerPair['_'.join([pair[0], pair[1]])] = pwScafBp
        is_seg_dict['_'.join([pair[0], pair[1]])] = is_seg

    return pd.DataFrame(res, columns= ('pair1', 'n1', 'pair2', 'n2', 'nSegAlleles', 'bjknife.wcFst', 'std.err.wcFst')), resPerPair, is_seg_dict



def plot_bjFst(df, fname):
    fig, ax = plt.subplots(figsize=(15,5))
    plt.errorbar(np.arange(1,df.shape[0]+1), df['bjknife.wcFst'], yerr=df['std.err.wcFst'], fmt='o', color='Black', elinewidth=2, capthick=1, errorevery=1, alpha=1, ms=4, capsize=5)
    plt.bar(np.arange(1,df.shape[0]+1), df['bjknife.wcFst'], tick_label=np.arange(1,df.shape[0]+1))
    plt.xlabel('pairs')
    plt.ylabel("$F_{ST}$")
    ax.set_xticklabels(df['pair1'] + ' - ' + df['pair2'], rotation= 40, ha= 'right', fontsize= 8)
    fig.savefig(os.path.join(fstfP, '.'.join([fname, 'bjFst_bar.pdf'])), bbox_inches='tight')




def getScafBp(varsIDX, selVars):
    if len(varsIDX) != len(selVars):
        raise ValueError('The two supplied files can not have different lengths')
    else:
        return varsIDX.astype('int64')[selVars]



def calcWCfst_per_site(ac, pairs, gtvars, idx1, idx2):
    acu = al.AlleleCountsArray(ac[pairs[0]][:] + ac[pairs[1]][:])
    is_seg = acu.is_segregating() & (acu.max_allele() == 1)
    gtmp = gtvars.compress(is_seg, axis=0)
    segSitesPos = scafbp[getScafBp(idx, is_seg)]
    # Weir & Cockerham's
    a, b, c = al.weir_cockerham_fst(gtmp, subpops=[ idx1, idx2 ], max_allele=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        snp_fst = (a / (a + b + c))[:,0]
    return pairs, np.count_nonzero(is_seg), snp_fst, segSitesPos, is_seg


def fst_per_site(ac, gtvars, poplvl, subsample= True):
    tmpDict = { key: ac[key] for key in ac if key not in ['all'] }
    combs = list(combinations(tmpDict, 2))
    resInf = list()
    resFst = list()
    resFstPerPair = dict()
    is_seg_dict = dict()
    for pair in combs:
        if subsample:
            minlen = min(len(poplvl[pair[0]]), len(poplvl[pair[1]]))
            idx1 = random.sample(poplvl[pair[0]], minlen)
            idx2 = random.sample(poplvl[pair[1]], minlen)
        else:
            idx1 = poplvl[pair[0]]
            idx2 = poplvl[pair[1]]
        pairs, nSegAlleles, fst, segSitesPos, is_seg = calcWCfst_per_site(ac, pair, gtvars, idx1, idx2)
        resInf.append((pair[0], len(idx1), pair[1], len(idx2), nSegAlleles))
        resFst.append(fst)
        pwFstScafBp = pd.DataFrame({'scafbp': segSitesPos, 'wcFst': fst})
        resFstPerPair['_'.join([pair[0], pair[1]])] = pwFstScafBp
        is_seg_dict['_'.join([pair[0], pair[1]])] = is_seg
    resFstDF = pd.DataFrame(resFst).T
    resFstDF.columns = combs
    return pd.DataFrame(resInf, columns= ('pair1', 'n1', 'pair2', 'n2', 'nSegAlleles')), resFstDF, resFstPerPair, is_seg_dict





# --------------------------------------------

if __name__ == "__main__":



    zarrPath = sys.argv[1]


    zarrname = zarrPath.strip('.zarr/')
    # create folders

    varname = zarrname.split('/')[1]

    statsP = os.path.join(zarrname, 'stats/al/')
    figsP = os.path.join(zarrname, 'figs/al/')
    pcafP = os.path.join(figsP, 'pca/')      # pca figs
    pcasP = os.path.join(statsP, 'pca/')     # pca stats
    varDenseP = os.path.join(figsP, 'varDense/')
    hetfP = os.path.join(figsP, 'hets/')
    sfsP = os.path.join(figsP, 'jsfs/')
    fstsP = os.path.join(zarrname, 'stats/al/fst/')
    fstfP = os.path.join(zarrname, 'figs/al/fst/')
    fstsPnests = os.path.join(zarrname, 'stats/al/fst/nests/')
    fstsPpops = os.path.join(zarrname, 'stats/al/fst/pops/')
    gemmasP = os.path.join(zarrname, 'stats/gemma/')



    folderList = [statsP, figsP, pcasP, pcafP, varDenseP, hetfP, sfsP, fstsPnests, fstsPpops, fstfP, gemmasP]
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





    ### allele counts, segregating alleles (for all inds) and heterozygosity

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

    # write pheno file (for gemma)
    phenos = ids[['started_aggression', 'reacted_aggressively', 'reacted_peacefully']]
    phenos.to_csv(os.path.join(gemmasP, '.'.join([ varname, 'pheno' ])), sep= ' ', index= False, header= False)


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
    print("-----------  Calculating pairwise distance matrix  -----------")

    dvar = al.pairwise_distance(gtvars.to_n_alt(), metric= 'cityblock')


    # heatmap with dendrogram

    condensedDvar = scipy.spatial.distance.squareform(dvar)

    n2col = dict(zip(ids['nest'].unique(), sns.color_palette()))
    rowCols = np.array(ids['nest'].map(n2col))

    cDdf = pd.DataFrame(condensedDvar, index= ids['nest'], columns=ids['id'])
    g = sns.clustermap(cDdf, row_colors= rowCols, cmap= 'jet')
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 12)     #ha= 'right', rotation= 40
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 12)

    g.savefig(os.path.join(figsP, 'pwDist.png'), bbox_inches='tight')

    # ---------------------------------------------------------------------

    #############   PCA   #############

    print("-----------  LD pruning, calculating PCAs and plotting figures  -----------")

    plot_ld(nAltVars[:1000], 'Pairwise LD after random downsampling', 'ld_1000.png')



    # random subsampling of loci
    n = 100000  # number of SNPs to choose randomly
    vidxVars = np.random.choice(nAltVars.shape[0], n, replace=False)
    vidxVars.sort()

    gnrVars = nAltVars.take(vidxVars, axis=0)

    plot_ld(gnrVars[:1000], 'Pairwise LD after random downsampling.', 'ld_100k_rand.png')


    segScafs = variants['variants/CHROM'][:][segAll_vars]
    segBP = variants['variants/POS'][:][segAll_vars]
    segVars = pd.DataFrame({'bp': segScafs, 'scaf': segBP})


    # LD pruning


    gnuVars, vars_ldPrd = ld_prune(nAltVars, segVars, size=50, step=20, threshold=.1, n_iter=4)


    plot_ld(gnuVars[:1000], 'Pairwise LD after LD pruning.', filename= 'ld_prune.png')


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

        # nests

    # PCA using SVD - LD-pruned data (59544 loci)
    coords1var, model1var = al.pca(gnuVars, n_components=10, scaler='patterson')


    fig_pca(coords1var, model1var, 'LD-pruned PCA', pops = ids['nest'], pcols= nest_cols, filename= 'pca_LDprune.png')


    ######
    # pca without LD pruning (random subset of 100000 loci)
    coords2var, model2var = al.pca(gnrVars, n_components=10, scaler='patterson')

    # pops
    #fig_pca(coords2var, model2var, 'Conventional PCA', pops = ids['pops'], pcols= pop_cols)

    # nests
    fig_pca(coords2var, model2var, 'Conventional PCA', pops = ids['nest'], pcols= nest_cols, filename= 'pca_100k_rand.png')


    # now for the full set (gtseg_vars)
    coords2allVars, model2allVars = al.pca(nAltVars, n_components=10, scaler='patterson')


    # pops
    #fig_pca(coords2allVars, model2allVars, 'Conventional PCA without LD pruning all variants', pops = ids['pops'], pcols= pop_cols)

    # nests
    fig_pca(coords2allVars, model2allVars, 'Conventional PCA without LD pruning', pops= ids['nest'], pcols= nest_cols, filename= 'pca_all.png')


    # pca with LD pruning, without Patterson's scaling
    coords3vars, model3vars = al.pca(gnuVars, n_components=10, scaler=None)


    # pops
    #fig_pca(coords3vars, model3vars, 'Conventional PCA LD-pruned variants without variance scaling', pops = ids['pops'], pcols= pop_cols)

    # nests
    fig_pca(coords3vars, model3vars, 'Conventional PCA LD-pruned variants without variance scaling.', pops = ids['nest'], pcols= nest_cols, filename= 'pca_LDprune_noPatterson.png')


    # randomized PCA with LD pruning
    coords5vars, model5vars = al.randomized_pca(gnuVars, n_components=10, scaler='patterson')


    # pops
    #fig_pca(coords5vars, model5vars, 'Randomized PCA', pops= ids['pops'], pcols= pop_cols)

    # nests
    fig_pca(coords5vars, model5vars, 'Randomized PCA LD-pruned variants', pops= ids['nest'], pcols= nest_cols, filename= 'pca_LDprune_rand.png')


    plotHeatPCs(coords1var, ids['nest'], PCs=5, filename= 'pca_LDprune_Heat.png')
    plotHeatPCs(coords2allVars, ids['nest'], PCs=5, filename= 'pca_all_Heat.png')

    ## get the Eigen values for PCAs

    # for all (segreg.) vars (n = 1,319,775)
    segScafs = variants['variants/CHROM'][:][segAll_vars]
    segBP = variants['variants/POS'][:][segAll_vars]
    segVars = pd.DataFrame({'bp': segScafs, 'scaf': segBP})
    ## for gemma:
    segVars.to_csv(os.path.join(gemmasP, 'vars_seg.gemma.scafbp'), sep= ' ', index= False, header= False)



    # write first 4 eigen values to file:
    pd.DataFrame([ x[:4] for x in coords1var ]).to_csv(os.path.join(pcasP, "vars_LDprune.eigen"), sep= ' ', index= False, header= False)

    pd.DataFrame([ x[:4] for x in coords2allVars ]).to_csv(os.path.join(pcasP, "vars_all.eigen"), sep= ' ', index= False, header= False)

    # write general shape info to file
    pd.Series([ gtvars.shape[0], gtseg_vars.shape[0], dvar.shape[0], gnuVars.shape[0] ], index= ['all', 'seg. alleles', 'pwDist shape', 'LD-pruned vars']).to_csv(os.path.join(statsP, "nVars.txt"), sep='\t', index= True, header= False)




    ######################

    # plot variant density

    print("--------  Plotting variant density and joint site frequency spectra  --------")

    pos = variants['variants/POS'][:]

    #pos = allel.SortedIndex(variants['variants/POS'])
    # probably have to subset individual scaffolds for this to have more meaning
    plot_windowed_variant_density(pos, window_size=100000, title='Raw variant density', filename= 'varDensity_window100k.png')



    plot_variant_hist('DP', bins= 100, filename= 'var_DP_hist.png')
    plot_variant_hist('AF', filename= 'var_AF_hist.png')


    # plot joint frequency distribution of two variables
    plot_variant_hist_2d('MCOV', 'MQM', downsample=10, filename= 'varDensity_MCOV_MQM_hist2d.png')

    # Ti vs Tv

    #mutations = np.char.add(variants['variants/REF'], variants['variants/ALT'][:, 0])


    #######################

    ## count the number (and proportion) of heterozygous calls

    #gtvars.count_het(axis=0)     # axis (0 = across loci, 1 = across samples)

    propHets = pd.Series(gtvars.count_het(axis= 0)/len(gtvars))
    #missing = gtsub.count_missing(axis=0)[:] / len(gtvars)

    # plot the proportion of heterozygous genotypes



    plotPropHets(propHets, ids, filename= 'propHets.png')



    #def plot_genotype_frequency_pops(pc, title):
    #    fig, ax = plt.subplots(figsize=(12, 4))
    #    sns.despine(ax=ax, offset=10)
    #    left = np.arange(len(pc))
    #    palette = sns.color_palette()
    #    pop2color = {'A': palette[0], 'N': palette[1], 'S': palette[2]}
    #    colors = [pop2color[p] for p in ids['pop'] ]
    #    ax.bar(left, pc, color=colors)
    #    ax.set_xlim(0, len(pc))
    #    ax.set_xlabel('Sample index')
    #    ax.set_ylabel('Percent calls')
    #    ax.set_title(title)
    #    handles = [mpl.patches.Patch(color=palette[0]),
    #               mpl.patches.Patch(color=palette[1]),
    #               mpl.patches.Patch(color=palette[2])]
    #    ax.legend(handles=handles, labels=list(np.unique(ids['pop'])), title='Population',
    #              bbox_to_anchor=(1, 1), loc='upper left')


    ##plot_genotype_frequency(missing, 'Missing')
    #plot_genotype_frequency_pops(propHets, 'Heterozygous')


    #def plot_genotype_frequency_nests(pc, title):
    #    fig, ax = plt.subplots(figsize=(12, 4))
    #    sns.despine(ax=ax, offset=10)
    #    left = np.arange(len(pc))
    #    palette = sns.color_palette()
    #    # change here
    #    pop2color = dict(zip(np.unique(ids['nest']), palette))
    #    colors = [pop2color[p] for p in ids['nest'] ]
    #    ax.bar(left, pc, color=colors)
    #    ax.set_xlim(0, len(pc))
    #    ax.set_xlabel('Sample index')
    #    ax.set_ylabel('Percent calls')
    #    ax.set_title(title)
    #    handles = [mpl.patches.Patch(color=palette[0]),
    #               mpl.patches.Patch(color=palette[1]),
    #               mpl.patches.Patch(color=palette[2]),
    #               mpl.patches.Patch(color=palette[3]),
    #               mpl.patches.Patch(color=palette[4]),
    #               mpl.patches.Patch(color=palette[5]),
    #               mpl.patches.Patch(color=palette[6]),
    #               mpl.patches.Patch(color=palette[7]),
    #               mpl.patches.Patch(color=palette[8])]
    #    ax.legend(handles=handles, labels=list(np.unique(ids['nest'])), title='Nests',
    #              bbox_to_anchor=(1, 1), loc='upper left')
    #
    #
    #plot_genotype_frequency_nests(propHets, 'Heterozygous')







    ############
    # ac_nests_vars
    #ac_seg = ac_subpops['all'].compress(segAll)
    #
    ##########################

    # plot joint SFS for pops


    plot_jsfs(ac_pops_vars, fname= 'jsfs')





    plot_jsfs_nests(ac_nests_vars, fname= 'jsfs_nests')

#       plot folded site frequency spectra for pairs of pops (see https://alimanfoo.github.io/2016/06/10/scikit-allel-tour.html)

########    fig, ax = plt.subplots(figsize=(8, 5))
########    sns.despine(ax=ax, offset=10)
########    sfs1 = allel.stats.sfs_folded_scaled(ac1)
########    allel.stats.plot_sfs_folded_scaled(sfs1, ax=ax, label='BFM', n=ac1.sum(axis=1).max())
########    sfs2 = allel.stats.sfs_folded_scaled(ac2)
########    allel.stats.plot_sfs_folded_scaled(sfs2, ax=ax, label='AOM', n=ac2.sum(axis=1).max())
########    ax.legend()
########    ax.set_title('Scaled folded site frequency spectra')
########    # workaround bug in scikit-allel re axis naming
########    ax.set_xlabel('minor allele frequency');





    # combine plots for each pop in one figure (we'll deal with this, once we've decided what will be used...

    #def plot_jsfs_nests(ac, fname):
    #    tmpDict = { key: ac[key] for key in ac if key not in ['all'] }
    #    for pop in ['A', 'N', 'S']:
    #        popDict = { k: v for k, v in tmpDict.items() if pop in k }
    #        nestCombs = list(combinations(popDict, 2))
    #        fig = plt.figure(figsize=(15, 5))
    #
    #        for i, nestPair in enumerate(nestCombs):
    #            print(i, nestPair)
    #            ax = fig.add_subplot(1,3,i+1)
    #            jsfs = al.stats.sf.joint_sfs(ac[nestPair[0]][:,1], ac[nestPair[1]][:,1])
    #            #fig, ax = plt.subplots(figsize=(6,6))
    #            al.stats.sf.plot_joint_sfs(jsfs, ax=ax)
    #            ax.set_ylabel(' '.join(['Alternate allele count,', nestPair[0] ]))
    #            ax.set_xlabel(' '.join(['Alternate allele count,', nestPair[1] ]))
    #            ax.cla()
    #        #fig.savefig(os.path.join(sfsP, '.'.join([fname, pop, 'png' ])), bbox_inches='tight')
    #        ##p.remove()
    #
    #
    #plot_jsfs_nests(ac_nests_vars, fname= 'jsfs_nests')


    # sample code for the above (from https://stackoverflow.com/questions/19053077/looping-over-data-and-creating-individual-figures#19053157 )

    #fig = plt.figure()
    #ax = fig.addsubplot(111)
    #
    ## Tinker with labels and spines
    #ax.set_xlabel(...)
    #ax.set_ylabel(...)
    #[a.label.set_color('black') for a in (ax.xaxis, ax.yaxis)]
    #...
    #
    ## Plot data and save figures
    #for mol in mols:
    #    for energy, density in data:
    #        ax.cla() # or ax.clear()
    #        p, = ax.plot(density, energy, 'ro')
    #
    #        fig.savefig(mol+".png")
    #        p.remove() # rem



    # -------------------------------------------------------- #

    # Fst

    print("---  Calculating Weir & Cockerham's Fst with std.err from block-jackknife  ---")

    bjFst_nests, nests_bjFstScafBp, nests_bjFst_is_seg = fst_bj_knife_pair(ac_nests_vars, gtvars, poplvl=nests, blen= 10000, subsample=True)
    bjFst_pops, pops_bjFstScafBp, pops_bjFst_is_seg = fst_bj_knife_pair(ac_pops_vars, gtvars, poplvl=pops, blen= 10000, subsample=True)

    bjFst_nests.to_csv(os.path.join(fstsP, 'nests.bjknife.wcFst.evenN.txt' ), header=True, index=False, sep= '\t')
    bjFst_pops.to_csv(os.path.join(fstsP, 'pops.bjknife.wcFst.evenN.txt' ), header=True, index=False, sep= '\t')

for pair, df in nests_bjFstScafBp.items():
    df.to_csv(os.path.join(fstsPnests, '.'.join(['nests', pair, 'wcFst_bjknife.scafbp'])), header= False, index= False, sep= '\t')

for pair, df in pops_bjFstScafBp.items():
    df.to_csv(os.path.join(fstsPpops, '.'.join(['pops', pair, 'wcFst_bjknife.scafbp'])), header= False, index= False, sep= '\t')

    plot_bjFst(bjFst_nests, 'nests')
    plot_bjFst(bjFst_pops, 'pops')



    print("--- Calculating per-site W & C Fst ---")

    pops_fstInf, pops_fstVal, pops_pwFstScafBp, pops_is_seg = fst_per_site(ac_pops_vars, gtvars, poplvl=pops, subsample= True)
    nests_fstInf, nests_fstVal, nests_pwFstScafBp, nests_is_seg  =  fst_per_site(ac_nests_vars, gtvars, poplvl=nests, subsample= True)

    pops_fstInf.to_csv(os.path.join(fstsP, 'pops.wcFst.inf.evenN.txt' ), header=True, index=False, sep= '\t')
    pops_fstVal.to_csv(os.path.join(fstsP, 'pops.wcFst.evenN.txt' ), header=True, index=False, sep= '\t')
    nests_fstInf.to_csv(os.path.join(fstsP, 'nests.wcFst.inf.evenN.txt' ), header=True, index=False, sep= '\t')
    nests_fstVal.to_csv(os.path.join(fstsP, 'nests.wcFst.evenN.txt' ), header=True, index=False, sep= '\t')


    for pair, df in pops_pwFstScafBp.items():
        df.to_csv(os.path.join(fstsPpops, '.'.join(['pops', pair, 'wcFst_persite.scafbp'])), header= True, index= False, sep= '\t')

    for pair, df in nests_pwFstScafBp.items():
        df.to_csv(os.path.join(fstsPnests, '.'.join(['nests', pair, 'wcFst_persite.scafbp'])), header= True, index= False, sep= '\t')





def plot_fst_per_site(df, dfInf, fname):
    fig, ax = plt.subplots(figsize= (7,5))
    df.boxplot()
    plt.xlabel("pairs")
    plt.ylabel("$F_{ST}$")
    ax.set_xticklabels(dfInf['pair1'] + ' - ' + dfInf['pair2'], rotation= 40, ha= 'right', fontsize= 8)
    fig.savefig(os.path.join(fstfP, '.'.join([fname, 'wcFst_perSite_box.pdf'])), bbox_inches='tight')


    plot_fst_per_site(pops_fstVal, pops_fstInf, fname= 'pops')
    plot_fst_per_site(nests_fstVal, nests_fstInf, fname= 'nests')



def getOutlier(dist):
    q1 = dist.quantile(0.25)
    q3 = dist.quantile(0.75)
    iqr = q3-q1 #Interquartile range
    #fence_low  = q1-1.5*iqr
    fence_high = q3+1.5*iqr
    #dist_in = dist.loc[(dist > fence_low) & (dist <= fence_high)]
    out = dist.loc[(dist > fence_high)]
    return out

def getOutlier_idx(dist):
    q1 = dist.quantile(0.25)
    q3 = dist.quantile(0.75)
    iqr = q3-q1 #Interquartile range
    #fence_low  = q1-1.5*iqr
    fence_high = q3+1.5*iqr
    #dist_in = dist.loc[(dist > fence_low) & (dist <= fence_high)]
    out = dist > fence_high
    return out

#pops_fstVal.apply(getOutlier)

pops_fstVal.apply(getOutlier).describe()    # number of outliers, etc.

## now returning the boolean for easier indexing
scafbp[pops_is_seg['A_S']][getOutlier_idx(pops_pwFstScafBp['A_S']['wcFst'])]
# index loci


#idx = al.UniqueIndex(np.arange(0, len(scafbp)))

#scafbp[getScafBp(idx, is_seg)]

print("--- Calculating correlations and plotting heatmaps ---")


fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,2,1)
sns.heatmap(xLDpruned.corr(), cmap= 'coolwarm', annot= True)
ax.set_title('correlation with PCs of LD pruned loci (185,858)')

ax = fig.add_subplot(1,2,2)
sns.heatmap(xSeg.corr(), cmap= 'coolwarm', annot= True)
ax.set_title('correlation with PCs of segregated loci (1,319,775)')
plt.tight_layout()
#['started_aggression', 'reacted_aggressively', 'reacted_peacefully']

fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1,2,1)
sns.heatmap(xLDpruned.cov(), cmap= 'coolwarm', annot= True)
ax.set_title('covariance with PCs of LD pruned loci (185,858)')
plt.tight_layout()

ax = fig.add_subplot(1,2,2)
sns.heatmap(xSeg.cov(), cmap= 'coolwarm', annot= True)
ax.set_title('covariance with PCs of segregating loci (1,319,775)')
plt.tight_layout()



#plt.errorbar( df['Model'], df['Mean'], yerr=df['SD'], fmt='o', color='Black', elinewidth=3,capthick=3,errorevery=1, alpha=1, ms=4, capsize = 5)
#plt.bar(df['Model'], df['Mean'],tick_label = df['Model'])##Bar plot
#plt.xlabel('Model') ## Label on X axis
#plt.ylabel('Average Performance') ##Label on Y axis
