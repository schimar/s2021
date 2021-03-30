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
    fig.savefig(os.path.join(fstfP, '.'.join([fname, 'bjFst_bar.png'])), bbox_inches='tight')




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


def plot_fst_per_site(df, dfInf, fname):
    fig, ax = plt.subplots(figsize= (7,5))
    df.boxplot()
    plt.xlabel("pairs")
    plt.ylabel("$F_{ST}$")
    ax.set_xticklabels(dfInf['pair1'] + ' - ' + dfInf['pair2'], rotation= 40, ha= 'right', fontsize= 8)
    fig.savefig(os.path.join(fstfP, '.'.join([fname, 'wcFst_perSite_box.png'])), bbox_inches='tight')


def plot_fsfs(gtvars, pops, fname, cols, subsample=True):
    popDict = { key: pops[key] for key in pops if key not in ['all'] }
    popDictlen = {key: len(value) for key, value in popDict.items()}
    if subsample:
        minIDlen = min(popDictlen.values())
        tmpPopDict = {key: random.sample(value, minIDlen) for key, value in popDict.items()}
        ac = gtvars.count_alleles_subpops(tmpPopDict, max_allele=1)
    else:
        tmpPopDict = popDict
        ac = gtvars.count_alleles_subpops(tmpPopDict, max_allele=1)
    fig, ax = plt.subplots(figsize=(8,8))
    acs = list()
    for key, alco in ac.items():
        sns.despine(ax=ax, offset=5)
        actmp = ac[key]
        sfs = al.stats.sf.sfs_folded_scaled(actmp)
        al.stats.sf.plot_sfs_folded(sfs, ax=ax, label=key, plot_kwargs={'color': cols[key]}, n=None)
        ax.legend()
        if subsample:
            ax.set_title('Folded site frequency spectra (' + str(minIDlen) + ' ind.)')
        else:
            ax.set_title('Folded site frequency spectra (uneven N inds)')
        fig.savefig(os.path.join(sfsP, '.'.join([fname, 'png' ])), bbox_inches='tight')


def garud_h_per_pop(gtvars, pops, downsamp= False):
    resDict = {}
    nds = min([ len(v) for k, v in pops.items() ])
    for pop in pops:
        if downsamp:
            ds = random.sample(pops[pop], nds)
            ht = gtvars.subset(sel0=ds, sel1=pops[pop]).to_haplotypes()
        else:
            ht = gtvars.subset(sel1=pops[pop]).to_haplotypes()
        h = al.garud_h(ht)
        resDict[pop] = h
    return pd.DataFrame(resDict, index= ('h1', 'h12', 'h123', 'h2_h1'))


def biplot2d(score,coeff,labels=None, pclabs=[1,2]):
    xs = score[:,0]
    ys = score[:,1]
    n = coeff.shape[0]
    scalex = 1.0/(xs.max() - xs.min())
    scaley = 1.0/(ys.max() - ys.min())
    fig = plt.figure(figsize=(12,8))
    plt.scatter(xs * scalex,ys * scaley,s=5)
    for i in range(n):
        plt.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
        if labels is None:
            plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, "Var"+str(i+1), color = 'green', ha = 'center', va = 'center')
        else:
            plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'g', ha = 'center', va = 'center')

    plt.xlabel("PC{}".format(pclabs[0]))
    plt.ylabel("PC{}".format(pclabs[1]))
    plt.grid()
    plt.tight_layout()
    fig.savefig(os.path.join(varpcafP,'.'.join(['biplot_pc', str(pclabs[0]), str(pclabs[1]), 'png'])), bbox_inches='tight')


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
    varpcafP = os.path.join(pcafP, 'varPca/')
    varpcasP = os.path.join(pcasP, 'varPca/')
    varDenseP = os.path.join(figsP, 'varDense/')
    hetfP = os.path.join(figsP, 'hets/')
    sfsP = os.path.join(figsP, 'sfs/')
    fstsP = os.path.join(zarrname, 'stats/al/fst/')
    fstfP = os.path.join(zarrname, 'figs/al/fst/')
    fstsPnests = os.path.join(zarrname, 'stats/al/fst/nests/')
    fstsPpops = os.path.join(zarrname, 'stats/al/fst/pops/')
    gemmasP = os.path.join(zarrname, 'stats/gemma/')
    selsP = os.path.join(zarrname, 'stats/al/sel/')
    selfP = os.path.join(zarrname, 'figs/al/sel/')


    folderList = [statsP, figsP, pcasP, pcafP, varpcafP, varpcasP, varDenseP, hetfP, sfsP, fstsPnests, fstsPpops, fstfP, selsP, selfP, gemmasP]
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



    # ------------------------------------



    #############   selective sweeps    #############

    print("--- Calculating and plotting Garud's H statistics ---")

            #  Garud's H (signatures of (soft vs hard) selective sweeps) (see Garud (2015) Recent Selective Sweeps in North American Drosophila melanogaster Show Signatures of Soft Sweeps)
    # (written to vars/ta{vartype}/{stats,figs}/al/sel/)


    # NOTE: for now, I've only written Garud's H stats and figs for the full sample size per pop (downsamp in garud_h_per_pop func)


    garud_h_per_pop(gtseg_vars, pops).to_csv(os.path.join(selsP, 'pops.garud_h.txt'), header= True, index= True, sep= '\t')

    garud_h_per_pop(gtseg_vars, nests).to_csv(os.path.join(selsP, 'nests.garud_h.txt'), header= True, index= True, sep= '\t')



    fig = plt.figure(figsize=(6, 5))
    ax = sns.heatmap(garud_h_per_pop(gtseg_vars, pops), cmap='coolwarm', annot= True)
    ax.set_title("Garud's H statistics - pops")
    fig.tight_layout()
    fig.savefig(os.path.join(selfP, 'pops.garud_h.png'), bbox_inches='tight')
    #fig.suptitle(title, y=1.02)

    fig = plt.figure(figsize=(10, 5))
    ax = sns.heatmap(garud_h_per_pop(gtseg_vars, nests), cmap='coolwarm', annot= True)
    ax.set_title("Garud's H statistics - nests")
    fig.tight_layout()
    fig.savefig(os.path.join(selfP, 'nests.garud_h.png'), bbox_inches='tight')







    #############   PCA gtvars  #############

    # (in vars/ta{vartype}/{stats,figs}/al/pca/)

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
    segVars = pd.DataFrame({'scaf': segScafs, 'bp': segBP})


    # LD pruning


    gnuVars, vars_ldPrd = ld_prune(nAltVars, segVars, size=50, step=20, threshold=.1, n_iter=4)

    vars_ldPrd.to_csv(os.path.join(pcasP, 'ld_prunedVars.scafbp.txt'), header= False, index= False, sep= '\t')

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
    # which one is the outlier in the LD-pruned PCA?
    ##np.where(coords1var[:,0] > 200)
    ##ids.iloc[59]    # 101a_S1

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

    # for all (segreg.) vars
    segScafs = variants['variants/CHROM'][:][segAll_vars]
    segBP = variants['variants/POS'][:][segAll_vars]
    segVars = pd.DataFrame({'bp': segScafs, 'scaf': segBP})
    ## for gemma:
    segVars.to_csv(os.path.join(gemmasP, 'vars_seg.gemma.scafbp'), sep= ' ', index= False, header= False)



    # write first 4 eigen values to file:
    pd.DataFrame([ x[:4] for x in coords1var ]).to_csv(os.path.join(pcasP, "vars_LDprune.loadings.txt"), sep= ' ', index= False, header= False)

    pd.DataFrame([ x[:4] for x in coords2allVars ]).to_csv(os.path.join(pcasP, "vars_all.loadings.txt"), sep= ' ', index= False, header= False)

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



    #
    ##########################

    # plot joint SFS for pops


    plot_jsfs(ac_pops_vars, fname= 'pops.jsfs')


    plot_jsfs_nests(ac_nests_vars, fname= 'nests.jsfs')



    # plot folded site frequency spectra for pairs of pops (see https://alimanfoo.github.io/2016/06/10/scikit-allel-tour.html)

    plot_fsfs(gtvars, pops, fname='pops.fsfs', cols= pop_cols, subsample=True)
    plot_fsfs(gtvars, pops, fname='pops.fsfs.unevenN', cols= pop_cols, subsample=False)

    plot_fsfs(gtvars, nests, fname='nests.fsfs', cols= nest_cols, subsample=True)
    plot_fsfs(gtvars, nests, fname='nests.fsfs.unevenN', cols= nest_cols, subsample=False)

    # when excess of rare variants, then suggesting a population expansion
    # when closer to neutral expectation, suggesting a more stable population size.





    #############   Fst   #############


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


    plot_fst_per_site(pops_fstVal, pops_fstInf, fname= 'pops')
    plot_fst_per_site(nests_fstVal, nests_fstInf, fname= 'nests')






    #############   PCA over variables (ids, PCs of gtvars and more)   #############




    df1 = ids[['elevation', 'agg', 'MMAI_worker', 'AI_worker', 'lat', 'lon']]
    df1.columns = ['elev', 'agg', 'MMAI', 'AI', 'lat', 'lon']

    df1.loc[:,'het'] = propHets
    df1.loc[:, 'nest'] = pd.factorize(ids['nest'])[0]
    df1.loc[:, 'pop'] = pd.factorize(ids['pop'])[0]
    df1.loc[:, 'gyn'] = 0
    df1['gyn'][df1['pop'] == 2] = 1


    df = pd.concat([df1, pd.DataFrame(coords1var[:,:4]), pd.DataFrame(coords2allVars[:,:4])], axis=1, join='inner')
    df.columns = ['elev', 'agg', 'MMAI', 'AI', 'lat', 'lon', 'het', 'nest', 'pop', 'gyn', 'pc1ldp', 'pc2ldp', 'pc3ldp', 'pc4ldp', 'pc1av', 'pc2av', 'pc3av', 'pc4av']



    df_st =  StandardScaler().fit_transform(df)

    pca_out = PCA(n_components=10).fit(df_st)
    # Proportion of Variance (from PC1 to PC6)
    pca_out.explained_variance_ratio_

    # Cumulative proportion of variance (from PC1 to PC6)
    np.cumsum(pca_out.explained_variance_ratio_)

    # component loadings (correlation coefficient between original variables and the component)
    # the squared loadings within the PCs always sums to 1
    loadings = pca_out.components_
    num_pc = pca_out.n_features_
    pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
    loadings_df = pd.DataFrame.from_dict(dict(zip(pc_list, loadings)))
    loadings_df['variable'] = df.columns.values
    loadings_df = loadings_df.set_index('variable')
    loadings_df.to_csv(os.path.join(varpcasP, 'varPca.loadings.txt'), header= True, index= False, sep= '\t')
    ## NOTE: write to file

    # positive and negative values in component loadings reflects the positive and negative
    # correlation of the variables with the PCs (they have a "positive projection" on first PC)


    # get correlation matrix plot for loadings
    fig = plt.figure(figsize=(14,12))
    sns.heatmap(loadings_df, annot=True, cmap='Spectral')
    plt.tight_layout()
    fig.savefig(os.path.join(varpcafP,'corrMat_loadPCA.png'), bbox_inches='tight')
    fig.clear()


    # screeplot with cumsum variance
    fig = plt.figure(figsize=(12,8))
    plt.bar(range(1,len(pca_out.explained_variance_ratio_ )+1),pca_out.explained_variance_ratio_ )
    plt.ylabel('Explained variance')
    plt.xlabel('Components')
    plt.plot(range(1,len(pca_out.explained_variance_ratio_ )+1),
             np.cumsum(pca_out.explained_variance_ratio_),
             c='red',
             label="Cumulative Explained Variance")
    plt.legend(loc='upper left')
    fig.savefig(os.path.join(varpcafP,'scree_cumsumVar.png'), bbox_inches='tight')
    fig.clear()

    # get PC scores
    pca_scores = PCA().fit_transform(df_st)


    # 2d biplot       https://ostwalprasad.github.io/machine-learning/PCA-using-python.html
    biplot2d(pca_scores[:,0:2],np.transpose(pca_out.components_[0:2, :]),list(df.columns), pclabs=[1,2])

    biplot2d(pca_scores[:,2:4],np.transpose(pca_out.components_[2:4, :]),list(df.columns), pclabs=[2,3])

    # NOTE: see https://www.reneshbedre.com/blog/principal-component-analysis.html#pca-loadings-plots
    # get eigenvalues (from PC1 to PC6)
    #pca_out.explained_variance_

    # get scree plot (for scree or elbow test)
    #from bioinfokit.visuz import cluster
    #cluster.screeplot(obj=[pc_list, pca_out.explained_variance_ratio_])

    # Scree plot will be saved in the same directory with name screeplot.png
    # get PCA loadings plots (2D and 3D)
    # 2D
    #cluster.pcaplot(x=loadings[0], y=loadings[1], labels=df.columns.values,
    #    var1=round(pca_out.explained_variance_ratio_[0]*100, 2),
    #    var2=round(pca_out.explained_variance_ratio_[1]*100, 2), show= False)
    #
    ## 3D
    #cluster.pcaplot(x=loadings[0], y=loadings[1], z=loadings[2],  labels=df.columns.values,
    #    var1=round(pca_out.explained_variance_ratio_[0]*100, 2), var2=round(pca_out.explained_variance_ratio_[1]*100, 2),
    #    var3=round(pca_out.explained_variance_ratio_[2]*100, 2))

    # get 2D biplot
    #cluster.biplot(cscore=pca_scores, loadings=loadings, labels=df.columns.values, var1=round(pca_out.explained_variance_ratio_[0]*100, 2),
    #    var2=round(pca_out.explained_variance_ratio_[1]*100, 2), valphadot=0.4, dotsize= 3, datapoints=True, colordot= 'lightgrey')#, colorlist=ids['nest'])
    #
    ## get 3D biplot
    #cluster.biplot(cscore=pca_scores, loadings=loadings, labels=df.columns.values,
    #    var1=round(pca_out.explained_variance_ratio_[0]*100, 2), var2=round(pca_out.explained_variance_ratio_[1]*100, 2),
    #    var3=round(pca_out.explained_variance_ratio_[2]*100, 2), valphadot=0.4, dotsize= 3, datapoints=True, colordot='lightgrey')#, colorlist=ids['nest'])





    #############   correlations for all variables   #############


    print("--- Calculating correlations and plotting heatmaps ---")

    fig.clear()
    fig = plt.figure(figsize=(14,12))
    sns.heatmap(df.corr(), cmap='coolwarm', annot=True)
    plt.tight_layout()
    plt.savefig(os.path.join(varpcafP,'df_corr_heat.png'), bbox_inches='tight')





# -----------------------------------------------------

#idx = al.UniqueIndex(np.arange(0, len(scafbp)))

#scafbp[getScafBp(idx, is_seg)]


#pops_fstVal.apply(getOutlier)

#pops_fstVal.apply(getOutlier).describe()    # number of outliers, etc.

## now returning the boolean for easier indexing
#scafbp[pops_is_seg['A_S']][getOutlier_idx(pops_pwFstScafBp['A_S']['wcFst'])]



    #############   Haplotype diversity   #############


# --------------------------------------------------------------
# NOTE: haplo div is 1, what does it mean?
#def haplo_div(gtvars, pops):
#    resDict = {}
#    for pop in pops:
#        ht = gtvars.subset(sel1=pops[pop]).to_haplotypes()
#        resDict[pop] = al.haplotype_diversity(ht)
#    return resDict
#
#haplo_div(gtvars, pops)


    #############   dxy (sequence_divergence)   #############
#    def get_scafs_iter:
#
#
#        ac =
#    for i, s in enumerate(np.unique(scaf)):
#        isScaf = scaf == s
#        tmp = gtvars.compress(isScaf, axis = 0)
#        actmp = tmp.count_alleles_subpops(
#        print(tmp.shape)
#
#
#    ac_pops_vars = gtvars.count_alleles_subpops(pops, max_allele=1)
#        segAll_vars = ac_pops_vars['all'].is_segregating()[:]
#        gtseg_vars = gtvars.compress(segAll_vars, axis=0)
#        nAltVars = gtseg_vars.to_n_alt()
#
#        # take subset of vars on scaf s
#        # get allele counts for pop pairs
#        # calc sequence_divergence (windowed?) for the two allele counts
#
#
#    gtvars.subset(sel0=sel1=pops['S'])






##########  misc  #########

#
#fig = plt.figure(figsize=(14,6))
#ax = fig.add_subplot(1,2,1)
#sns.heatmap(xLDpruned.corr(), cmap= 'coolwarm', annot= True)
#ax.set_title('correlation with PCs of LD pruned loci (185,858)')
#
#ax = fig.add_subplot(1,2,2)
#sns.heatmap(xSeg.corr(), cmap= 'coolwarm', annot= True)
#ax.set_title('correlation with PCs of segregated loci (1,319,775)')
#plt.tight_layout()
##['started_aggression', 'reacted_aggressively', 'reacted_peacefully']
#
#fig = plt.figure(figsize=(14,6))
#ax = fig.add_subplot(1,2,1)
#sns.heatmap(xLDpruned.cov(), cmap= 'coolwarm', annot= True)
#ax.set_title('covariance with PCs of LD pruned loci (185,858)')
#plt.tight_layout()
#
#ax = fig.add_subplot(1,2,2)
#sns.heatmap(xSeg.cov(), cmap= 'coolwarm', annot= True)
#ax.set_title('covariance with PCs of segregating loci (1,319,775)')
#plt.tight_layout()


