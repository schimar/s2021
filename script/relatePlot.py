#! /usr/bin/python

# This script plots the vcftools relatedness2 values


from sys import argv
import os
import numpy as np
import pandas as pd
import allel as al
import zarr
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
sns.set_style('whitegrid')




if __name__ == "__main__":

    zarrPath = argv[1]

    zarrname = zarrPath.strip('.zarr/')
    # create folders

    vcffP = os.path.join(zarrname, 'figs/vcftools/')

    folderList = [vcffP]
    for folder in folderList:
        if not os.path.exists(folder):
            os.makedirs(folder)



    variants = zarr.open_group(zarrPath, mode='r')
    fpath = argv[2]

    rel2 = pd.read_csv(fpath, sep= '\t')

    ids = pd.read_table('samples109.txt', sep='\t', index_col=False)
    ids['id_nest'] = ids['id'] + '_' + ids['nest']
    ids = ids.sort_values(by='nest')

    nInds = ids.groupby(by= ['nest', 'pop']).count()['sample']
    #np.all(list(variants['samples']) == ids['id_nest'].values)

    samples = list(variants['samples'])
    subsIndex = [samples.index(s) for s in ids['id_nest']]
    ids['subsIndex'] = subsIndex
    ids.sort_values(by=['subsIndex'], inplace= True)

    # Manichaikul relatedness (see Manichaikul2010.pdf table 1)
    ## based on the KING inference. You can interpret the relatedness_phi as the probability to find identical alleles when randomly sampling one allele from each heterozygous individual. So for one individual AB, and the parent AC, there is p=0.25 to choose A from both individuals. That probability is 0.5 when AB is compared to AB.
    ## -->> an estimated kinship coefficient range >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] corresponds to duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree relationships respectively


    m = np.zeros(shape=(109,109))
    m = pd.DataFrame(m, index= ids['id_nest'], columns= ids['id_nest'])
    for i, ind1 in enumerate(m.index):
        for j, ind2 in enumerate(m.columns):
            mani = rel2['RELATEDNESS_PHI'][(rel2['INDV1'] == ind1) & (rel2['INDV2'] == ind2)]
            m.iloc[i,j] = np.float(mani)

    # truncate values < 0 to 0
    m[m < 0] = 0
    rel01 = rel2['RELATEDNESS_PHI'].astype('float')
    rel01[rel01 < 0] = 0

    fig, ax = plt.subplots(figsize=(12,6))
    sns.distplot(rel2['RELATEDNESS_PHI'])
    ax.set_xlabel('Manichaikul relatedness')
    fig.savefig(os.path.join(vcffP, 'rel2.distplot.untruncated.png'), bbox_inches='tight')


    fig, ax = plt.subplots(figsize=(12,6))
    sns.distplot(rel01) #rel2['RELATEDNESS_PHI'])
    ax.set_xlabel('Manichaikul relatedness')
    #plt.show()
    fig.savefig(os.path.join(vcffP, 'rel2.distplot.png'), bbox_inches='tight')

    fig, ax = plt.subplots(figsize=(16,16))
    sns.heatmap(m, cmap='viridis', annot=False)
    ax.set_title("Relatedness matrix (Manichaikul - KING)")
    fig.tight_layout()
    fig.savefig(os.path.join(vcffP, 'rel2.heat.png'), bbox_inches='tight')



