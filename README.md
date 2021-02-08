# ta_dna_snakemake_pbs

======================================================

Snakemake workflow for *Tetramorium alpestre* DNA analysis on the mach2 HPC cluster with the PBS-torque batch job submission system. This repository was initially based on the [Snakemake Cluster Tutorial](https://github.com/SchlossLab/snakemake_cluster_tutorial.git) and the [Software Carpentry lesson repository](https://hpc-carpentry.github.io/hpc-python/17-cluster/). If you want to run your own custom workflow, check out the [insert link here]() generic workflow. 


======================================================

## conda and other [dependencies](https://github.com/schimar/ta_dna_snakemake_pbs/blob/main/envs/s21.yaml)   

create environment from yaml file (in envs/):
```
# run these two once, to create the environment:
conda init bash
conda env create -f envs/s21.yaml

# with this, you can activate the environment with all [dependencies](https://github.com/schimar/ta_dna_snakemake_pbs/blob/main/envs/s21.yaml):
conda activate ta

# (also, when ssh'ing onto mach2, you can activate the env and then do a dry-run of your workflow) 
```

## how to submit the main snakemake job:
```
qsub code/clusterSnakemake.pbs
```




