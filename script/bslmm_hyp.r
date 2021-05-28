#! /usr/bin/R


#setwd("/home/schimar/ui/ta/ta_dna_snakemake_pbs/vars/taSubInDel/stats/gemma/output/")
## Load hyperparameters
## ========================================================================

argv <- commandArgs()
infile <- argv[6]

hyp <- read.table(infile, header=T)
## ========================================================================

# from http://romainvilloutreix.alwaysdata.net/romainvilloutreix/wp-content/uploads/2017/01/gwas_gemma-2017-01-17.pdf

#path <- dirname(infile)

outfile <- paste0(substr(infile, 1, nchar(infile)-3), 'pdf')
print(outfile)

# plot traces and distributions of hyperparameters
# ==============================================================================
pdf(file=outfile, width=8.3,height=11.7)
layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))

# PVE
# ------------------------------------------------------------------------------
plot(hyp$pve, type="l", ylab="PVE", main="PVE - trace")
hist(hyp$pve, main="PVE - posterior distribution", xlab="PVE")
plot(density(hyp$pve), main="PVE - posterior distribution", xlab="PVE")
# ------------------------------------------------------------------------------

# PGE
# ------------------------------------------------------------------------------
plot(hyp$pge, type="l", ylab="PGE", main="PGE - trace")
hist(hyp$pge, main="PGE - posterior distribution", xlab="PGE")
plot(density(hyp$pge), main="PGE - posterior distribution", xlab="PGE")
# ------------------------------------------------------------------------------

# pi
# ------------------------------------------------------------------------------
plot(hyp$pi, type="l", ylab="pi", main="pi")
hist(hyp$pi, main="pi", xlab="pi")
plot(density(hyp$pi), main="pi", xlab="pi")
# ------------------------------------------------------------------------------

# No gamma
# ------------------------------------------------------------------------------
plot(hyp$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
hist(hyp$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
plot(density(hyp$pi), main="n_gamma - posterior distribution", xlab="n_gamma")
# ------------------------------------------------------------------------------
dev.off()
# ==============================================================================






