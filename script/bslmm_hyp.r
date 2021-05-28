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
mainDir <- dirname(infile)
subDir <- unlist(strsplit(basename(infile), '.', fixed= T))[1]
#substr(basename(infile), 1, nchar(basename(infile))-4)
outDir <- file.path(mainDir, subDir)

if (!dir.exists(outDir)){
	dir.create(outDir)
} else {
    print(paste0(outDir, " already exists"))
}

outfile <- paste0(file.path(outDir, subDir), '.hyp.pdf')

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

# get table of hyperparameters

# Get mean, median, and 95% ETPI of hyperparameters
# ========================================================================
# h-> approximation to proportion of phenotypic variance
# explained by variants (PVE)
h <- c("h",mean(hyp$h),quantile(hyp$h, probs=c(0.5,0.025,0.975)))
# pve -> PVE
pve <- c("PVE", mean(hyp$pve),quantile(hyp$pve,
probs=c(0.5,0.025,0.975)))
# rho-> approximation to proportion of genetic variance explained by variants
# with major effect (PGE)
# rho=0 -> pure LMM, highly polygenic basis
# rho=1 -> pure BVSR, few major effect loci
rho <- c("rho",mean(hyp$rho),quantile(hyp$rho, probs=c(0.5,0.025,0.975)))
# pge -> PGE
pge <- c("PGE",mean(hyp$pge),quantile(hyp$pge, probs=c(0.5,0.025,0.975)))
# pi -> proportion of variants with non-zero effects
pi <- c("pi",mean(hyp$pi),quantile(hyp$pi, probs=c(0.5,0.025,0.975)))
# n.gamma -> number of variants with major effect
n.gamma <- c("n.gamma",mean(hyp$n_gamma),quantile(hyp$n_gamma,
probs=c(0.5,0.025,0.975)))

# get table of hyperparameters and save it to a file
# ==============================================================================
hyp.table <- as.data.frame(rbind(h,pve,rho,pge,pi,n.gamma),row.names=F)
colnames(hyp.table) <- c("hyperparam", "mean","median","2.5%", "97.5%")
# show table
hyp.table

# write table to file
write.table(hyp.table, file= file.path(outDir, "hyperparameters.dsv"), sep="\t", quote=F)



