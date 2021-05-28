#! /usr/bin/R

argv <- commandArgs()

print(argv)
# from http://romainvilloutreix.alwaysdata.net/romainvilloutreix/wp-content/uploads/2017/01/gwas_gemma-2017-01-17.pdf

#setwd("/home/schimar/ui/ta/ta_dna_snakemake_pbs/vars/taSubInDel/stats/gemma/output/")
## Load hyperparameters
## ========================================================================
#hyp <- read.table("bslmm.hyp.txt",header=T)
## ========================================================================
#
## Get mean, median, and 95% ETPI of hyperparameters
## ========================================================================
## h-> approximation to proportion of phenotypic variance
## explained by variants (PVE)
#h <- c("h",mean(hyp$h),quantile(hyp$h, probs=c(0.5,0.025,0.975)))
## pve -> PVE
#pve <- c("PVE", mean(hyp$pve),quantile(hyp$pve,
#probs=c(0.5,0.025,0.975)))
## rho-> approximation to proportion of genetic variance explained by variants
## with major effect (PGE)
## rho=0 -> pure LMM, highly polygenic basis
## rho=1 -> pure BVSR, few major effect loci
#rho <- c("rho",mean(hyp$rho),quantile(hyp$rho, probs=c(0.5,0.025,0.975)))
## pge -> PGE
#pge <- c("PGE",mean(hyp$pge),quantile(hyp$pge, probs=c(0.5,0.025,0.975)))
## pi -> proportion of variants with non-zero effects
#pi <- c("pi",mean(hyp$pi),quantile(hyp$pi, probs=c(0.5,0.025,0.975)))
## n.gamma -> number of variants with major effect
#n.gamma <- c("n.gamma",mean(hyp$n_gamma),quantile(hyp$n_gamma,
#probs=c(0.5,0.025,0.975)))
## ==============================================================================
## Analysing BSLMM 
#
## get table of hyperparameters and save it to a file
## ==============================================================================
#hyp.table <- as.data.frame(rbind(h,pve,rho,pge,pi,n.gamma),row.names=F)
#colnames(hyp.table) <- c("hyperparam", "mean","median","2.5%", "97.5%")
## show table
#hyp.table
## write table to file
#write.table(hyp.table, file="hyperparameters.dsv", sep="\t", quote=F)
## ==============================================================================
## Table should look like this:
#hyp.table
#
## plot traces and distributions of hyperparameters
## ==============================================================================
## set up layout
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
## h
## ------------------------------------------------------------------------------
#plot(hyp$h, type="l", ylab="h", main="h")
#hist(hyp$h, main="", xlab="h")
#plot(density(hyp$h), main="", xlab="h")
## ------------------------------------------------------------------------------
#
#
## plot traces and distributions of hyperparameters
## ==============================================================================
## set up layout
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
## h
## ------------------------------------------------------------------------------
#plot(hyp$h, type="l", ylab="h", main="h")
#hist(hyp$h, main="", xlab="h")
#plot(density(hyp$h), main="", xlab="h")
## ------------------------------------------------------------------------------
## PVE
## ------------------------------------------------------------------------------
#plot(hyp$pve, type="l", ylab="PVE", main="PVE")
#hist(hyp$pve, main="", xlab="PVE")
#plot(density(hyp$pve), main="", xlab="PVE")
## ------------------------------------------------------------------------------
## rho
## ------------------------------------------------------------------------------
#plot(hyp$rho, type="l", ylab="rho", main="rho")
#hist(hyp$rho, main="", xlab="rho")
#plot(density(hyp$rho), main="", xlab="rho")
## ------------------------------------------------------------------------------
#
#
#
#
#
#
#
## plot traces and distributions of hyperparameters
## ==============================================================================
#pdf(file="hyperparameters.pdf", width=8.3,height=11.7)
#layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))
#
## PVE
## ------------------------------------------------------------------------------
#plot(hyp$pve, type="l", ylab="PVE", main="PVE - trace")
#hist(hyp$pve, main="PVE - posterior distribution", xlab="PVE")
#plot(density(hyp$pve), main="PVE - posterior distribution", xlab="PVE")
## ------------------------------------------------------------------------------
#
## PGE
## ------------------------------------------------------------------------------
#plot(hyp$pge, type="l", ylab="PGE", main="PGE - trace")
#hist(hyp$pge, main="PGE - posterior distribution", xlab="PGE")
#plot(density(hyp$pge), main="PGE - posterior distribution", xlab="PGE")
## ------------------------------------------------------------------------------
#
## pi
## ------------------------------------------------------------------------------
#plot(hyp$pi, type="l", ylab="pi", main="pi")
#hist(hyp$pi, main="pi", xlab="pi")
#plot(density(hyp$pi), main="pi", xlab="pi")
## ------------------------------------------------------------------------------
#
## No gamma
## ------------------------------------------------------------------------------
#plot(hyp$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
#hist(hyp$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
#plot(density(hyp$pi), main="n_gamma - posterior distribution", xlab="n_gamma")
## ------------------------------------------------------------------------------
#dev.off()
## ==============================================================================
#
#
#
#
#
#
