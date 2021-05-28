# library to speed up loading of big tables
library(data.table)

argv <- commandArgs()
infile <- argv[6]

#infile <- "vars/taSubInDel/stats/gemma/output/tst_enq.param.txt"



## Provide the dir name(i.e sub dir) that you want to create under main dir:
mainDir <- dirname(infile)
subDir <- unlist(strsplit(basename(infile), '.', fixed= T))[1]
#substr(basename(infile), 1, nchar(basename(infile))-4)
outDir <- file.path(mainDir, subDir)

if (!dir.exists(outDir)){
	dir.create(outDir)
} else {
    print(paste0(outDir, " already exists"))
}

ofs <- file.path(outDir, subDir)

#basename(infile)
#dirname(infile) 

# Load parameters
# ========================================================================
params <- as.data.frame(fread(infile, header=T, sep="\t"))

chr <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(params$rs, "_"), '[[', 2)), ":"), "[[", 1)))
ps <- unlist(lapply(strsplit(params$rs, ":"), '[[', 3))
params["chr"] <- chr
params$ps <- as.numeric(ps)

# ========================================================================
# Get variants with sparse effect size on phenotypes
# ==============================================================================

# add sparse effect size (= beta * gamma) to data frame
params["eff"] <- abs(params$beta*params$gamma)

# get variants with effect size > 0
pareffs <- params[params$eff>0, ]

# show number of variants with measurable effect
## nrow(pareffs)

# sort by decreasing effect size
pareffs.sort <- pareffs[order(-pareffs$eff), ]

# show top 10 variants with highest effect
##head(pareffs.sort, 10)

# variants with the highest sparse effects
# ------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1 <- pareffs.sort[pareffs.sort$eff > quantile(pareffs.sort$eff, 0.99),]
# top 0.1% variants (above 99.9% quantile)
top01 <- pareffs.sort[pareffs.sort$eff > quantile(pareffs.sort$eff, 0.999),]
# top 0.01% variants (above 99.99% quantile)
top001 <- pareffs.sort[pareffs.sort$eff > quantile(pareffs.sort$eff, 0.9999),]
# ------------------------------------------------------------------------

# write tables
write.table(top1, file= paste0(outDir, "/top1eff.dsv"), quote=F, row.names=F, sep="\t")
write.table(top01, file= paste0(outDir, "/top0.1eff.dsv"), quote=F, row.names=F, sep="\t")
write.table(top001, file= paste0(outDir, "/top0.01eff.dsv"), quote=F, row.names=F, sep="\t")

# ==========================================================================

# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC
# (the number of times it appears in the MCMC with a non-zero sparse effect)
# sort variants by descending PIP
params.pipsort <- params[order(-params$gamma), ]
# Show top 10 variants with highest PIP
head(params.pipsort, 10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01 <- params.pipsort[params.pipsort$gamma >= 0.01,]
# variants with effect in 10% MCMC samples or more
pip10 <- params.pipsort[params.pipsort$gamma >= 0.10,]
# variants with effect in 25% MCMC samples or more
pip25 <- params.pipsort[params.pipsort$gamma >= 0.25,]
# variants with effect in 50% MCMC samples or more
pip50 <- params.pipsort[params.pipsort$gamma >= 0.50,]
# write tables
write.table(pip01, file= paste0(outDir, "/pip01.dsv"), quote=F, row.names=F, sep="\t")
write.table(pip10, file= paste0(outDir, "/pip10.dsv"), quote=F, row.names=F, sep="\t")
write.table(pip25, file= paste0(outDir, "/pip25.dsv"), quote=F, row.names=F, sep="\t")
write.table(pip50, file= paste0(outDir, "/pip50.dsv"), quote=F, row.names=F, sep="\t")
# -----------------------------------------------------------------------

# plot variants PIPs across linkage groups/chromosomes
# ==============================================================================
# Prepare data
# ------------------------------------------------------------------------------
# add linkage group column (chr)
#chr <- gsub("lg|_.+","",params$rs)

# sort by linkage group and position
params.sort <- params[order(as.numeric(params$chr), params$ps),]
# get list of linkage groups/chromosomes
chrs <- sort(as.numeric(unique(chr)))
# ------------------------------------------------------------------------------
# Plot to a png file because the number of dots is very high
# drawing this kind of plot over the network is very slow
# also opening vectorial files with many objects is slow
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
png(file= paste0(outDir, "/pip_plot.png"), width=11.7,height=8.3,units="in",res=200)
# set up empty plot
plot(-1,-1,xlim=c(0,nrow(params.sort)),ylim=c(0,1),ylab="PIP",xlab="linkage group",
xaxt="n")

# plot grey bands for chromosome/linkage groups
# ------------------------------------------------------------------------------
chrs <- sort(as.numeric(unique(chr)))
start <- 1
lab.pos <- vector()
for (ch in chrs){
size <- nrow(params.sort[params.sort$chr==ch,])
cat ("CH: ", ch, "\n")
colour <- "light grey"
if (ch%%2 > 0){
polygon(c(start,start,start+size,start+size,start), c(0,1,1,0,0), col=colour,
border=colour)
}
cat("CHR: ", ch, " variants: ", size, "(total: ", (start+size), ")\n")
txtpos <- start+size/2
lab.pos <- c(lab.pos, txtpos)
start <- start + size
}# Add variants outside linkage groups
chrs <- c(chrs,"NA")
size <- nrow(params.sort[params.sort$chr == "NA", ])
lab.pos <- c(lab.pos, start+size/2)
# ------------------------------------------------------------------------------
# Add x axis labels
axis(side=1, at= lab.pos, labels= chrs, tick= F)

# plot PIP for all variants
# ------------------------------------------------------------------------------
# rank of variants across linkage groups
x <- seq(1,length(params.sort$gamma), 1)
# PIP
y <- params.sort$gamma
# sparse effect size, used for dot size
z <- params.sort$eff
# log-transform to enhance visibility
z[z == 0] <- 0.00000000001
z <- 1/abs(log(z))
# plot
symbols(x, y, circles= z, bg= "black", inches=1/5, fg=NULL, add=T)
# ------------------------------------------------------------------------------

# highligh high PIP variants (PIP>=0.25)
# ------------------------------------------------------------------------------
# plot threshold line
abline(h= 0.25, lty= 3, col= "dark grey")
# rank of high PIP variants across linkage groups
x <- match(params.sort$gamma[params.sort$gamma >= 0.25], params.sort$gamma)
# PIP
y <- params.sort$gamma[params.sort$gamma >= 0.25]
# sparse effect size, used for dot size
z <- params.sort$eff[params.sort$gamma >= 0.25]
z <- 1/abs(log(z))

if(length(z) > 0) {
	symbols(x, y, circles= z, bg="red", inches=1/5, fg=NULL, add=T)
	# add labels for high PIP variants
	text(x, y, labels= params.sort$rs[params.sort$gamma >= 0.25], adj= c(0,0), cex=0.5)
} else {
	print("no variants with PIP >= 0.25")
}
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# close device
dev.off()
# ==============================================================================
#### This is to be done outside the current R session. Launch another interactive session
#### in Iceberg and execute:
# $ display -resize 1920x1080 output/pip_plot.png


