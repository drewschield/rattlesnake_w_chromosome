############################################################################
# GC and CpG content of Z(X), W(Y), and autosomes for amniotes
############################################################################

### Goal: compare GC and CpG content between sex chromsomes and autosomes for
### the prairie rattlesnake and a sample of bird and mammal genomes.

### Working directory and dependencies--------------------------------------

setwd('./gc')

install.packages("ppcor")
library(ppcor)
library(car)

### Read in data------------------------------------------------------------

## GC content 
cvv.gc.w <- read.table('Cviridis_CV0650_candidate_W.GC.txt',header=T)

cvv.gc <- read.table('CroVir_genome.GC.10kb.txt',header=T)
cvv.gc.a <- cvv.gc[which(cvv.gc$chrom!='scaffold-Z'),]
cvv.gc.z <- cvv.gc[which(cvv.gc$chrom=='scaffold-Z'),]
cvv.gc.z <- cvv.gc.z[which(cvv.gc.z$start<106800000),]

gal.gc <- read.table('genome_gallus.GC.10kb.txt',header=T)
# Gallus Z and W chromosome scaffolds are NC_006127.5 and NC_006126.5
gal.gc.a <- gal.gc[which(gal.gc$chrom!='NC_006127.5' & gal.gc$chrom!='NC_006126.5'),]
gal.gc.z <- gal.gc[which(gal.gc$chrom=='NC_006127.5'),]
gal.gc.w <- gal.gc[which(gal.gc$chrom=='NC_006126.5'),]

tae.gc <- read.table('genome_taeniopygia.GC.10kb.txt',header=T)
# Taeniopygia Z and W chromosome scaffolds are NC_045027.1 and NC_045028.1
tae.gc.a <- tae.gc[which(tae.gc$chrom!='NC_045027.1' & tae.gc$chrom!='NC_045028.1'),]
tae.gc.z <- tae.gc[which(tae.gc$chrom=='NC_045027.1'),]
tae.gc.w <- tae.gc[which(tae.gc$chrom=='NC_045028.1'),]

hom.gc <- read.table('genome_homo.GC.10kb.txt',header=T)
# Homo X and Y chromosome scaffolds are NC_000023.11 and NC_000024.10
hom.gc.a <- hom.gc[which(hom.gc$chrom!='NC_000023.11' & hom.gc$chrom!='NC_000024.10'),]
hom.gc.x <- hom.gc[which(hom.gc$chrom=='NC_000023.11'),]
hom.gc.y <- hom.gc[which(hom.gc$chrom=='NC_000024.10'),]

mus.gc <- read.table('genome_mus.GC.10kb.txt',header=T)
# Mus X and Y chromosome scaffolds are NC_000086.8 and NC_000087.8
mus.gc.a <- mus.gc[which(mus.gc$chrom!='NC_000086.8' & mus.gc$chrom!='NC_000087.8'),]
mus.gc.x <- mus.gc[which(mus.gc$chrom=='NC_000086.8'),]
mus.gc.y <- mus.gc[which(mus.gc$chrom=='NC_000087.8'),]

## CpG content

cvv.cpg.w <- read.table('Cviridis_CV0650_candidate_W.CpG.txt',header=T)

cvv.cpg <- read.table('CroVir_genome.CpG.10kb.txt',header=T)
cvv.cpg.a <- cvv.cpg[which(cvv.cpg$chrom!='scaffold-Z'),]
cvv.cpg.z <- cvv.cpg[which(cvv.cpg$chrom=='scaffold-Z'),]
cvv.cpg.z <- cvv.cpg.z[which(cvv.cpg.z$start<106800000),]

gal.cpg <- read.table('genome_gallus.CpG.10kb.txt',header=T)
# Gallus Z and W chromosome scaffolds are NC_006127.5 and NC_006126.5
gal.cpg.a <- gal.cpg[which(gal.cpg$chrom!='NC_006127.5' & gal.cpg$chrom!='NC_006126.5'),]
gal.cpg.z <- gal.cpg[which(gal.cpg$chrom=='NC_006127.5'),]
gal.cpg.w <- gal.cpg[which(gal.cpg$chrom=='NC_006126.5'),]

tae.cpg <- read.table('genome_taeniopygia.CpG.10kb.txt',header=T)
# Taeniopygia Z and W chromosome scaffolds are NC_045027.1 and NC_045028.1
tae.cpg.a <- tae.cpg[which(tae.cpg$chrom!='NC_045027.1' & tae.cpg$chrom!='NC_045028.1'),]
tae.cpg.z <- tae.cpg[which(tae.cpg$chrom=='NC_045027.1'),]
tae.cpg.w <- tae.cpg[which(tae.cpg$chrom=='NC_045028.1'),]

hom.cpg <- read.table('genome_homo.CpG.10kb.txt',header=T)
# Homo X and Y chromosome scaffolds are NC_000023.11 and NC_000024.10
hom.cpg.a <- hom.cpg[which(hom.cpg$chrom!='NC_000023.11' & hom.cpg$chrom!='NC_000024.10'),]
hom.cpg.x <- hom.cpg[which(hom.cpg$chrom=='NC_000023.11'),]
hom.cpg.y <- hom.cpg[which(hom.cpg$chrom=='NC_000024.10'),]

mus.cpg <- read.table('genome_mus.CpG.10kb.txt',header=T)
# Mus X and Y chromosome scaffolds are NC_000086.8 and NC_000087.8
mus.cpg.a <- mus.cpg[which(mus.cpg$chrom!='NC_000086.8' & mus.cpg$chrom!='NC_000087.8'),]
mus.cpg.x <- mus.cpg[which(mus.cpg$chrom=='NC_000086.8'),]
mus.cpg.y <- mus.cpg[which(mus.cpg$chrom=='NC_000087.8'),]

### Plot distributions of GC and CpG on autosomes and sex chromosomes-------

par(mfrow=c(2,5))
boxplot(cvv.gc.a$GC_content,cvv.gc.z$GC_content,cvv.gc.w$GC_content,outline=F,col=c('grey','seagreen','navy'),names=c('Auto','Z','W'),ylim=c(0,0.6),ylab='GC Content',main='prairie rattlesnake')
boxplot(tae.gc.a$GC_content,tae.gc.z$GC_content,tae.gc.w$GC_content,outline=F,col=c('grey','seagreen','navy'),names=c('Auto','Z','W'),ylim=c(0,0.6),ylab='GC Content',main='zebra finch')
boxplot(gal.gc.a$GC_content,gal.gc.z$GC_content,gal.gc.w$GC_content,outline=F,col=c('grey','seagreen','navy'),names=c('Auto','Z','W'),ylim=c(0,0.6),ylab='GC Content',main='chicken')
boxplot(hom.gc.a$GC_content,hom.gc.x$GC_content,hom.gc.y$GC_content,outline=F,col=c('grey','seagreen','navy'),names=c('Auto','X','Y'),ylim=c(0,0.6),ylab='GC Content',main='human')
boxplot(mus.gc.a$GC_content,mus.gc.x$GC_content,mus.gc.y$GC_content,outline=F,col=c('grey','seagreen','navy'),names=c('Auto','X','Y'),ylim=c(0,0.6),ylab='GC Content',main='house mouse')

boxplot(cvv.cpg.a$CpG_content,cvv.cpg.z$CpG_content,cvv.cpg.w$CpG_content,outline=F,col=c('grey','seagreen','navy'),names=c('Auto','Z','W'),ylim=c(0,0.07),ylab='CpG Content',main='prairie rattlesnake')
boxplot(tae.cpg.a$CpG_content,tae.cpg.z$CpG_content,tae.cpg.w$CpG_content,outline=F,col=c('grey','seagreen','navy'),names=c('Auto','Z','W'),ylim=c(0,0.07),ylab='CpG Content',main='zebra finch')
boxplot(gal.cpg.a$CpG_content,gal.cpg.z$CpG_content,gal.cpg.w$CpG_content,outline=F,col=c('grey','seagreen','navy'),names=c('Auto','Z','W'),ylim=c(0,0.07),ylab='CpG Content',main='chicken')
boxplot(hom.cpg.a$CpG_content,hom.cpg.x$CpG_content,hom.cpg.y$CpG_content,outline=F,col=c('grey','seagreen','navy'),names=c('Auto','X','Y'),ylim=c(0,0.07),ylab='CpG Content',main='human')
boxplot(mus.cpg.a$CpG_content,mus.cpg.x$CpG_content,mus.cpg.y$CpG_content,outline=F,col=c('grey','seagreen','navy'),names=c('Auto','X','Y'),ylim=c(0,0.07),ylab='CpG Content',main='house mouse')

### Summary statistics of GC and CpG distributions--------------------------

mean(cvv.gc.a$GC_content)
mean(tae.gc.a$GC_content)
mean(gal.gc.a$GC_content)
mean(hom.gc.a$GC_content)
mean(mus.gc.a$GC_content)

sd(cvv.gc.a$GC_content)
sd(tae.gc.a$GC_content)
sd(gal.gc.a$GC_content)
sd(hom.gc.a$GC_content)
sd(mus.gc.a$GC_content)

mean(cvv.gc.z$GC_content)
mean(tae.gc.z$GC_content)
mean(gal.gc.z$GC_content)
mean(hom.gc.x$GC_content)
mean(mus.gc.x$GC_content)

sd(cvv.gc.z$GC_content)
sd(tae.gc.z$GC_content)
sd(gal.gc.z$GC_content)
sd(hom.gc.x$GC_content)
sd(mus.gc.x$GC_content)

mean(cvv.gc.w$GC_content)
mean(tae.gc.w$GC_content)
mean(gal.gc.w$GC_content)
mean(hom.gc.y$GC_content)
mean(mus.gc.y$GC_content)

sd(cvv.gc.w$GC_content)
sd(tae.gc.w$GC_content)
sd(gal.gc.w$GC_content)
sd(hom.gc.y$GC_content)
sd(mus.gc.y$GC_content)

### Statistical comparisons of GC & CpG distributions-----------------------

## GC & CpG - Shapiro-Wilk normality tests
shapiro.test(cvv.gc.a$GC_content[1:5000])
shapiro.test(cvv.gc.w$GC_content)

shapiro.test(cvv.cpg.a$CpG_content[1:5000])
shapiro.test(cvv.cpg.w$CpG_content)

## GC - Mann-Whitney U tests (not normally distributed)
wilcox.test(cvv.gc.a$GC_content,cvv.gc.w$GC_content)
wilcox.test(cvv.gc.z$GC_content,cvv.gc.w$GC_content)
wilcox.test(tae.gc.a$GC_content,tae.gc.w$GC_content)
wilcox.test(tae.gc.z$GC_content,tae.gc.w$GC_content)
wilcox.test(gal.gc.a$GC_content,gal.gc.w$GC_content)
wilcox.test(gal.gc.z$GC_content,gal.gc.w$GC_content)
wilcox.test(hom.gc.a$GC_content,hom.gc.y$GC_content)
wilcox.test(hom.gc.x$GC_content,hom.gc.y$GC_content)
wilcox.test(mus.gc.a$GC_content,mus.gc.y$GC_content)
wilcox.test(mus.gc.x$GC_content,mus.gc.y$GC_content)

## CpG - Mann-Whitney U tests (not normally distributed)
wilcox.test(cvv.cpg.a$CpG_content,cvv.cpg.w$CpG_content)
wilcox.test(cvv.cpg.z$CpG_content,cvv.cpg.w$CpG_content)
wilcox.test(tae.cpg.a$CpG_content,tae.cpg.w$CpG_content)
wilcox.test(tae.cpg.z$CpG_content,tae.cpg.w$CpG_content)
wilcox.test(gal.cpg.a$CpG_content,gal.cpg.w$CpG_content)
wilcox.test(gal.cpg.z$CpG_content,gal.cpg.w$CpG_content)
wilcox.test(hom.cpg.a$CpG_content,hom.cpg.y$CpG_content)
wilcox.test(hom.cpg.x$CpG_content,hom.cpg.y$CpG_content)
wilcox.test(mus.cpg.a$CpG_content,mus.cpg.y$CpG_content)
wilcox.test(mus.cpg.x$CpG_content,mus.cpg.y$CpG_content)

## Test for equal variances using Levene's tests
leveneTest(cvv.gc.a$GC_content[1:2139],cvv.gc.z$GC_content[1:2139])
leveneTest(cvv.gc.w$GC_content[1:2139],cvv.gc.z$GC_content[1:2139])
leveneTest(cvv.gc.w$GC_content[1:2139],cvv.gc.a$GC_content[1:2139])

## Test for equal variances using Levene's tests
leveneTest(cvv.cpg.a$CpG_content[1:2139],cvv.cpg.z$CpG_content[1:2139])
leveneTest(cvv.cpg.w$CpG_content[1:2139],cvv.cpg.z$CpG_content[1:2139])
leveneTest(cvv.cpg.w$CpG_content[1:2139],cvv.cpg.a$CpG_content[1:2139])

