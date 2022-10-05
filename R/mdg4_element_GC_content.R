############################################################################
# Analysis of mdg4 element GC content on the W chromosome
############################################################################

### Goal: Examine the distribution of GC content for known snake mdg4
### elements used to annotate the prairie rattlesnake W chromosome. Compare
### the density of mdg4 elements and GC content for W-linked scaffolds.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./repeats')

install.packages("ppcor")
library(ppcor)
library(scales)

### GC content of known mdg4 elements compared to other prevalent TEs-------

## These are data from the 'known' TE sequences in the snake repeat library
## used to annotate repeats on the W chromosome.

## Read in data

gc.rep <- read.table('Snakes_Known_TElib.GC.txt',header=T)
gc.mdg <- read.table('Snakes_Known_TElib.mdg4.GC.txt',header=T)

## Note: comparisons with L1 and CR1 elements can be done by parsing the
## known elements for a specific family, e.g.,:

## gc.cr1 <- gc.rep[which("CR1" in gc.rep$chrom),]

### Basic statistics

mean(gc.rep$GC_content)
sd(gc.rep$GC_content)

mean(gc.mdg$GC_content)
sd(gc.mdg$GC_content)

### Basic plotting

boxplot(gc.rep$GC_content,gc.mdg$GC_content,outline=F,
        names=c('All','mdg4'))


### GC content of annotated mdg4 elements compared to other prevalent TEs---

## These are data for actual sequences from annotated elements on the W.

## Read in data

gc.rep.seq <- read.table('Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_all.GC.txt',header=T)
gc.mdg.seq <- read.table('Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_mdg4.GC.txt',header=T)
gc.l1.seq <- read.table('Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_L1-CIN4.GC.txt',header=T)
gc.l2.seq <- read.table('Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_L2-CR1-Rex.GC.txt',header=T)

gc.nomdg.seq <- read.table('Cviridis_CV0650_candidate_W.rescaffold.rename.full_mask.reformat.repeats_no-mdg4.GC.txt',header=T)

## Basic statistics

mean(gc.rep.seq$GC_content)
sd(gc.rep.seq$GC_content)

mean(gc.mdg.seq$GC_content)
sd(gc.mdg.seq$GC_content)

mean(gc.l1.seq$GC_content)
sd(gc.l1.seq$GC_content)

mean(gc.l2.seq$GC_content)
sd(gc.l2.seq$GC_content)

mean(gc.nomdg.seq$GC_content)
sd(gc.nomdg.seq$GC_content)

## Statistical comparison of GC content

# Test of normality
shapiro.test(gc.mdg.seq$GC_content)
# Not normally distributed, use Mann-whitney U tests

wilcox.test(gc.mdg.seq$GC_content,gc.l1.seq$GC_content)
wilcox.test(gc.mdg.seq$GC_content,gc.l2.seq$GC_content)
wilcox.test(gc.mdg.seq$GC_content,gc.nomdg.seq$GC_content)

## Basic plotting

boxplot(gc.nomdg.seq$GC_content,gc.l1.seq$GC_content,gc.l2.seq$GC_content,gc.mdg.seq$GC_content,outline=F,
        names=c('All','L1','L2','mdg4'))

### Relationship between GC content and mdg4 element density on the W-------

## Read in data
w.mdg <- read.table('Cviridis_CV0650_candidate_W.repeat_content_mdg4.txt',header=T)
w.gc <- read.table('../gc/Cviridis_CV0650_candidate_W.rescaffold.rename.GC.txt',header=T)
w.cpg <- read.table('../gc/Cviridis_CV0650_candidate_W.rescaffold.rename.CpG.txt',header=T)

## Basic plotting

par(mfrow=c(2,1))
plot(w.mdg$prop_repeat,w.gc$GC_content,pch=20,col=alpha('grey25',0.25))
abline(lm(w.gc$GC_content~w.mdg$prop_repeat),lty=2,lwd=2,col='navy')
plot(w.mdg$prop_mdg4,w.gc$GC_content,pch=20,col=alpha('grey25',0.25))
abline(lm(w.gc$GC_content~w.mdg$prop_mdg4),lty=2,lwd=2,col='navy')

## Correlation coefficients

cor.test(w.mdg4$prop_mdg4,w.gc$GC_content,method='spearman')
cor.test(w.mdg4$prop_repeat,w.gc$GC_content,method='spearman')

## Perform partial correlation, controlling for CpG & total repeats
comp.data <- data.frame(w.mdg4$prop_mdg4,w.gc$GC_content,w.mdg$prop_repeat,w.cpg$CpG_content)
pcor(comp.data, method = 'spearman')

## Extract scaffolds lacking Gypsy elements

w.no.mdg <- w.mdg[which(w.mdg$prop_mdg4==0),]
w.no.gc <- w.gc[which(w.gyp$prop_mdg4==0),]

## Plot scatterplots without/with mdg4 elements

par(mfrow=c(2,1))
plot(w.no.mdg$prop_repeat,w.no.gc$GC_content,pch=20,col=alpha('grey25',0.25),xlab='% repeats (no mdg4)',ylab='GC content')
abline(lm(w.no.gc$GC_content~w.no.mdg$prop_repeat),lty=2,lwd=2,col='navy')
plot(w.mdg$prop_mdg4,w.gc$GC_content,pch=20,col=alpha('grey25',0.25),xlab='% mdg4',ylab='GC content')
abline(lm(w.gc$GC_content~w.mdg$prop_mdg4),lty=2,lwd=2,col='navy')

cor.test(w.no.mdg$prop_repeat,w.no.gc$GC_content,method='spearman')
cor.test(w.mdg$prop_mdg4,w.gc$GC_content,method='spearman')
