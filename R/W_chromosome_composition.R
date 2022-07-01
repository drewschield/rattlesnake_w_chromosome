############################################################################
# W Chromosome composition and comparison to Z and autosomes
############################################################################

### Goal: examine measures of composition (GC, CpG, repeat content) and make
### comparisons between the W, Z, and Autosomes.

### Working directory and dependencies--------------------------------------

setwd('./')

install.packages("ppcor")
library(ppcor)

### Read in data------------------------------------------------------------

# GC content 
cvv.gc.w <- read.table('./gc/Cviridis_CV0650_candidate_W.GC.txt',header=T)

cvv.gc <- read.table('./gc/CroVir_genome.GC.10kb.txt',header=T)
cvv.gc.a <- cvv.gc[which(cvv.gc$chrom!='scaffold-Z'),]
cvv.gc.z <- cvv.gc[which(cvv.gc$chrom=='scaffold-Z'),]
cvv.gc.z <- cvv.gc.z[which(cvv.gc.z$start<106800000),]

# CpG content

cvv.cpg.w <- read.table('./gc/Cviridis_CV0650_candidate_W.CpG.txt',header=T)

cvv.cpg <- read.table('./gc/CroVir_genome.CpG.10kb.txt',header=T)
cvv.cpg.a <- cvv.cpg[which(cvv.cpg$chrom!='scaffold-Z'),]
cvv.cpg.z <- cvv.cpg[which(cvv.cpg$chrom=='scaffold-Z'),]
cvv.cpg.z <- cvv.cpg.z[which(cvv.cpg.z$start<106800000),]

# Repeat content

cvv.rep.w <- read.table('./repeats/Cviridis_CV0650_candidate_W.repeat_content_mdg4.txt',header=T)

cvv.rep <- read.table('./repeats/CroVir_genome.mdg4.10kb.txt',header=T)
cvv.rep.a <- cvv.rep[which(cvv.rep$chrom!='scaffold-Z'),]
cvv.rep.z <- cvv.rep[which(cvv.rep$chrom=='scaffold-Z'),]
cvv.rep.z <- cvv.rep.z[which(cvv.rep.z$start<106800000),]

# GC content of snake repeat library specific elements
snake.mdg <- read.table('./repeats/snake_repeat_library/Snakes_Known_TElib.GC.mdg4.txt',header=T)
snake.l1 <- read.table('./repeats/snake_repeat_library/Snakes_Known_TElib.GC.L1-CIN4.filter.txt',header=T)
snake.CR1 <- read.table('./repeats/snake_repeat_library/Snakes_Known_TElib.GC.L2-CR1-Rex.filter.txt',header=T)

# Repeat ages (from RpeatMasker output)
z.age <- read.csv('./repeats/CroVir_Z_18Snake.csv',header=T)
w.age <- read.csv('./repeats/CroVir_W_18Snake.csv',header=T)

### Analysis of GC content---------------------------------------------------

mean(cvv.gc.w$GC_content)
sd(cvv.gc.w$GC_content)

mean(cvv.gc.z$GC_content)
sd(cvv.gc.z$GC_content)

mean(cvv.gc.a$GC_content)
sd(cvv.gc.a$GC_content)

# Does W have higher GC content than the Z?
wilcox.test(cvv.gc.w$GC_content,cvv.gc.z$GC_content)
# Yes.

# Does W have higher GC content than autosomes?
wilcox.test(cvv.gc.w$GC_content,cvv.gc.a$GC_content)
# Yes.

# Plot density distributions
par(mfrow=c(1,2))
plot(density(cvv.gc.a$GC_content),col='grey',main=NA,xlab='GC Content')
polygon(density(cvv.gc.a$GC_content),col=alpha('grey',0.6))
lines(density(cvv.gc.z$GC_content),col='seagreen')
polygon(density(cvv.gc.z$GC_content),col=alpha('seagreen',0.6))
lines(density(cvv.gc.w$GC_content),col='navy')
polygon(density(cvv.gc.w$GC_content),col=alpha('navy',0.6))

# Or, make a boxplot
boxplot(cvv.gc.a$GC_content,cvv.gc.z$GC_content,cvv.gc.w$GC_content,medcol=c('grey','seagreen','navy'),lwd=2,names=c('Autosomes','Z','W'),pch=20,ylab='GC Content')
boxplot(cvv.gc.a$GC_content,cvv.gc.z$GC_content,cvv.gc.w$GC_content,medcol='black',col=alpha(c('grey','seagreen','navy'),0.7),lwd=1.5,names=c('Autosomes','Z','W'),pch=20,ylab='GC Content')

### Analysis of CpG content---------------------------------------------------

mean(cvv.cpg.w$CpG_content)
sd(cvv.cpg.w$CpG_content)

mean(cvv.cpg.z$CpG_content)
sd(cvv.cpg.z$CpG_content)

mean(cvv.cpg.a$CpG_content)
sd(cvv.cpg.a$CpG_content)

# Does W have higher CpG content than the Z?
t.test(cvv.cpg.w$CpG_content,cvv.cpg.z$CpG_content)
# Yes.

# Does W have higher CpG content than autosomes?
t.test(cvv.cpg.w$CpG_content,cvv.cpg.a$CpG_content)
# Yes.

# Plot density distributions
par(mfrow=c(1,1))
plot(density(cvv.cpg.a$CpG_content),col='grey',main=NA,xlab='CpG Content')
polygon(density(cvv.cpg.a$CpG_content),col=alpha('grey',0.5))
lines(density(cvv.cpg.z$CpG_content),col='seagreen')
polygon(density(cvv.cpg.z$CpG_content),col=alpha('seagreen',0.5))
lines(density(cvv.cpg.w$CpG_content),col='darkblue')
polygon(density(cvv.cpg.w$CpG_content),col=alpha('darkblue',0.5))

# Or, make a boxplot
boxplot(cvv.cpg.a$CpG_content,cvv.cpg.z$CpG_content,cvv.cpg.w$CpG_content,medcol=c('grey','seagreen','navy'),lwd=2,names=c('Autosomes','Z','W'),pch=20,ylab='CpG Content')
boxplot(cvv.cpg.a$CpG_content,cvv.cpg.z$CpG_content,cvv.cpg.w$CpG_content,medcol='black',col=alpha(c('grey','seagreen','navy'),0.7),lwd=1.5,names=c('Autosomes','Z','W'),pch=20,ylab='CpG Content')

### Plot distributions of GC and CpG content-------------------------------------

par(mfrow=c(1,2))
boxplot(cvv.gc.a$GC_content,cvv.gc.z$GC_content,cvv.gc.w$GC_content,medcol='black',col=alpha(c('grey','seagreen','navy'),0.7),lwd=1.5,names=c('Autosomes','Z','W'),pch=20,ylab='GC Content')
boxplot(cvv.cpg.a$CpG_content,cvv.cpg.z$CpG_content,cvv.cpg.w$CpG_content,medcol='black',col=alpha(c('grey','seagreen','navy'),0.7),lwd=1.5,names=c('Autosomes','Z','W'),pch=20,ylab='CpG Content')

# Great job.

### Analysis of repeat elements--------------------------------------------------

# UPDATE/REFORMAT WHILE WORKING ON REPEAT SECTIONS ON GITHUB
# UPDATE/REFORMAT WHILE WORKING ON REPEAT SECTIONS ON GITHUB
# UPDATE/REFORMAT WHILE WORKING ON REPEAT SECTIONS ON GITHUB

mean(cvv.rep.w$prop_Gypsy)
sd(cvv.rep.w$prop_Gypsy.DIRS)

mean(cvv.rep.z$DIRS.Gypsy_proportion)
sd(cvv.rep.z$DIRS.Gypsy_proportion)

mean(cvv.rep.a$DIRS.Gypsy_proportion)
sd(cvv.rep.a$DIRS.Gypsy_proportion)

# Does W have higher DIRS/Gypsy content than the Z?
wilcox.test(cvv.rep.w$prop_Gypsy.DIRS,cvv.rep.z$DIRS.Gypsy_proportion)
# Holy smokes, yes.

# Does W have higher DIRS/Gypsy content than autosomes?
wilcox.test(cvv.rep.w$prop_Gypsy.DIRS,cvv.rep.a$DIRS.Gypsy_proportion)
# You bet your ass it does.

# Does Z have higher DIRS/Gypsy content than autosomes?
wilcox.test(cvv.rep.z$DIRS.Gypsy_proportion,cvv.rep.a$DIRS.Gypsy_proportion)
# It sure does.

# Make a boxplot
boxplot(cvv.rep.a$DIRS.Gypsy_proportion,cvv.rep.z$DIRS.Gypsy_proportion,cvv.rep.w$prop_Gypsy.DIRS,medcol='black',col=alpha(c('grey','seagreen','navy'),0.7),lwd=1.5,names=c('Autosomes','Z','W'),pch=20,ylab='Proportion DIRS/Gypsy')

# Holy shit the W is nuts.

### Comparison of GC content and DIRS/Gypsy elements------------------------------

# Is GC content on the W associated with repeats overall?
plot(cvv.gc.w$GC_content,cvv.rep.w$prop_repeat,pch=20)
cor.test(cvv.gc.w$GC_content,cvv.rep.w$prop_repeat, method='spearman')
# Yes, weakly.

# Is GC content on the W associated with Gypsy elements?
par(mfrow=c(1,1))
plot(cvv.gc.w$GC_content,cvv.rep.w$prop_Gypsy,pch=20,xlab='GC Content',ylab='Percent Gypsy elements',col=alpha('navy',0.5))
abline(lm(cvv.rep.w$prop_Gypsy~cvv.gc.w$GC_content),lty=2,lwd=2,col='red')
cor.test(cvv.gc.w$GC_content,cvv.rep.w$prop_Gypsy, method='spearman')
# Yes.

# test partial correlation of GC content and Gypsy elements, controlling for CpG and total repeat content
comp.data <- data.frame(cvv.gc.w$GC_content,cvv.cpg.w$CpG_content,cvv.rep.w$prop_repeat,cvv.rep.w$prop_Gypsy)
pcor(comp.data, method = 'spearman')

# Examine GC content versus DIRS/Gypsy on Z and autosomes
par(mfrow=c(1,3))
plot(cvv.gc.a$GC_content,cvv.rep.a$DIRS.Gypsy_proportion,pch=20,xlab='GC Content',ylab='Percent Gypsy elements',col=alpha('grey50',0.5),xlim=c(0.15,0.65),ylim=c(0,0.6))
text(x=0.25,y=0.55,"Autosome",col='grey50')
plot(cvv.gc.z$GC_content,cvv.rep.z$DIRS.Gypsy_proportion,pch=20,xlab='GC Content',ylab='Percent Gypsy/DIRS elements',col=alpha('seagreen',0.5),xlim=c(0.15,0.65),ylim=c(0,0.6))
text(x=0.25,y=0.55,"Z Chromosome",col='seagreen')
plot(cvv.gc.w$GC_content,cvv.rep.w$prop_Gypsy,pch=20,xlab='GC Content',ylab='Percent Gypsy/DIRS elements',col=alpha('navy',0.5),xlim=c(0.15,0.65),ylim=c(0,0.6))
text(x=0.25,y=0.55,"W Chromosome",col='navy')


cor.test(cvv.gc.z$GC_content,cvv.rep.z$DIRS.Gypsy_proportion,method='spearman')
cor.test(cvv.gc.a$GC_content,cvv.rep.a$DIRS.Gypsy_proportion,method='spearman')

### Analysis of GC content of specific elements from snake repeat library----------

mean(snake.gyp$GC_content)
sd(snake.gyp$GC_content)

mean(snake.l1$GC_content)
sd(snake.l1$GC_content)

mean(snake.CR1$GC_content)
sd(snake.CR1$GC_content)

t.test(snake.gyp$GC_content,snake.l1$GC_content)
t.test(snake.gyp$GC_content,snake.CR1$GC_content)

### Analysis of Z and W repeat age distributions-----------------------------------

par(mfrow=c(2,1))
plot(w.age$Div,w.age$LTR_Gypsy,pch=20,ylim=c(0,1))
plot(z.age$Div,z.age$LTR_Gypsy,pch=20,ylim=c(0,1))

plot(w.age$LTR_Gypsy,z.age$LTR_Gypsy,xlim=c(0,1),ylim=c(0,1),pch=20)
cor.test(w.age$LTR_Gypsy,z.age$LTR_Gypsy)

w.sum <- sum(w.age$LTR_Gypsy)
z.sum <- sum(z.age$LTR_Gypsy)

w.weight <- (w.age$LTR_Gypsy/w.sum)
z.weight <- (z.age$LTR_Gypsy/z.sum)

plot(w.age$Div,w.weight)
plot(z.age$Div,z.weight)

plot(w.weight,z.weight,pch=20)
cor.test(w.weight,z.weight)

t.test(w.weight,z.weight)
