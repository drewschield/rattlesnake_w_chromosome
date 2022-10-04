############################################################################
# Comparative mapping study in other species
############################################################################

### Goal: compare distributions of male and female read coverage, and
### identify putative W-linked scaffolds in the female assembly.

### Working directory and dependencies--------------------------------------

setwd('./comparative_W_coverage/gallus/')

install.packages("car")
install.packages("Rmisc")
install.packages("pbkrtest")
library(scales)
library(zoo)
library(Rmisc)
library(car)

### -------------------------------------------------------------------------
###                             Gallus
### -------------------------------------------------------------------------

### Read in data-------------------------------------------------------------

### These data are 'global' (i.e., contig-wide) mean read depths for female
### and male resquencing data against the Gallus gallus v6 genome.

gal.m <- read.table('gallus_male_GRCg6a.10kb.regions.bed',header=T)
gal.f <- read.table('gallus_female_GRCg6a.10kb.regions.bed',header=T)

### Intersect male and female data frames by scaffolds with data

gal.common <- merge(gal.m,gal.f, by=c("chrom","start","end"))
write.table(gal.common,'./gallus_all_common_scaffolds.txt',sep='\t',quote=F)

### Histograms---------------------------------------------------------------

hist(gal.common$mean.x,breaks=100000,xlim=c(0,100))
hist(gal.common$mean.y,breaks=100000,xlim=c(0,100))

### Log2 Normalized ratios---------------------------------------------------

## Set median per sex
gal.mmed <- median(gal.common$mean.x)
gal.fmed <- median(gal.common$mean.y)
## Normalize by median per sex
gal.mnorm <- gal.common$mean.x/gal.mmed
gal.fnorm <- gal.common$mean.y/gal.fmed

## Calculate ratio of normalized coverage; output scaffolds with log2 F:M >= 1
gal.fmnorm <- log2(gal.fnorm/gal.mnorm)
plot(gal.fmnorm)
hist(gal.fmnorm,breaks=1000,xlim=c(1,10),ylim=c(0,500))
count(gal.fmnorm>=1)
count(gal.fmnorm>=0.75)
gal.common["fmnorm"] <- gal.fmnorm

gal.fmnorm.thresh <- gal.common[which(gal.common$fmnorm>=1),]
write.table(gal.fmnorm.thresh,'./gallus_log2FM_thresh_scaffolds.txt',sep='\t',quote=F)

sum(gal.fmnorm.thresh$end-gal.fmnorm.thresh$start)

### Calculate proportion of female coverage; output scaffolds > Q3 + 1.5*IQR

## Calculate total normalized coverage per scaffold
gal.tnorm <- gal.mnorm + gal.fnorm
## Calculate proportion per sex of normalized coverage
gal.fpropnorm <- gal.fnorm/gal.tnorm
gal.mpropnorm <- gal.mnorm/gal.tnorm
# Calculate quantiles and interquartile range for female proportion
quantile(gal.fpropnorm,c(0.25,0.5,0.75),na.rm=T)
IQR(gal.fpropnorm,na.rm=T)
# Q3 = 0.5141, IQR = 0.0212
gal.IQR.thresh <- 0.5141 + (1.5*0.0212)
# Determine how many scaffolds have female-biased coverage
count(gal.fpropnorm >= gal.IQR.thresh)

# Basic plots for normalized proportions
hist(gal.fpropnorm,breaks=100)
hist(gal.mpropnorm,breaks=100)
plot(gal.fpropnorm,gal.mpropnorm)

# Add female normalized proportion to common dataframe
gal.common["fpropnorm"] <- gal.fpropnorm
write.table(gal.common,'./gallus_all_common_scaffolds_female_cov_proportion.txt',quote=F,sep='\t')

# Extract table of scaffolds with female proportion of coverage > Q3 + 1.5*IQR
gal.fpropnorm.thresh <- gal.common[which(gal.common$fpropnorm>=gal.IQR.thresh),]
write.table(gal.fpropnorm.thresh,'./gallus_IQR_thresh_scaffolds.txt',sep='\t',quote=F)

sum(gal.fpropnorm.thresh$end - gal.fpropnorm.thresh$start)
