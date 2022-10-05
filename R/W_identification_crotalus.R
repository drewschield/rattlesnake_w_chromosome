############################################################################
# Comparative mapping study in Crotalus
############################################################################

### Goal: compare distributions of male and female read coverage, and
### identify putative W-linked scaffolds in the female assembly.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

setwd('./W_chromosome_identification/')
library(scales)
library(zoo)
library(Rmisc)
install.packages("car")
library(car)

### Read in data-------------------------------------------------------------

### These data are 'global' (i.e., contig-wide) mean read depths for female
### and male resquencing data against female assembly pseudohaplotypes (1 and 2).

m1 <- read.table('./coverage_exp_crotalus/CV0007_CV0650_pseudohap1.mosdepth.summary.txt',header=T)
f1 <- read.table('./coverage_exp_crotalus/CV0011_CV0650_pseudohap1.mosdepth.summary.txt',header=T)
## Remove last 'total' row
m1 <- head(m1,-1)
f1 <- head(f1,-1)

# Calculate mean of means for both sexes for average depth reporting in paper!
mean(m1$mean)
mean(f1$mean)

m2 <- read.table('./coverage_exp_crotalus/CV0007_CV0650_pseudohap2.mosdepth.summary.txt',header=T)
f2 <- read.table('./coverage_exp_crotalus/CV0011_CV0650_pseudohap2.mosdepth.summary.txt',header=T)
## Remove last 'total' row
m2 <- head(m2,-1)
f2 <- head(f2,-1)

### Identify scaffolds not in male mapping (i.e., potential W)--------------

only.fem1 <- setdiff(f1$chrom,m1$chrom)
write.table(only.fem1,'./pseudohap1_female-only_scaffolds.txt',quote=F,row.names=FALSE)

only.fem2 <- setdiff(f2$chrom,m2$chrom)
write.table(only.fem2,'./pseudohap2_female-only_scaffolds.txt',quote=F,row.names=FALSE)

### Intersect male and female data frames by scaffolds with data

common1 <- merge(m1,f1, by.x = "chrom", by.y = "chrom")
write.table(common1,'./pseudohap1_all_common_scaffolds.txt',quote=F,row.names=FALSE)

common2 <- merge(m2,f2, by.x = "chrom", by.y = "chrom")
write.table(common2, './pseudohap2_all_common_scaffolds.txt',quote=F,row.names=FALSE)

# 'x' and 'y' in joined common data equal 'male' and 'female', respectively.

### Basic plotting-----------------------------------------------------------

hist(common1$mean.x,breaks=100000,xlim=c(0,100))
hist(common1$mean.y,breaks=100000,xlim=c(0,100))

hist(common2$mean.x,breaks=100000,xlim=c(0,100))
hist(common2$mean.y,breaks=100000,xlim=c(0,100))

### Normalized ratios--------------------------------------------------------

### Pseudohaplotype 1 results:

## Set median per sex
mmed <- median(common1$mean.x)
fmed <- median(common1$mean.y)
## Normalize by median per sex
mnorm <- common1$mean.x/mmed
fnorm <- common1$mean.y/fmed

## Calculate ratio of normalized coverage; output scaffolds with log2 F:M >= 1
fmnorm <- log2(fnorm/mnorm)
plot(fmnorm)
hist(fmnorm,breaks=1000,xlim=c(1,10),ylim=c(0,500))
count(fmnorm>=1)

common1["fmnorm"] <- fmnorm

fmnorm.thresh <- common1[which(common1$fmnorm>=1),]
write.table(fmnorm.thresh,'./pseudohap1_log2FM_thresh_scaffolds.txt',sep='\t',quote=F,row.names=FALSE)

sum(fmnorm.thresh$length.x)


## Calculate proportion of female coverage; output scaffolds > Q3 + 1.5*IQR
# Calculate total normalized coverage per scaffold
tnorm <- mnorm + fnorm
# Calculate proportion per sex of normalized coverage
fpropnorm <- fnorm/tnorm
mpropnorm <- mnorm/tnorm
# Calculate quantiles and interquartile range for female proportion
quantile(fpropnorm,c(0.25,0.5,0.75))
IQR(fpropnorm)
# Q3 = 0.5155, IQR = 0.032
thresh <- 0.5155 + (1.5*0.0132)
# Determine how many scaffolds have female-biased coverage
count(fpropnorm >= thresh)

# Basic plots for normalized proportions
hist(fpropnorm,breaks=100)
hist(mpropnorm,breaks=100)
plot(fpropnorm,mpropnorm)

# Add female normalized proportion to common dataframe
common1["fpropnorm"] <- fpropnorm
write.table(common1,'./pseudohap1_all_common_scaffolds_female_cov_proportion.txt',quote=F,sep='\t')

# Extract table of scaffolds with female proportion of coverage > Q3 + 1.5*IQR
fpropnorm.thresh <- common1[which(common1$fpropnorm>=thresh),]
write.table(fpropnorm.thresh,'./pseudohap1_IQR_thresh_scaffolds.txt',sep='\t',quote=F)

sum(fpropnorm.thresh$length.x)


### Pseudohaplotype 2 results:

## Set median per sex
mmed <- median(common2$mean.x)
fmed <- median(common2$mean.y)
## Normalize by median per sex
mnorm <- common2$mean.x/mmed
fnorm <- common2$mean.y/fmed

fmnorm <- log2(fnorm/mnorm)
plot(fmnorm)
hist(fmnorm,breaks=1000,xlim=c(1,10),ylim=c(0,500))
count(fmnorm>=1)

common2["fmnorm"] <- fmnorm

fmnorm.thresh <- common2[which(common2$fmnorm>=1),]
write.table(fmnorm.thresh,'./pseudohap2_log2FM_thresh_scaffolds.txt',sep='\t',quote=F,row.names=FALSE)

sum(fmnorm.thresh$length.x)

## Calculate total normalized
tnorm <- mnorm + fnorm
## Calculate proportion per sex of normalized coverage
fpropnorm <- fnorm/tnorm
mpropnorm <- mnorm/tnorm
## Calculate quantiles and interquartile range for female proportion
quantile(fpropnorm,c(0.25,0.5,0.75))
IQR(fpropnorm)
### Q3 = 0.5155, IQR = 0.032
thresh <- 0.5155 + (1.5*0.032)
## Determine how many scaffolds have female-biased coverage
count(fpropnorm >= thresh)

hist(fpropnorm)
hist(mpropnorm)
plot(fpropnorm,mpropnorm)

## Add female normalized proportion to common dataframe
common2["fpropnorm"] <- fpropnorm
write.table(common2,'./pseudohap2_all_common_scaffolds_female_cov_proportion.txt',quote=F,sep='\t')

# Extract table of scaffolds with female proportion of coverage > Q3 + 1.5*IQR
fpropnorm.thresh <- common2[which(common2$fpropnorm>=thresh),]
write.table(fpropnorm.thresh,'./pseudohap2_IQR_thresh_scaffolds.txt',sep='\t',quote=F)

sum(fpropnorm.thresh$length.x)
