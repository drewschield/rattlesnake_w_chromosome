############################################################################
# Comparative recent stratum coverage in colubroids
############################################################################

### Goal: compare relative female coverage in thamnophis to patterns in
### pitviper species in the recent stratum to test whether recombination
### suppression occurred since the split of viperids and colubrids.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Load dependencies-------------------------------------------------------

library(Rmisc)
library(scales)

### Read in data------------------------------------------------------------

tham.m <- read.table('./coverage/mosdepth_results/thamnophis_mDNA.regions.bed',header=F)
tham.f <- read.table('./coverage/mosdepth_results/thamnophis_fDNA.regions.bed',header=F)
dein.m <- read.table('./coverage/mosdepth_results/deinagkistrodon_mDNA.regions.bed',header=F)
dein.f <- read.table('./coverage/mosdepth_results/deinagkistrodon_fDNA.regions.bed',header=F)
sist.m <- read.table('./coverage/mosdepth_results/sistrurus_mDNA.regions.bed',header=F)
sist.f <- read.table('./coverage/mosdepth_results/sistrurus_fDNA.regions.bed',header=F)
crot.m <- read.table('./coverage/mosdepth_results/crotalus_mDNA.regions.bed',header=F)
crot.f <- read.table('./coverage/mosdepth_results/crotalus_fDNA.regions.bed',header=F)

### Calculate log2FM per species--------------------------------------------

# The last command per block appends the log2FM values to the female data
# frame.

tham.m.a <- tham.m[which(tham.m$V1!="scaffold-Z"),]
tham.f.a <- tham.f[which(tham.f$V1!="scaffold-Z"),]
tham.m.z <- tham.m[which(tham.m$V1=="scaffold-Z"),]
tham.f.z <- tham.f[which(tham.f$V1=="scaffold-Z"),]
tham.m.a.med <- median(tham.m.a$V4,na.rm=T)
tham.f.a.med <- median(tham.f.a$V4,na.rm=T)
tham.m.norm <- tham.m.z$V4/tham.m.a.med
tham.f.norm <- tham.f.z$V4/tham.f.a.med
tham.fm.norm <- log2(tham.f.norm/tham.m.norm)
tham.f.z$fmnorm <- tham.fm.norm

dein.m.a <- dein.m[which(dein.m$V1!="scaffold-Z"),]
dein.f.a <- dein.f[which(dein.f$V1!="scaffold-Z"),]
dein.m.z <- dein.m[which(dein.m$V1=="scaffold-Z"),]
dein.f.z <- dein.f[which(dein.f$V1=="scaffold-Z"),]
dein.m.a.med <- median(dein.m.a$V4,na.rm=T)
dein.f.a.med <- median(dein.f.a$V4,na.rm=T)
dein.m.norm <- dein.m.z$V4/dein.m.a.med
dein.f.norm <- dein.f.z$V4/dein.f.a.med
dein.fm.norm <- log2(dein.f.norm/dein.m.norm)
dein.f.z$fmnorm <- dein.fm.norm

sist.m.a <- sist.m[which(sist.m$V1!="scaffold-Z"),]
sist.f.a <- sist.f[which(sist.f$V1!="scaffold-Z"),]
sist.m.z <- sist.m[which(sist.m$V1=="scaffold-Z"),]
sist.f.z <- sist.f[which(sist.f$V1=="scaffold-Z"),]
sist.m.a.med <- median(sist.m.a$V4,na.rm=T)
sist.f.a.med <- median(sist.f.a$V4,na.rm=T)
sist.m.norm <- sist.m.z$V4/sist.m.a.med
sist.f.norm <- sist.f.z$V4/sist.f.a.med
sist.fm.norm <- log2(sist.f.norm/sist.m.norm)
sist.f.z$fmnorm <- sist.fm.norm

crot.m.a <- crot.m[which(crot.m$V1!="scaffold-Z"),]
crot.f.a <- crot.f[which(crot.f$V1!="scaffold-Z"),]
crot.m.z <- crot.m[which(crot.m$V1=="scaffold-Z"),]
crot.f.z <- crot.f[which(crot.f$V1=="scaffold-Z"),]
crot.m.a.med <- median(crot.m.a$V4,na.rm=T)
crot.f.a.med <- median(crot.f.a$V4,na.rm=T)
crot.m.norm <- crot.m.z$V4/crot.m.a.med
crot.f.norm <- crot.f.z$V4/crot.f.a.med
crot.fm.norm <- log2(crot.f.norm/crot.m.norm)
crot.f.z$fmnorm <- crot.fm.norm

### Plot log2FM per species-------------------------------------------------

par(mfrow=c(4,1))
plot(tham.fm.norm,pch=20,col=alpha('black',0.075),ylim=c(-2,1))
abline(h=-1,lty=2,lwd=1.5,col='black')
abline(h=0,lty=2,lwd=1.5,col='black')
abline(v=9800,col='seagreen')
abline(v=10680,col='seagreen')

plot(dein.fm.norm,pch=20,col=alpha('black',0.075),ylim=c(-2,1))
abline(h=-1,lty=2,lwd=1.5,col='black')
abline(h=0,lty=2,lwd=1.5,col='black')
abline(v=9800,col='seagreen')
abline(v=10680,col='seagreen')

plot(sist.fm.norm,pch=20,col=alpha('black',0.075),ylim=c(-2,1))
abline(h=-1,lty=2,lwd=1.5,col='black')
abline(h=0,lty=2,lwd=1.5,col='black')
abline(v=9800,col='seagreen')
abline(v=10680,col='seagreen')

plot(crot.fm.norm,pch=20,col=alpha('black',0.075),ylim=c(-2,1))
abline(h=-1,lty=2,lwd=1.5,col='black')
abline(h=0,lty=2,lwd=1.5,col='black')
abline(v=9800,col='seagreen')
abline(v=10680,col='seagreen')

### Break out into strata and plot boxplots----------------------------------

tham.f.os <- tham.f.z[which(tham.f.z$V3<98000000),]
tham.f.rs <- tham.f.z[which(tham.f.z$V2>=98000000 & tham.f.z$V3 < 106800000),]
tham.f.par <- tham.f.z[which(tham.f.z$V2>=106800000),]

dein.f.os <- dein.f.z[which(dein.f.z$V3<98000000),]
dein.f.rs <- dein.f.z[which(dein.f.z$V2>=98000000 & dein.f.z$V3 < 106800000),]
dein.f.par <- dein.f.z[which(dein.f.z$V2>=106800000),]

sist.f.os <- sist.f.z[which(sist.f.z$V3<98000000),]
sist.f.rs <- sist.f.z[which(sist.f.z$V2>=98000000 & sist.f.z$V3 < 106800000),]
sist.f.par <- sist.f.z[which(sist.f.z$V2>=106800000),]

crot.f.os <- crot.f.z[which(crot.f.z$V3<98000000),]
crot.f.rs <- crot.f.z[which(crot.f.z$V2>=98000000 & crot.f.z$V3 < 106800000),]
crot.f.par <- crot.f.z[which(crot.f.z$V2>=106800000),]

par(mfrow=c(4,1))
boxplot(tham.f.os$fmnorm,tham.f.rs$fmnorm,tham.f.par$fmnorm,col=c('seagreen','aquamarine','grey'),names=c('OS','RS','PAR'),outline=F,ylim=c(-2,1))
boxplot(dein.f.os$fmnorm,dein.f.rs$fmnorm,dein.f.par$fmnorm,col=c('seagreen','aquamarine','grey'),names=c('OS','RS','PAR'),outline=F,ylim=c(-2,1))
boxplot(sist.f.os$fmnorm,sist.f.rs$fmnorm,sist.f.par$fmnorm,col=c('seagreen','aquamarine','grey'),names=c('OS','RS','PAR'),outline=F,ylim=c(-2,1))
boxplot(crot.f.os$fmnorm,crot.f.rs$fmnorm,crot.f.par$fmnorm,col=c('seagreen','aquamarine','grey'),names=c('OS','RS','PAR'),outline=F,ylim=c(-2,1))
