library(Rmisc)

auto <- read.table('./coverage_exp_crotalus/female.all.CroVir.autosome.gene.regions.bed',header=T)
chrw <- read.table('./coverage_exp_crotalus/female.all.chrW_agouti.gene.regions.bed',header=T)

auto.CV0011 <- median(auto$CV0011)*0.5
auto.CV0629 <- median(auto$CV0629)*0.5
auto.CV0646 <- median(auto$CV0646)*0.5
auto.CV0650 <- median(auto$CV0650)*0.5

chrw.CV0011 <- chrw$CV0011/auto.CV0011
chrw.CV0629 <- chrw$CV0629/auto.CV0629
chrw.CV0646 <- chrw$CV0646/auto.CV0646
chrw.CV0650 <- chrw$CV0650/auto.CV0650

chrw.mean <- (chrw.CV0011+chrw.CV0629+chrw.CV0646+chrw.CV0650)/4

chrw$norm.cov <- chrw.mean

# Get quick stats on autosomal normalized distribution
auto.CV0011.norm <- (auto$CV0011*0.5)/auto.CV0011
auto.CV0629.norm <- (auto$CV0629*0.5)/auto.CV0629
auto.CV0646.norm <- (auto$CV0646*0.5)/auto.CV0646
auto.CV0650.norm <- (auto$CV0650*0.5)/auto.CV0650
auto.norm.mean <- (auto.CV0011.norm+auto.CV0629.norm+auto.CV0646.norm+auto.CV0650.norm)/4

mean(auto.norm.mean)
sd(auto.norm.mean)
par(mfrow=c(1,1))
plot(density(auto.norm.mean),xlim=c(0,10))
count(auto.norm.mean>1.75)
208/(208+15695) # Less than 2.5% of autosomal genes have coverage greater than 1.75 - seems like a reasonable benchmark for assigning two or more copies.

write.table(chrw,file='./coverage_exp_crotalus/female.all.chrW_agouti.gene.regions.norm.txt',quote=F,row.names=F,sep = "\t")

# Basic statistics

count(chrw$norm.cov > 1.75)
58/228
count(chrw$norm.cov > 10)

# Read in copy number inferences with Z coordinates (pre-formatted)

coord <- read.table('./resources/W_copy-number.txt',header=T)
coord.noZ <- coord[which(coord$z.chrom.position==0),]

# output at 4 x 6
par(mfrow=c(1,1))
plot(coord$z.chrom.position,coord$copy.number,pch=20,ylim=c(0,10),ylab='Copy number',xlab='Z chromosome position')
abline(h=1,lty=2)
plot(1,coord.noZ$copy.number)

plot(coord$z.chrom.position,coord$copy.number,pch=20,ylim=c(0,30))
plot(coord$z.chrom.position,coord$copy.number,pch=20)


max(coord$Z_position,na.rm=T)
cor.test(coord$Z_position,coord$norm.cov,method='spearman')

plot(chrw.mean,pch=20)
abline(h=0.5,lty=2)

boxplot(chrw.mean,pch=20,outline=F)

library(Rmisc)
mean(chrw.mean)
sd(chrw.mean)
count(chrw.mean>0.75)

plot(hist(chrw.mean,breaks=10))
