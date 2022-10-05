############################################################################
# Divergence between Crotalus and Anolis autosomal orthologs
############################################################################

### Goal: examine distribution of dS between Crotalus and Anolis orthologs
### from autosomes to determine a lineage-specific mutation rate estimate
### to apply to divergence between Z and W gametologs.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Load dependencies-------------------------------------------------------

library(Rmisc)
library(scales)

### Read in data------------------------------------------------------------

# dN/dS results
ca.div <- read.table('./crotalus_anolis_ortholog.autosome.dnds.txt',header=T)

### Filter data-------------------------------------------------------------

ca.div$ds[ca.div$ds<0] <- NA
ca.div <- ca.div[which(ca.div$ds<2),]

ca.4fds[ca.4fds$Nei<0] <- NA
ca.4fds <- merge(ca.4fds,ca.div,by='ortholog')

boxplot(ca.div$ds,ca.4fds$Nei)

### Calculate median dS-----------------------------------------------------

ds.med <- median(ca.div$ds,na.rm=T)

### Determine lineage-specific mutation rate--------------------------------

# We can divide the median divergence estimate by a 'known' divergence time
# between Crotalus and Anolis to get a representative mutation rate.

# The median divergence time for these species on timetree.org is 161 MYA.

# The mutation rate calculation formulat is (dS / divergence time in years) / 2

time <- 167000000
rate <- (ds.med/time)/2

rate <- 2.4e-09
# To get a sex-linked divergence rate, we need to account for 1) male-biased
# mutation and 2) female-specific mutation. The first is known (uZ/uA) = 1.1
# (Schield et al. 2021). The female rate is unknown, but we can conservatively
# assume it is one-fourth the autosomal rate.

rate.z <- rate*1.1
rate.w <- rate/4    # This was divided by 3, but 4 makes more sense I think based on alphaM estimates from Z and Autosomes

# Finally, the combined sex-linked mutation rate is:

rate.sex <- rate.z+rate.w

### Transform ZW gametolog dS estimates to time-----------------------------

# Read in data
div <- read.table('./ZW_gametolog_pairwise_Nei_d_dn_ds.txt',header=T)

# Split gametologs into inferred strata

s2 <- div[which(div$start<70000000),]
s1 <- div[which(div$start>70000000 & div$start<98000000),]
s3 <- div[which(div$start>98000000),]

# Calculate median divergence time per stratum

mean(s1$dS)/rate.sex
mean(s2$dS)/rate.sex
mean(s3$dS)/rate.sex

wilcox.test(s1$dS,s2$dS)
wilcox.test(s1$dS,s3$dS)
wilcox.test(s2$dS,s3$dS)
