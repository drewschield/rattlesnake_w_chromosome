############################################################################
# Testing the 'refugium hypothesis' on autosomes, Z, and W chromosomes
############################################################################

### Goal: test the hypothesis that sex chromosomes, and the W in particular,
### have an excess of repeat elements (and ERVs/mdg4, specifically). 

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set RI values (calculated in Supp. Tables-------------------------------

ri.tot.a <- -0.044
ri.tot.z <- 0.21 
ri.tot.w <- 1.8
  
ri.erv.a <- -0.42
ri.erv.z <- 0.7
ri.erv.w <- 16.64
  
ri.mdg.a <- -0.24
ri.mdg.z <- 0.87
ri.mdg.w <- 7.5

ri.l1.a <- -0.16
ri.l1.z <- 0.7
ri.l1.w <- 4.5

ri.cr1.a <- -0.006
ri.cr1.z <- 0.059
ri.cr1.w <- 0.028

### Plot RI values for repeat categories------------------------------------

par(mfrow=c(1,5))
barplot(c(ri.tot.a,ri.tot.z,ri.tot.w),col=c('grey','seagreen','navy'),ylab='Repeat refugium index (RI)',names=c('Auto','Z','W'))
barplot(c(ri.erv.a,ri.erv.z,ri.erv.w),col=c('grey','seagreen','navy'),ylab='ERV refugium index (RI)',names=c('Auto','Z','W'))
barplot(c(ri.mdg.a,ri.mdg.z,ri.mdg.w),col=c('grey','seagreen','navy'),ylab='mdg4 refugium index (RI)',names=c('Auto','Z','W'))
barplot(c(ri.l1.a,ri.l1.z,ri.l1.w),col=c('grey','seagreen','navy'),ylab='L1 refugium index (RI)',names=c('Auto','Z','W'))
barplot(c(ri.cr1.a,ri.cr1.z,ri.cr1.w),col=c('grey','seagreen','navy'),ylab='CR1 refugium index (RI)',names=c('Auto','Z','W'))

### Pie charts of sequence bp and repeat bp---------------------------------

## sequence bp proportions
prop.a <- 0.8940894391
prop.z <- 0.08691597094
prop.w <- 0.01899458998

par(mfrow=c(1,6))
slices.bp <- c(prop.a,prop.z,prop.w)
pie(slices.bp,col=c('grey','seagreen','navy'),labels=c('Auto','Z','W'),main='Sequence bp')

slices.rep <- c(437748675,54041723,20210532)
pie(slices.rep,col=c('grey','seagreen','navy'),labels=c('Auto','Z','W'),main='Repeat bp')

slices.erv <- c(825697,235896,535169)
pie(slices.erv,col=c('grey','seagreen','navy'),labels=c('Auto','Z','W'),main='ERV bp')

slices.mdg <- c(9846145,2379345,2355232)
pie(slices.mdg,col=c('grey','seagreen','navy'),labels=c('Auto','Z','W'),main='mgd4 bp')

slices.l1 <- c(10015653,1997066,1459196)
pie(slices.l1,col=c('grey','seagreen','navy'),labels=c('Auto','Z','W'),main='L1 bp')

slices.cr1 <- c(36542683,3786887,803511)
pie(slices.cr1,col=c('grey','seagreen','navy'),labels=c('Auto','Z','W'),main='CR1 bp')

### chi-square tests of uniform distribution of repeats---------------------

## Here, we're testing if the observed bp of repeats in whichever category
## on autosomes, Z, and W chromosomes matches a uniform distribution, given
## the lengths of the chromosomes.

## The observed and expected values have been calculated elsewhere.

## set expected and observed values
tot.obs.a <- 437748675
tot.obs.z <- 54041723
tot.obs.w <- 20210532

erv.obs.a <- 825697
erv.obs.z <- 235896
erv.obs.w <- 535169

mdg.obs.a <- 9846145
mdg.obs.z <- 2379345
mdg.obs.w <- 2355232

l1.obs.a <- 10015653
l1.obs.z <- 1997066
l1.obs.w <- 1459196

cr1.obs.a <- 36542683
cr1.obs.z <- 3786887
cr1.obs.w <- 803511

flLTR.obs.a <- 13603 
flLTR.obs.z <- 1786
flLTR.obs.w <- 572

prop.a <- 0.8940894391
prop.z <- 0.08691597094
prop.w <- 0.01899458998

## run chi-square tests

# all repeats
tot.obsfreq <- c(tot.obs.a,tot.obs.z,tot.obs.w)
tot.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(tot.obsfreq,p=tot.nullprobs)

# ERVs
erv.obsfreq <- c(erv.obs.a,erv.obs.z,erv.obs.w)
erv.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(erv.obsfreq,p=erv.nullprobs)

# mdg4 elements
mdg.obsfreq <- c(mdg.obs.a,mdg.obs.z,mdg.obs.w)
mdg.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(mdg.obsfreq,p=mdg.nullprobs)

# L1 elements
l1.obsfreq <- c(l1.obs.a,l1.obs.z,l1.obs.w)
l1.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(l1.obsfreq,p=l1.nullprobs)

# CR1 elements
cr1.obsfreq <- c(cr1.obs.a,cr1.obs.z,cr1.obs.w)
cr1.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(cr1.obsfreq,p=cr1.nullprobs)

# fl-LTRs
flLTR.obsfreq <- c(flLTR.obs.a,flLTR.obs.z,flLTR.obs.w)
flLTR.nullprobs <- c(prop.a,prop.z,prop.w)

chisq.test(flLTR.obsfreq,p=flLTR.nullprobs)


