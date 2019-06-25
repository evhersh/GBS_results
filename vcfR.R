############
# Packages
############
library("vcfR")
library("poppr")
library("ape")
library("RColorBrewer")
library(gdsfmt)
library(SNPRelate)
library("ggplot2")
library(reshape2)
library(plotly)
library(adegenet)

######
# vcfR
######

#first assembly
# vcf <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/SNP.DP3g95p5maf05.HWE.recode.vcf")
# read data, one ploidy at a time, and then combine

# strata (sample, pop, ms)
mystrata <- read.csv("~/Google Drive/GitHub/Hookeri-GBS/popmap_all.csv")

# vcf to vcfR
vcf <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/final.filtered.snps.vcf", verbose = TRUE)
# vcf.dips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.dips.vcf")
# vcf.trips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.trips.vcf")
# vcf.tets <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.tets.vcf")

# mess with CNV
gt<-extract.gt(vcf)
hets<-is_het(gt)
is.na(vcf@gt[,-1][!hets]) <- TRUE
ad <- extract.gt(vcf, element = "AD")
ad1 <- masplit(ad, record = 1)
ad2 <- masplit(ad, record = 2)
freq1 <- ad1/(ad1+ad2)
freq2 <- ad2/(ad1+ad2)
myPeaks1 <- freq_peak(freq1, getPOS(vcf))
is.na(myPeaks1$peaks[myPeaks1$counts < 20]) <- TRUE
myPeaks2 <- freq_peak(freq2, getPOS(vcf), lhs = FALSE)
is.na(myPeaks2$peaks[myPeaks2$counts < 20]) <- TRUE
myPeaks1
peak_to_ploid(myPeaks1)
freq_peak_plot(pos = getPOS(vcf), ab1 = freq1, ab2 = freq2, fp1 = myPeaks1, fp2=myPeaks2)

# vcfR Knaus
# determine ploidy
knitr::kable(vcf@gt[c(1:2,11,30),1:4])
ad <- extract.gt(vcf, element = 'AD')
knitr::kable(ad[c(1:2,11,30),1:4])
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)

ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)



# diploid !!!!
hist(ad2[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# diploid
hist(ad2[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n", main="Diploid", xlab="")
hist(ad1[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))

tiff("ABplotDip.tiff", res=300, units="in", width=8.6, height=5.8)
hist(ad2[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n", main="Diploid", xlab="")
hist(ad1[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
dev.off()

# triploid !!!!
hist(ad2[,"L17-A_1"], breaks = seq(0,1,by=0.02), col = "indianred", xaxt="n", main="Triploid")
hist(ad1[,"L17-A_1"], breaks = seq(0,1,by=0.02), col = "tan2", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

tiff("ABplotTrip.tiff", res=300, units="in", width=8.6, height=5.8)
hist(ad2[,"L17-A_1"], breaks = seq(0,1,by=0.02), col = "indianred", xaxt="n", main="Triploid")
hist(ad1[,"L17-A_1"], breaks = seq(0,1,by=0.02), col = "tan2", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))
dev.off()

# Dan's
hist(ad2[,"L62-A_5"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"L62-A_5"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# L39
hist(ad2[,"L39-A_4"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"L39-A_4"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# tetraploid???
tiff("ABplotTet.tiff", res=300, units="in", width=8.6, height=5.8)
hist(ad2[,"S03-A_4"], breaks = seq(0,1,by=0.02), col = "indianred", xaxt="n", main="Tetraploid", xlab="Allele balance")
hist(ad1[,"S03-A_4"], breaks = seq(0,1,by=0.02), col = "tan2", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))
dev.off()

# L16-A_2
hist(ad2[,"L16-A_5"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"L16-A_5"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# L45-A_1 (diploid in a triploid maze...)
hist(ad2[,"L45-A_1"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"L45-A_1"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# C59 (ind 3 looks diploid, the others look weird...maybe remove from analysis?)
hist(ad2[,"C59-S_2"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"C59-S_2"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

#### Clean-up? ####
# ad <- extract.gt(vcf, element = 'AD')
# #ad[1:3,1:4]
# 
# allele1 <- masplit(ad, record = 1)
# allele2 <- masplit(ad, record = 2)
# 
# # Subset to a vector for manipulation.
# tmp <- allele1[,"B42-S_1"]
# sum(tmp == 0, na.rm = TRUE)
# tmp <- tmp[tmp > 0]
# tmp <- tmp[tmp <= 100]
# 
# hist(tmp, breaks=seq(0,100,by=1), col="#808080", main = "B42-S_1")
# 
# sums <- apply(allele1, MARGIN=2, quantile, probs=c(0.15, 0.95), na.rm=TRUE)
# sums[,"B42-S_1"]
# 
# abline(v=sums[,"B42-S_1"], col=2, lwd=2)
# 
# tmp <- allele2[,"B42-S_1"]
# tmp <- tmp[tmp>0]
# tmp <- tmp[tmp<=100]
# 
# hist(tmp, breaks=seq(1,100,by=1), col="#808080", main="B42-S_1")
# 
# sums[,"B42-S_1"]
# 
# abline(v=sums[,"B42-S_1"], col=2, lwd=2)
# 
# ad <- extract.gt(vcf, element = 'AD')
# allele1 <- masplit(ad, record = 1)
# allele2 <- masplit(ad, record = 2)
# boxplot(allele1, las=3)
# 
# sums <- apply(allele1, MARGIN=2, quantile, probs=c(0.15, 0.95), na.rm=TRUE)
# # Allele 1
# dp2 <- sweep(allele1, MARGIN=2, FUN = "-", sums[1,])
# #allele1[dp2 < 0] <- NA
# vcf@gt[,-1][ dp2 < 0 & !is.na(vcf@gt[,-1]) ] <- NA
# dp2 <- sweep(allele1, MARGIN=2, FUN = "-", sums[2,])
# #allele1[dp2 > 0] <- NA
# vcf@gt[,-1][dp2 > 0] <- NA
# # Allele 2
# dp2 <- sweep(allele2, MARGIN=2, FUN = "-", sums[1,])
# vcf@gt[,-1][ dp2 < 0 & !is.na(vcf@gt[,-1]) ] <- NA
# dp2 <- sweep(allele2, MARGIN=2, FUN = "-", sums[2,])
# vcf@gt[,-1][dp2 > 0] <- NA
# 
# ad <- extract.gt(vcf, element = 'AD')
# allele1 <- masplit(ad, record = 1)
# allele2 <- masplit(ad, record = 2)
# boxplot(allele1, las=3)
# 
# gt <- extract.gt(vcf, element = 'GT')
# hets <- is_het(gt)
# is.na( ad[ !hets ] ) <- TRUE
# 
# allele1 <- masplit(ad, record = 1)
# allele2 <- masplit(ad, record = 2)
# 
# ad1 <- allele1 / (allele1 + allele2)
# ad2 <- allele2 / (allele1 + allele2)





#### random ####
# diploid assembly (no HWE)
vcf <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/random.vcf")
is.biallelic(vcf)
my_genlight <- vcfR2genlight(vcf, n.cores=2)
my_genind <- vcfR2genind(vcf, ploidy=2)
plot(my_genlight)
vcf
queryMETA(vcf)
queryMETA(vcf, element = 'DP')
queryMETA(vcf, element = 'FORMAT=<ID=DP')

myMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- myMiss/nrow(vcf)
library(RColorBrewer)
palette(brewer.pal(n=12, name = 'Set3'))

par(mar = c(12,4,4,2))
barplot(myMiss, las = 2, col = 1:12)
title(ylab = "Missingness (%)")

myMiss <- apply(dp, MARGIN = 1, function(x){ sum(is.na(x)) })
myMiss <- myMiss/ncol(vcf@gt[,-1])

hist(myMiss, col = "#8DD3C7", xlab = "Missingness (%)")

dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
par(mar=c(12,4,4,2))
boxplot(dp, col=2:8, las=3)
title(ylab = "Depth (DP)")