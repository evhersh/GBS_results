## @knitr load

# packages
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
library(knitr)
library(rmarkdown)


# Data import 

# strata (sample, pop, ms)
mystrata <- read.csv("~/Google Drive/GitHub/Hookeri-GBS/popmap_all.csv")

# vcf to vcfR
vcf <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/final.filtered.snps.vcf")
vcf.dips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.dips.vcf")
vcf.trips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.trips.vcf")
vcf.tets <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.tets.vcf")

# vcfR to genind to genclone
dips.gi <- vcfR2genind(vcf.dips, sep = "/", ploidy=2)
dips.gc <- as.genclone(dips.gi)
sampleorder <- match(indNames(dips.gc), mystrata$id)
strata(dips.gc) <- mystrata[sampleorder,]
setPop(dips.gc) <- ~pop

trips.gi <- vcfR2genind(vcf.trips, sep = "/", ploidy=3)
trips.gc <- as.genclone(trips.gi)
sampleorder <- match(indNames(trips.gc), mystrata$id)
strata(trips.gc) <- mystrata[sampleorder,]
setPop(trips.gc) <- ~pop

tets.gi <- vcfR2genind(vcf.tets, sep = "/", ploidy=4)
tets.gc <- as.genclone(tets.gi)
sampleorder <- match(indNames(tets.gc), mystrata$id)
strata(tets.gc) <- mystrata[sampleorder,]
setPop(tets.gc) <- ~pop

# combine genclones
dipsNtripsNtets.gc <- repool(dips.gc,trips.gc,tets.gc)
#dipsNtripsNtets.gc$pop <- factor(dipsNtripsNtets.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))

AllPops.gc <-as.genclone(dipsNtripsNtets.gc)
AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))

mll(AllPops.gc)


## @knitr MLGs

# calculate raw euclidian distance
dist <- dist(AllPops.gc)

# assign MLG's using raw euclidian distance from dist() [above]
fstats <- filter_stats(AllPops.gc, distance=dist, plot=TRUE)

# looks like this gives the same clone mlg assignments as my IBS stuff
mlg.filter(AllPops.gc, distance=dist) <- 100

mlgtab <- mlg.table(AllPops.gc)





