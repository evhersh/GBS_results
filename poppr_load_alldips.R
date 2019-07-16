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
library(knitr)
library(ggpubr)
library(hierfstat)

###############
# Data import #
###############

# strata (sample, pop, ms)
mystrata <- read.csv("~/Google Drive/GitHub/Hookeri-GBS/popmap_all.csv")


# vcf to vcfR
vcf.alldips <- read.vcfR("filtered.snps.alldips.vcf")

# vcfR to genind to genclone. Keep return.alleles=TRUE change from 0/1 to A/T 
alldips.gi <- vcfR2genind(vcf.alldips, sep = "/", ploidy=2, return.alleles = TRUE)
#dips.gi <- vcfR2genind(vcf.dips, sep = "/", ploidy=2)
alldips.gc <- as.genclone(alldips.gi)
sampleorder <- match(indNames(alldips.gc), mystrata$id)
strata(alldips.gc) <- mystrata[sampleorder,]
setPop(alldips.gc) <- ~pop

alldips.gc$pop <- factor(alldips.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))

save(alldips.gc, file="alldips.gc.RData")

# looks at multi-locus lineages
mll(AllPops.gc)


####################
# make more colors #
####################

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

cols <- brewer.pal(n = nPop(AllPops.gc), name = "Paired")


# CLONE ID

dist.alldips <- dist(alldips.gc)
#dist.trips <- dist(trips.gc)

# assign MLG's using raw euclidian distance from dist() [above]
fstats.alldips <- filter_stats(alldips.gc, distance=dist.alldips, plot=TRUE)
#fstats.trips <- filter_stats(trips.gc, distance=dist.trips, plot=TRUE)

# looks like this gives the same clone mlg assignments as my IBS stuff
mlg.filter(alldips.gc, distance=dist.alldips) <- 75
mlg.table(alldips.gc)

