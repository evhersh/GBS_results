##### packages
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
library(ggpubr)
library(pegas)
library(hierfstat)

##### load
load("diploid_stats.RData")

##### Summary

#dips.summary <- summary(dips.gc)
dips.summary

##### HW test and F-statistics

dips.gp <- genind2genpop(dips.gc)

#dips.hwt <- hw.test(dips.gc)
dips.hwt

#fstat.dips <- fstat(dips.gc)
fstat.dips
# low inbreeding (ind/pop) = Fis

fstat.c59 <- fstat(seppop.dips$`C59-S`)

#Gtest.dips <- gstat.randtest(dips.gc,nsim=99)
Gtest.dips
plot(Gtest.dips)

matFst.dips <- pairwise.fst(dips.gc)
matFst.dips
colnames(matFst.dips) <- popNames(dips.gc)
rownames(matFst.dips) <- popNames(dips.gc)

matFst.dips2 <- as.matrix(matFst.dips)
colnames(matFst.dips2) <- popNames(dips.gc)
rownames(matFst.dips2) <- popNames(dips.gc)


dips.tree.fst <- nj(matFst.dips)
plot(dips.tree.fst, type="unr", tip.col="coral3", font=2)
annot <- round(dips.tree.fst$edge.length,2)
edgelabels(annot[annot>0], which(annot>0), frame="n")
add.scale.bar()

##### inbreeding
seppop.dips <- seppop(dips.gc)
seppop.dips

temp <- inbreeding(seppop.dips, N=100)

Fbar <- sapply(temp, mean)

hist(Fbar, col="firebrick", main="Average inbreeding in sexual individuals")

dips.stats <- basic.stats(dips.gc)

fdens <- inbreeding(dips.gc, res.type = "function")

Fest <- inbreeding(dips.gc, res.type = "estimate")
mostInbred <- which.max(Fest)
