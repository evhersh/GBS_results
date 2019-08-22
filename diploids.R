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
load("dips.gc.RData")

##### Summary

dips.summary <- summary(dips.gc)

##### HW test and F-statistics

dips.hwt <- hw.test(dips.gc, B=0)
dips.hwt

fstat.dips <- fstat(dips.gc)

Gtest.dips <- gstat.randtest(dips.gc,nsim=99)
Gtest.dips
plot(Gtest.dips)

matFst.dips <- pairwise.fst(dips.gc)

##### inbreeding
seppop.dips <- seppop(dips.gc)
seppop.dips

temp <- inbreeding(seppop.dips, N=100)

Fbar <- sapply(temp, mean)

hist(Fbar, col="firebrick", main="Average inbreeding in sexual individuals")

