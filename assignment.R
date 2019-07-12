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

# load data
load("AllPops.gc.RData")

dips.gc <- popsub(AllPops.gc, sublist=sex.list)
apos.gc <- popsub(AllPops.gc, sublist=apo.list)

all.repool<-repool(dips.gc,apos.gc, list=TRUE)
# dapc optimization for sexuals
sextab <- tab(dips.gc, NA.method="mean")
grp <- pop(dips.gc)

sexval <- xvalDapc(sextab, grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

dapc.dips <- dapc(dips.gc, group=dips.gc$pop, n.pca= 10, n.da = 100)
scatter(dapc.dips, grp = dips.gc$pop, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)

pred.apo <- predict.dapc(dapc.dips, newdata = apos.gc)
  
  
  
