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
library(tidyverse)

# load data
load("AllPops.gc.RData")

dips.gc <- popsub(AllPops.gc, sublist=sex.list, drop=FALSE)
apos.gc <- popsub(AllPops.gc, sublist=apo.list, drop=FALSE)

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


col <- rainbow(length(levels(pop(dips.gc))))
col.points <- transp(col[as.integer(pop(dips.gc))],.2)
scatter(dapc.dips, col=col, bg="white", scree.da=0, pch="",
        cstar=0, clab=0, xlim=c(-10,10), legend=TRUE)
par(xpd=TRUE)
points(dapc.dips$ind.coord[,1], dapc.dips$ind.coord[,2], pch=20,
       col=col.points, cex=5)
col.sup <- col[as.integer(pop(apos.gc))]
points(pred.apo$ind.scores[,1], pred.apo$ind.scores[,2], pch=15,
       col=transp(col.sup,.7), cex=2)
add.scatter.eig(dapc.dips$eig,15,1,2, posi="bottomright", inset=.02)
  
table.value(table(pred.apo$assign, pop(apos.gc)), col.lab=levels(pop(apos.gc)))
