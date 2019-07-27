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

load("AllPops.gc.RData")
load("trips.gc.RData")

############
# Clone ID #
############

# calculate raw euclidian distance
dist <- dist(AllPops.gc)
dist2 <- dist(apos.gc)
dist2.mat <- as.matrix(dist(apos.gc))
ddist.apos <- prevosti.dist(apos.gc)

#dist.trips <- dist(trips.gc)

# assign MLG's using raw euclidian distance from dist() [above]
fstats <- filter_stats(AllPops.gc, distance=dist, plot=TRUE)
fstats2 <- filter_stats(apos.gc, distance=ddist.apos, plot=TRUE)
#fstats.trips <- filter_stats(trips.gc, distance=dist.trips, plot=TRUE)

# looks like this gives the same clone mlg assignments as my IBS stuff
mlg.filter(AllPops.gc, distance=ddist) <- 0.1
mlg.table(AllPops.gc)

mlg.filter(apos.gc, distance=ddist.apos) <- 0.13
mlg.table(apos.gc)

apos.gc <- addStrata(apos.gc, apo.vec, name="mlg")

apos.subYK.gc <- popsub(apos.gc, blacklist="SM-A")

# mlg stats
mll(AllPops.gc)

apo.list <- c("C23-A", "C27-A", "C43-A", "C85-A", "C86-A", "C87-A", "C88-A", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A", "L45-A", "L62-A", "S03-A", "SM-A")
sex.list <- c("B42-S", "B46-S", "B49-S", "B53-S", "B60-S", "C59-S", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L45-S", "L62-S")

ApoPops.gc <- popsub(AllPops.gc, sublist=apo.list, drop=FALSE)
mlg.filter(ApoPops.gc, distance=ddist.apos) <- 0.13


apo.crosspop <- mlg.crosspop(AllPops.gc, sublist=apo.list)
mlg.table(AllPops.gc, sublist=apo.list)

apo.vec <- mlg.vector(ApoPops.gc)
mlgIDs <- mlg.id(ApoPops.gc)
mlg57inds <- mlgIDs[as.character(57)]

# mlg 57
mlg57.gc <- ApoPops.gc[apo.vec==57]
mlg57.dist <- provesti.dist(mlg57.gc)
hist(mlg57.dist)
mlg57.mat <- as.matrix(dist(mlg57.gc))
mlg57.mat["C27-A_1",]

# mlg 103
mlg103.gc <- ApoPops.gc[apo.vec==103]
mlg103.dist <- dist(mlg103.gc)
hist(mlg103.dist)
mlg103.mat <- as.matrix(dist(mlg103.gc))
mlg103.mat["L06-A_1",]

# mlg 4
mlg4.gc <- ApoPops.gc[apo.vec==4]
mlg4.dist <- dist(mlg4.gc)
hist(mlg4.dist)
mlg4.mat <- as.matrix(dist(mlg4.gc))
mlg4.mat["S03-A_1",]

# mlg 12
mlg12.gc <- ApoPops.gc[apo.vec==12]
mlg12.dist <- dist(mlg12.gc)
mlg12.mat <- as.matrix(dist(mlg12.gc))
mlg12.mat["SM-A_1",]

# mlg 55
mlg55.gc <- ApoPops.gc[apo.vec==55]
mlg55.dist <- dist(mlg55.gc)
mlg55.mat <- as.matrix(dist(mlg55.gc))
mlg55.mat["C43-A_1",]

# mlg 91
mlg91.gc <- ApoPops.gc[apo.vec==91]
mlg91.dist <- dist(mlg91.gc)
mlg91.mat <- as.matrix(dist(mlg91.gc))
mlg91.mat["L16-A_1",]

# mlg 104
mlg104.gc <- ApoPops.gc[apo.vec==104]
mlg104.dist <- dist(mlg104.gc)
mlg104.mat <- as.matrix(dist(mlg104.gc))
mlg104.mat["L62-A_5",]

# mlg 108
mlg108.gc <- ApoPops.gc[apo.vec==108]
mlg108.dist <- dist(mlg108.gc)
mlg108.mat <- as.matrix(dist(mlg108.gc))
mlg108.mat["C23-A_1",]

# mlg 114
mlg114.gc <- ApoPops.gc[apo.vec==114]
mlg114.dist <- dist(mlg114.gc)
mlg114.mat <- as.matrix(dist(mlg114.gc))
mlg114.mat["L39-A_1",]

# mlg 105
mlg105.gc <- ApoPops.gc[apo.vec==105]
mlg105.dist <- dist(mlg105.gc)
mlg105.mat <- as.matrix(dist(mlg105.gc))
mlg105.mat["L39-A_4",]

all.mlg.dists <- c(mlg57.dist,mlg103.dist,mlg4.dist,mlg12.dist,mlg55.dist,mlg91.dist,mlg104.dist,mlg108.dist,mlg114.dist)

hist(all.mlg.dists)

dist2.mat["L62-A_5","L39-A_4"]
dist2.mat["C43-A_1","S03-A_1"]

##### genclone object with pop defined as MLG

addStrata(mlg57.gc, name=mlg) <- "mlg57"
