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
ddist2 <- prevosti.dist(AllPops.gc)
dist2.mat <- as.matrix(dist(apos.gc))
ddist.apos <- prevosti.dist(apos.gc)
MULTI.ddist.apos <- prevosti.dist(MULTI.apos.gc)

#dist.trips <- dist(trips.gc)

# assign MLG's using raw euclidian distance from dist() [above]
fstats <- filter_stats(AllPops.gc, distance=ddist2, plot=TRUE, stats="THRESHOLD")
sapply(fstats, cutoff_predictor)
fstats2 <- filter_stats(apos.gc, distance=ddist.apos, plot=TRUE)
fstats3 <- filter_stats(MULTI.apos.gc, distance=MULTI.ddist.apos, plot=TRUE)
#fstats.trips <- filter_stats(trips.gc, distance=dist.trips, plot=TRUE)

# looks like this gives the same clone mlg assignments as my IBS stuff
mlg.filter(AllPops.gc, distance=ddist2) <- 0.1
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

sapply(pthresh, cutoff_predictor)

apo.crosspop <- mlg.crosspop(AllPops.gc, sublist=apo.list, df=TRUE)
apo.table <- mlg.table(AllPops.gc, sublist=apo.list)

table.value(apo.table, col.labels=colnames(apo.table))

apo.table.melt <- melt(apo.table)

ggplot(subset(apo.table.melt, value>0), aes(x=Var1, y=Var2))+
  geom_point(aes(size=value), shape=21, colour="black", fill="cornflowerblue")+
  labs(x="population", y="MLG")+
  scale_size_area(max_size=15)

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

# mlg 7
mlg7.gc <- ApoPops.gc[apo.vec==7]
mlg7.dist <- dist(mlg7.gc)
hist(mlg7.dist)
mlg7.mat <- as.matrix(dist(mlg7.gc))
mlg7.mat["S03-A_1",]

# mlg 12
mlg12.gc <- ApoPops.gc[apo.vec==12]
mlg12.dist <- dist(mlg12.gc)
mlg12.mat <- as.matrix(dist(mlg12.gc))
mlg12.mat["SM-A_1",]

# mlg 52
mlg52.gc <- ApoPops.gc[apo.vec==52]
mlg52.dist <- dist(mlg52.gc)
mlg52.mat <- as.matrix(dist(mlg52.gc))
mlg52.mat["C43-A_1",]

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

# mlg 110
mlg110.gc <- ApoPops.gc[apo.vec==110]
mlg110.dist <- dist(mlg110.gc)
mlg110.mat <- as.matrix(dist(mlg110.gc))
mlg110.mat["C23-A_1",]

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

