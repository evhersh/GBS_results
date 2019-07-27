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

# unspit apo replicate distances
mean(rep.dists) # 58.15389
hist(rep.dists)

# within-pop apo distances
load("dips.gc.RData")
load("trips.gc.RData") # trips.gc
distmat.trips <- as.matrix(dist(trips.gc))
labels(distmat.trips)
distmat.trips["C43-A_3","S03-A_1"]
############################
# within-MLG apo distances #
############################
apo.vec <- mlg.vector(ApoPops.gc)
mlgIDs <- mlg.id(ApoPops.gc)
mlg57inds <- mlgIDs[as.character(57)]

# mlg 57
mlg57.gc <- ApoPops.gc[apo.vec==57]
mlg57.dist <- dist(mlg57.gc)
hist(mlg57.dist)

# mlg 103
mlg103.gc <- ApoPops.gc[apo.vec==103]
mlg103.dist <- dist(mlg103.gc)
hist(mlg103.dist)

# mlg 4
mlg4.gc <- ApoPops.gc[apo.vec==4]
mlg4.dist <- dist(mlg4.gc)
hist(mlg4.dist)

# mlg 12
mlg12.gc <- ApoPops.gc[apo.vec==12]
mlg12.dist <- dist(mlg12.gc)

# mlg 55
mlg55.gc <- ApoPops.gc[apo.vec==55]
mlg55.dist <- dist(mlg55.gc)

# mlg 91
mlg91.gc <- ApoPops.gc[apo.vec==91]
mlg91.dist <- dist(mlg91.gc)

# mlg 104
mlg104.gc <- ApoPops.gc[apo.vec==104]
mlg104.dist <- dist(mlg104.gc)

# mlg 108
mlg108.gc <- ApoPops.gc[apo.vec==108]
mlg108.dist <- dist(mlg108.gc)

# mlg 114
mlg114.gc <- ApoPops.gc[apo.vec==114]
mlg114.dist <- dist(mlg114.gc)

all.mlg.dists <- c(mlg57.dist,mlg103.dist,mlg4.dist,mlg12.dist,mlg55.dist,mlg91.dist,mlg104.dist,mlg108.dist,mlg114.dist)

hist(all.mlg.dists)
apo.mlg.dists <-as.data.frame(all.mlg.dists)
apo.mlg.dists$ms <- rep("apo", nrow(apo.mlg.dists))
apo.mlg.dists

colnames(apo.mlg.dists) <- c("dist","ms")
############
# sex pops #
############

dips.sep.gc <- seppop(dips.gc)

dist.b42 <- dist(dips.sep.gc$`B42-S`)
dist.b46 <- dist(dips.sep.gc$`B46-S`)
dist.b49 <- dist(dips.sep.gc$`B49-S`)
dist.b53 <- dist(dips.sep.gc$`B53-S`)
dist.b60 <- dist(dips.sep.gc$`B60-S`)
dist.C59 <- dist(dips.sep.gc$`C59-S`)
dist.l05 <- dist(dips.sep.gc$`L05-S`)
dist.l08 <- dist(dips.sep.gc$`L08-S`)
dist.l10 <- dist(dips.sep.gc$`L10-S`)
dist.l11 <- dist(dips.sep.gc$`L11-S`)
dist.l12 <- dist(dips.sep.gc$`L12-S`)
dist.l13 <- dist(dips.sep.gc$`L13-S`)
dist.l45 <- dist(dips.sep.gc$`L45-S`)
dist.l62 <- dist(dips.sep.gc$`L62-S`)

sex.mlg.dists<-c(dist.b42,dist.b46,dist.b49,dist.b53,dist.b60,dist.C59,dist.l05,dist.l08,dist.l10,dist.l11,dist.l12,dist.l13,dist.l45,dist.l62)
hist(sex.mlg.dists)

sex.mlg.dists <-as.data.frame(sex.mlg.dists)
sex.mlg.dists$ms <- rep("sex", nrow(sex.mlg.dists))
sex.mlg.dists

colnames(sex.mlg.dists) <- c("dist","ms")

wInPopDists <- rbind(apo.mlg.dists,sex.mlg.dists)

###### ggplot ######

gg.wInPopDists <- ggplot(data=wInPopDists, aes(x=dist, fill=ms))+
  geom_histogram(aes(y=..density..), bins=50, alpha=0.6, color="black", position="identity", size=.75)+
  theme_classic()+
  scale_fill_manual(values=c("red", "blue"), name="Mating System", labels=c("Apomictic", "Sexual"))+
  geom_density(alpha=.3)+
  xlim(0,160)+
  geom_vline(aes(xintercept=54.67), color="black", linetype="dashed")+
  labs(x="Pairwise euclidian distance", y="Density")

gg.wInPopDists <- ggplot(data=wInPopDists, aes(x=dist, fill=ms))+
  geom_histogram(aes(y=..density..), bins=50, alpha=0.6, color="black", position="identity", size=.75)+
  theme_classic()+
  scale_fill_manual(values=c("red", "blue"), name="Mating System", labels=c("Apomictic", "Sexual"))+
  geom_density(alpha=.3)+
  xlim(0,160)+
  geom_vline(aes(xintercept=54.67), color="black", linetype="dashed")+
  labs(x="Pairwise euclidian distance", y="Density")

#scale_x_continuous(breaks=seq(0.70, 1, 0.025))+



#### different distances ####

hist(dist)

ddist <- prevosti.dist(AllPops.gc)
ddist <- prevosti.dist(AllPops.gc)
ddist.dips <- prevosti.dist(dips.gc)

hist(ddist)

filter_stats(AllPops.gc, distance=ddist, plot=TRUE)
filter_stats(AllPops.gc, distance=dist, plot=TRUE)
filter_stats(dips.gc, distance=ddist.dips, plot=TRUE)


# 0.09 - 0.13 = 49 MLGs
mlg.filter(AllPops.gc, distance=ddist) <- 0.1
mlg.table(AllPops.gc)

mlg.filter(AllPops.gc, distance=ddist) <- 
mlg.table(AllPops.gc)

ddist.apos.df <- as.data.frame(ddist.apos)

ggplot(data=ddist.apos.df)+
  geom_histogram()
