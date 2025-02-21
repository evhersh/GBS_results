##### Packages #####

library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(reshape2)
library(plotly)
library(adegenet)
library(knitr)
library(ggpubr)
library(hierfstat)
#library(ggtree)
library(ggnewscale)
library(dplyr)
library(here)

##### Data import #####

# strata (sample, pop, ms)
mystrata <- read.csv(here("data", "popmap_all.csv"))

# vcf to vcfR
vcf <- read.vcfR(here("data", "final.filtered.snps.vcf"))
vcf.dips <- read.vcfR(here("data", "filtered.dips.vcf"))
vcf.trips <- read.vcfR(here("data", "filtered.trips.vcf"))
vcf.tets <- read.vcfR(here("data", "filtered.tets.vcf"))

# vcfR to genind to genclone
dips.gi <- vcfR2genind(vcf.dips, sep = "/", ploidy=2, return.alleles = TRUE)
dips.gc <- as.genclone(dips.gi)
sampleorder <- match(indNames(dips.gc), mystrata$id)
strata(dips.gc) <- mystrata[sampleorder,]
setPop(dips.gc) <- ~pop

trips.gi <- vcfR2genind(vcf.trips, sep = "/", ploidy=3, return.alleles = TRUE)
trips.gl <- vcfR2genlight(vcf.trips)
trips.gc <- as.genclone(trips.gi)
sampleorder <- match(indNames(trips.gc), mystrata$id)
strata(trips.gc) <- mystrata[sampleorder,]
setPop(trips.gc) <- ~pop

tets.gi <- vcfR2genind(vcf.tets, sep = "/", ploidy=4, return.alleles = TRUE)
tets.gc <- as.genclone(tets.gi)
sampleorder <- match(indNames(tets.gc), mystrata$id)
strata(tets.gc) <- mystrata[sampleorder,]
setPop(tets.gc) <- ~pop

# combine genclones
dipsNtripsNtets.gc <- repool(dips.gc,trips.gc,tets.gc)

AllPops.gc <-as.genclone(dipsNtripsNtets.gc)
# AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("B53-S", "B60-S", "B42-S", 
#                                                   "B46-S", "B49-S", "L62-S", 
#                                                   "L62-A", "L05-S", "L08-S", 
#                                                   "L10-S", "L11-S", "L12-S", 
#                                                   "L13-S", "L06-A", "L16-A", 
#                                                   "L17-A", "L39-A", "L41-A",
#                                                   "L45-S", "L45-A", "C87-A", 
#                                                   "C86-A", "C88-A", "C85-A", 
#                                                   "C27-A", "C23-A", "C43-A", 
#                                                   "S03-A", "SM-A", "C59-S"))

AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("CO1-S", "CO2-S", "CO3-S", 
                                                  "CO4-S", "CO5-S", "WY1-S", 
                                                  "WY1-A", "WY2-S", "WY3-S", 
                                                  "WY4-S", "WY6-S", "WY7-S", 
                                                  "WY8-S", "WY9-A", "WY10-A", 
                                                  "WY11-A", "WY12-A", "WY13-A",
                                                  "WY14-S", "WY14-A", "WY15-A", 
                                                  "MT1-A", "MT2-A", "MT3-A", 
                                                  "MT4-A", "ND-A", "BC-A", 
                                                  "SK-A", "YK1-A", "YK2-S"))

apo.list <- c("WY1-A", "WY9-A", "WY10-A", 
              "WY11-A", "WY12-A", "WY13-A",
              "WY14-A", "WY15-A", 
              "MT1-A", "MT2-A", "MT3-A", 
              "MT4-A", "ND-A", "BC-A", 
              "SK-A", "YK1-A")
sex.list <- c("CO1-S", "CO2-S", "CO3-S", 
              "CO4-S", "CO5-S", "WY1-S", 
              "WY2-S", "WY3-S", 
              "WY4-S", "WY6-S", "WY7-S", 
              "WY8-S",  
              "WY14-S", "YK2-S")

dips.gc <- popsub(AllPops.gc, sublist=sex.list, drop=FALSE)
apos.gc <- popsub(AllPops.gc, sublist=apo.list, drop=FALSE)

saveRDS(AllPops.gc, here("data", "AllPops_gc.rds"))

############
# Distance #
############

ddist2 <- prevosti.dist(AllPops.gc)
ddist.apos <- prevosti.dist(apos.gc)
ddist.dips <- prevosti.dist(dips.gc)

mlg.filter(AllPops.gc, distance=ddist2) <- 0.1
mlg.table(AllPops.gc)

mlg.filter(apos.gc, distance=ddist.apos) <- 0.1
mlg.table(apos.gc)

apo.vec <- mlg.vector(apos.gc)
mlgIDs <- mlg.id(apos.gc)
apos.gc <- addStrata(apos.gc, apo.vec, name="mlg")

# plot distances

ddist.apos.df <- melt(as.matrix(ddist.apos), varnames=c("row", "col"))
apodist.df <- ddist.apos.df[as.numeric(ddist.apos.df$row) > as.numeric(ddist.apos.df$col),]

apodist.df %>% filter(between(value, 0.055, 0.1))
# apodist.df %>% filter(between(value, 0.07, 0.1))
apodist.df %>% filter(between(value, 0.1, 0.15))
# apodist.df %>% filter(between(value, 0.25, 0.4))

gg.MLGdists <- ggplot(data=apodist.df, aes(x=value))+
  geom_histogram(aes(y=..density..), bins=100, color="black", fill="grey50", size=.75)+
  theme_classic()+
  geom_density(fill="darkgrey",alpha=.3)+
  #geom_jitter(data=error.df, aes(x=rep.dists, y=20), pch=21, fill="slateblue", alpha=0.5)+
  #geom_point(data=error.mean, aes(x=mean, y=21), pch=21, fill="white")+
  #geom_errorbarh(data=error.mean, aes(xmin=lower.ci, xmax=upper.ci, y=21), height=0, inherit.aes = FALSE)+
  geom_vline(aes(xintercept=0.02647497), color="red", linetype="twodash", size=.75)+
  geom_vline(aes(xintercept=0.1), color="black", linetype="dashed", size=.75)+
  #geom_vline(aes(xintercept=0.055), color="black", size=.75)+
  labs(x="", y="Density")+
  theme(legend.position = "none")+
  xlim(0, 0.35)

#sexual pairwise dist prevosti

ddist.dips.df <- melt(as.matrix(ddist.dips), varnames=c("row", "col"))
sexdist.df <- ddist.dips.df[as.numeric(ddist.dips.df$row) > as.numeric(ddist.dips.df$col),]

gg.sexdists <- ggplot(data=sexdist.df, aes(x=value))+
  geom_histogram(aes(y=..density..), bins=100, color="black",fill="grey80", size=.75)+
  theme_classic()+
  geom_density(fill="white", alpha=.2)+
  labs(x="Pairwise prevosti distance", y="Density")+
  theme(legend.position = "none")+
  xlim(0, 0.35)

png("distance_apos.png", height=5, width=8, res=300, units="in")
gg.MLGdists
dev.off()

png("./figures/distance.png", height=7, width=8, res=300, units="in")
ggarrange(gg.MLGdists, gg.sexdists, nrow=2, ncol=1, labels=c("A", "B"))
dev.off()

##### Set up groups #####
# apos given their MLG designation
# sexuals are based on general location / NJ tree
AllPops.gc@other$group <- "blank"
AllPops.gc@other$group[1:18] <- "CO-S"
AllPops.gc@other$group[7:9] <- "Laramie-S"
AllPops.gc@other$group[16:18] <- "YK-S"
AllPops.gc@other$group[19:41] <- "Laramie-S"
AllPops.gc@other$group[37] <- "L45-S"
AllPops.gc@other$group[42:46] <- "ND-A"
AllPops.gc@other$group[47:51] <- "MT.big-A"
AllPops.gc@other$group[52:56] <- "BC-A"
AllPops.gc@other$group[57:61] <- "MT.big-A"
AllPops.gc@other$group[62:76] <- "MT.big-A"
AllPops.gc@other$group[77:81] <- "Laramie.big-A"
AllPops.gc@other$group[c(82,83,85)] <- "Laramie.small-A"
AllPops.gc@other$group[c(84,86,87,88,89,90,91)] <- "Laramie.big-A"
AllPops.gc@other$group[c(92,93,94,96)] <- "L39.big-A"
AllPops.gc@other$group[95] <- "L39.small-A"
AllPops.gc@other$group[97:100] <- "MT.big-A"
AllPops.gc@other$group[101:104] <- "MT.big-A"
AllPops.gc@other$group[105] <- "L62-A"
AllPops.gc@other$group[106:109] <- "SK-A"
AllPops.gc@other$group[110:114] <- "YK-A"

AllPops.gc@other$group <- factor(AllPops.gc@other$group, levels = c(
  "CO-S",
  "Laramie-S",
  "L62-A",
  "Laramie.big-A",
  "Laramie.small-A",
  "L45-S",
  "L39.big-A",
  "L39.small-A",
  "MT.big-A",
  "ND-A",
  "SK-A",
  "BC-A",
  "YK-A",
  "YK-S"
))

AllPops.gc@other$group2 <- factor(AllPops.gc@other$group, levels = c(
  "CO-S",
  "Laramie-S",
  "L45-S",
  "YK-S",
  "L62-A",
  "Laramie.big-A",
  "Laramie.small-A",
  "L39.big-A",
  "L39.small-A",
  "MT.big-A",
  "ND-A",
  "SK-A",
  "BC-A",
  "YK-A"
))

