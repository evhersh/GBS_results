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
  geom_histogram(aes(y=..density..), bins=100, color="black", fill="darkgrey", size=.75)+
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
  geom_histogram(aes(y=..density..), bins=100, color="black",fill="white", size=.75)+
  theme_classic()+
  geom_density(fill="white", alpha=.2)+
  labs(x="Pairwise prevosti distance", y="Density")+
  theme(legend.position = "none")+
  xlim(0, 0.35)

png("distance_apos.png", height=5, width=8, res=300, units="in")
gg.MLGdists
dev.off()

png("distance.png", height=7, width=8, res=300, units="in")
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

##### set up colors #####
group.cols2 <- c("CO-S" = "#E31A1C",
                 "Laramie-S" = "sienna2",
                 "L62-A"="#CAB2D6",
                 "Laramie.big-A"="#A6CEE3",
                 "Laramie.small-A"="cyan",
                 "L45-S"="khaki",
                 "L39.big-A"="darkorchid3",
                 "L39.small-A"="deeppink3",
                 "MT.big-A"="#1F78B4",
                 "ND-A"="olivedrab",
                 "SK-A"="lawngreen",
                 "BC-A"="darkgreen",
                 "YK-A"="lightseagreen",
                 "YK-S"="goldenrod")

group.cols3 <- c("CO-S" = "#E65100",
                 "Laramie-S" = "#FFB300",
                 
                 "L62-A"="#E1BEE7",
                 "Laramie.big-A"="#7E57C2",
                 "Laramie.small-A"="#4A14BC",
                 "L45-S"="#9c27b0",
                 "L39.big-A"="#AD1457",
                 "L39.small-A"="#EC407A",
                 "MT.big-A"="#1565C0",
                 "ND-A"="#64B5F6",
                 
                 "SK-A"="#004d40",
                 "BC-A"="#26A69A",
                 "YK-A"="#558B2F",
                 "YK-S"="#a5d6a7")

##### NJ TREE #####
hookeri.nj1000 <- aboot(AllPops.gc, dist = provesti.dist, sample = 1000, tree = "nj", cutoff = 50, quiet = TRUE)

saveRDS(hookeri.nj1000, here("data", "hookeri_nj1000.rds"))

hookeri.nj1000$node.label2 <- as.numeric(hookeri.nj1000$node.label)
hookeri.nj1000$node.label2[hookeri.nj1000$node.label<90] <- 1
hookeri.nj1000$node.label2[hookeri.nj1000$node.label>90] <- 2
hookeri.nj1000$node.label2[is.na(hookeri.nj1000$node.label2)] <- 1

tree.pch <- group.cols3[AllPops.gc@other$group]
tree.pch[1:41] <- 24
tree.pch[42:114] <- 21
tree.pch <- as.numeric(tree.pch)

tree.legend.pch <- c(24, 24, 21, 21, 21, 24, 21, 21, 21, 21, 21, 21, 21, 24)

png("njtree_all_v2.png", height=7, width=7, res=300, units="in")
plot.phylo(hookeri.nj1000, type="unrooted", cex=0.6, lab4ut = "axial", font=2, show.tip.label = FALSE, no.margin = TRUE)
tiplabels(pch=tree.pch, col="black", bg=group.cols3[AllPops.gc@other$group], cex=1.3)
add.scale.bar()
legend("bottomright", legend=levels(AllPops.gc@other$group), cex=0.8, pch=tree.legend.pch, col="black", pt.bg=group.cols3, pt.cex=1.7, inset=0.01)
dev.off()


##### Observed heterozygosity #####

ind.Hobs <- data.frame(matrix(NA, nrow=114, ncol=3))
colnames(ind.Hobs) <- c("ind", "ploidy", "Hobs")
ind.Hobs$ms <- "NA"
ind.Hobs$pop <- "NA"
ind.Hobs$mlg <- "NA"
AllPops.gc[1]@mlg

for(i in 1:114) {
  ind.Hobs[i,5] <- as.character(AllPops.gc[i]@strata$pop)
  ind.Hobs[i,4] <- as.character(AllPops.gc[i]@strata$ms)
  ind.Hobs[i,3] <- mean(summary(AllPops.gc[i])$Hobs, na.rm=TRUE)
  ind.Hobs[i,2] <- ploidy(AllPops.gc[i])
  ind.Hobs[i,1] <- indNames(AllPops.gc[i])
  print(i)
}

ind.Hobs$ms <- factor(ind.Hobs$ms, levels=c("S", "A"))
ind.Hobs$group <- AllPops.gc@other$group

mean.Hobs <- ind.Hobs %>%
  group_by(ms) %>%
  summarise(mean=mean(Hobs), sd=sd(Hobs))

gg.Hobs <- ggplot()+
  geom_boxplot(data=ind.Hobs, aes(y=Hobs, x=ms), alpha=0, width=0.2, outlier.shape=NA)+
  geom_jitter(data=ind.Hobs, aes(x=ms, y=Hobs, fill=group), alpha=1, shape=21, width=0.05, size=2.5)+
  scale_fill_manual(values=group.cols3, name="Mating system")+
  #guides(fill="none")+
  new_scale_fill()+
  geom_errorbar(data=mean.Hobs, aes(ymax=mean+sd, ymin=mean-sd, x=as.numeric(as.factor(ms))+0.4), width=0)+
  geom_point(data=mean.Hobs, aes(x=as.numeric(as.factor(ms))+0.4, y=mean, fill=ms), shape=21, size=3.5, show.legend = FALSE)+
  scale_fill_manual(values=c("white", "white"))+
  labs(x="Mating system", y="Mean observed heterozygosity")+
  theme_bw()+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))+
  theme(legend.title=element_blank(), 
        legend.text = element_text(size=10), 
        legend.key.size = unit(0.3, 'cm'),
        legend.position=c(0.9, 0.2), 
        legend.background = element_blank(), 
        legend.box.background = element_rect(colour = "black"))


png("ind_het_ms.png", height=7, width=8, res=300, units="in")
gg.Hobs
dev.off()


