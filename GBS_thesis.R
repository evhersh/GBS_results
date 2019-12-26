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
library(ggtree)
library(ggnewscale)


###############
# Data import #
###############

# strata (sample, pop, ms)
mystrata <- read.csv("~/Google Drive/GitHub/Hookeri-GBS/popmap_all.csv")

# vcf to vcfR
# vcf <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/final.filtered.snps.vcf")
# vcf.dips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.dips.vcf")
# vcf.trips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.trips.vcf")
# vcf.tets <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.tets.vcf")

# vcfR to genind to genclone
dips.gi <- vcfR2genind(vcf.dips, sep = "/", ploidy=2, return.alleles = TRUE)
#dips.gi <- vcfR2genind(vcf.dips, sep = "/", ploidy=2)
dips.gc <- as.genclone(dips.gi)
sampleorder <- match(indNames(dips.gc), mystrata$id)
strata(dips.gc) <- mystrata[sampleorder,]
setPop(dips.gc) <- ~pop

trips.gi <- vcfR2genind(vcf.trips, sep = "/", ploidy=3, return.alleles = TRUE)
#trips.gi <- vcfR2genind(vcf.trips, sep = "/", ploidy=3)
trips.gl <- vcfR2genlight(vcf.trips)
#ploidy(trips.gi) <- 3
trips.gc <- as.genclone(trips.gi)
sampleorder <- match(indNames(trips.gc), mystrata$id)
strata(trips.gc) <- mystrata[sampleorder,]
setPop(trips.gc) <- ~pop

tets.gi <- vcfR2genind(vcf.tets, sep = "/", ploidy=4, return.alleles = TRUE)
#tets.gi <- vcfR2genind(vcf.tets, sep = "/", ploidy=4)
tets.gc <- as.genclone(tets.gi)
sampleorder <- match(indNames(tets.gc), mystrata$id)
strata(tets.gc) <- mystrata[sampleorder,]
setPop(tets.gc) <- ~pop

# combine genclones
dipsNtripsNtets.gc <- repool(dips.gc,trips.gc,tets.gc)
#dipsNtripsNtets.gc$pop <- factor(dipsNtripsNtets.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))

AllPops.gc <-as.genclone(dipsNtripsNtets.gc)
AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))

apo.list <- c("C23-A", "C27-A", "C43-A", "C85-A", "C86-A", "C87-A", "C88-A", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A", "L45-A", "L62-A", "S03-A", "SM-A")
sex.list <- c("B42-S", "B46-S", "B49-S", "B53-S", "B60-S", "C59-S", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L45-S", "L62-S")

dips.gc <- popsub(AllPops.gc, sublist=sex.list, drop=FALSE)
apos.gc <- popsub(AllPops.gc, sublist=apo.list, drop=FALSE)

GBS.cols <- 
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

# apodist.df %>% filter(between(value, 0.05, 0.1))
# apodist.df %>% filter(between(value, 0.07, 0.1))
# apodist.df %>% filter(between(value, 0.1, 0.15))
# apodist.df %>% filter(between(value, 0.25, 0.4))

gg.MLGdists <- ggplot(data=apodist.df, aes(x=value))+
  geom_histogram(aes(y=..density..), bins=100, color="black", fill="slateblue", size=.75)+
  theme_classic()+
  geom_density(alpha=.3)+
  geom_vline(aes(xintercept=0.02948973), color="black", linetype="twodash", size=.75)+
  labs(x="", y="Density")+
  theme(legend.position = "none")+
  xlim(0, 0.35)

#sexual pairwise dist prevosti

ddist.dips.df <- melt(as.matrix(ddist.dips), varnames=c("row", "col"))
sexdist.df <- ddist.dips.df[as.numeric(ddist.dips.df$row) > as.numeric(ddist.dips.df$col),]

gg.sexdists <- ggplot(data=sexdist.df, aes(x=value))+
  geom_histogram(aes(y=..density..), bins=100, color="black",fill="gray48", size=.75)+
  theme_classic()+
  geom_density(alpha=.3)+
  labs(x="Pairwise prevosti distance", y="Density")+
  theme(legend.position = "none")+
  xlim(0, 0.35)

#png("distance.png", height=7, width=8, res=300, units="in")
ggarrange(gg.MLGdists, gg.sexdists, nrow=2, ncol=1, labels=c("A", "B"))
#dev.off()

###########
# K-means #
###########

grp2 <- find.clusters(apos.gc, max.n.clust=30) # looks like lowest BIC is 13 clusters...
names(grp2)
table(pop(apos.gc), grp2$grp)
table.value(table(pop(apos.gc), grp2$grp), col.lab=paste("inf", 1:8),
            row.lab=paste("ori", 1:8))

## from vcfR -- run k-means multiple times
# library(vcfR)
# vcf <- read.vcfR("prubi_gbs.vcf.gz")
# pop.data <- read.table("population_data.gbs.txt", sep = "\t", header = TRUE)
# all(colnames(vcf@gt)[-1] == pop.data$AccessID)
# ## [1] TRUE
# gl_rubi <- vcfR2genlight(vcf)

# library(adegenet)
maxK2 <- 30
myMat2 <- matrix(nrow=50, ncol=maxK2)
colnames(myMat2) <- 1:ncol(myMat2)
for(i in 1:nrow(myMat2)){
  grp2 <- find.clusters(apos.gc, n.pca = 200, choose.n.clust = FALSE,  max.n.clust = maxK2)
  myMat2[i,] <- grp2$Kstat
}

my_df3 <- melt(myMat2)
colnames(my_df3)[1:3] <- c("Group", "K", "BIC")
my_df3$K <- as.factor(my_df3$K)
head(my_df3)

p01 <- ggplot(my_df3, aes(x = K, y = BIC))
p01 <- p01 + geom_boxplot()
p01 <- p01 + theme_bw()
p01 <- p01 + xlab("Number of groups (K)")
p01 <- p01 + theme(axis.text=element_text(size=6))
p01 # 12 clusters kind of looks right?

# plot 2 #
my_k00 <- 8

grp20_l <- vector(mode = "list", length = length(my_k00))
dapc00_l <- vector(mode = "list", length = length(my_k00))

for(i in 1:length(dapc00_l)){
  set.seed(70)
  grp20_l[[i]] <- find.clusters(apos.gc, n.pca = 200, n.clust = my_k00[i])
  dapc00_l[[i]] <- dapc(apos.gc, pop = grp20_l[[i]]$grp, n.pca = 8, n.da = my_k00[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp2_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df00 <- as.data.frame(dapc00_l[[ length(dapc00_l) ]]$ind.coord)
my_df00$Group <- dapc00_l[[ length(dapc00_l) ]]$grp
head(my_df00)

my_df00$Group <- as.character(my_df00$Group)

subset(my_df01, Posterior==1 & K==8)
# group A = L06, L16 (laramie big)
# group B = C27, MT big
# group C = L39-A_4, L62-A_5
# group D = L39-A_1,2,3,5
# group E = L16-A_1,2,4
# group F = C43
# group G = C23
# group H = SM
# k8
# G-B
my_df00[my_df00$Group == 1, "Group"] <- "B"
my_df00[my_df00$Group == 2, "Group"] <- "A"
my_df00[my_df00$Group == 3, "Group"] <- "D"
my_df00[my_df00$Group == 4, "Group"] <- "E"
my_df00[my_df00$Group == 5, "Group"] <- "H"
my_df00[my_df00$Group == 6, "Group"] <- "F"
my_df00[my_df00$Group == 7, "Group"] <- "G"
my_df00[my_df00$Group == 8, "Group"] <- "C"



p02 <- ggplot(my_df00, aes(x = LD1, y = LD2, fill = Group))
p02 <- p02 + geom_point(size = 3, shape = 21)
p02 <- p02 + theme_bw()
p02 <- p02 + scale_color_brewer(palette="Paired")
p02 <- p02 + scale_fill_brewer(palette="Paired")
#p02 <- p02 + scale_color_manual(values=c(col_vector))
#p02 <- p02 + scale_fill_manual(values=c(paste(col_vector, "66", sep = "")))
p02


# plot 3
my_k01 <- 8:12

grp2_l <- vector(mode = "list", length = length(my_k01))
dapc01_l <- vector(mode = "list", length = length(my_k01))

for(i in 1:length(dapc01_l)){
  set.seed(10)
  grp2_l[[i]] <- find.clusters(apos.gc, n.pca = 200, n.clust = my_k01[i])
  dapc01_l[[i]] <- dapc(apos.gc, pop = grp2_l[[i]]$grp, n.pca = 8, n.da = my_k01[i])
  #  dapc01_l[[i]] <- dapc(gl_rubi, pop = grp2_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df01 <- as.data.frame(dapc01_l[[ length(dapc01_l) ]]$ind.coord)
my_df01$Group <- dapc01_l[[ length(dapc01_l) ]]$grp
head(my_df01)

tmp01 <- as.data.frame(dapc01_l[[1]]$posterior)
tmp01$K <- my_k01[1]
tmp01$Sample <- rownames(tmp01)
tmp01 <- melt(tmp01, id = c("Sample", "K"))
names(tmp01)[3:4] <- c("Group", "Posterior")
tmp01$pop <- mystrata$pop[match(tmp01$Sample, mystrata$id)]
my_df01 <- tmp01

for(i in 2:length(dapc01_l)){
  tmp01 <- as.data.frame(dapc01_l[[i]]$posterior)
  tmp01$K <- my_k01[i]
  tmp01$Sample <- rownames(tmp01)
  tmp01 <- melt(tmp01, id = c("Sample", "K"))
  names(tmp01)[3:4] <- c("Group", "Posterior")
  tmp01$pop <- mystrata$pop[match(tmp01$Sample, mystrata$id)]
  my_df01 <- rbind(my_df01, tmp01)
}

grp2.labs <- paste("K =", my_k01)
names(grp2.labs) <- my_k01

my_df01$pop <- factor(my_df01$pop, levels=c("L62-A","L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A"))

# run first with both uncommented, then comment top one (and run second line) when re-making the colors
#my_df001 <- my_df01
my_df01 <- my_df001


my_df01$Group <- as.character(my_df01$Group)

# G-B, C-G
my_df01[ my_df01$K == 8 & my_df01$Group == 1, "Group"] <- "A"
my_df01[ my_df01$K == 8 & my_df01$Group == 2, "Group"] <- "C"
my_df01[ my_df01$K == 8 & my_df01$Group == 3, "Group"] <- "G"
my_df01[ my_df01$K == 8 & my_df01$Group == 4, "Group"] <- "D"
my_df01[ my_df01$K == 8 & my_df01$Group == 5, "Group"] <- "E"
my_df01[ my_df01$K == 8 & my_df01$Group == 6, "Group"] <- "F"
my_df01[ my_df01$K == 8 & my_df01$Group == 7, "Group"] <- "B"
my_df01[ my_df01$K == 8 & my_df01$Group == 8, "Group"] <- "H"

# G B A F C D I H E
# A-C, C-I, C-E, I-G, G-B, C-G
my_df01[ my_df01$K == 9 & my_df01$Group == 1, "Group"] <- "B"
my_df01[ my_df01$K == 9 & my_df01$Group == 2, "Group"] <- "C"
my_df01[ my_df01$K == 9 & my_df01$Group == 3, "Group"] <- "A"
my_df01[ my_df01$K == 9 & my_df01$Group == 4, "Group"] <- "F"
my_df01[ my_df01$K == 9 & my_df01$Group == 5, "Group"] <- "G"
my_df01[ my_df01$K == 9 & my_df01$Group == 6, "Group"] <- "D"
my_df01[ my_df01$K == 9 & my_df01$Group == 7, "Group"] <- "I"
my_df01[ my_df01$K == 9 & my_df01$Group == 8, "Group"] <- "H"
my_df01[ my_df01$K == 9 & my_df01$Group == 9, "Group"] <- "E"

# J A G B E C D I F H
# B-A, D-G, G-B, J-H, G-C, I-J, C-J, C-F, G-B, C-G
my_df01[ my_df01$K == 10 & my_df01$Group == 1, "Group"] <- "J"
my_df01[ my_df01$K == 10 & my_df01$Group == 2, "Group"] <- "A"
my_df01[ my_df01$K == 10 & my_df01$Group == 3, "Group"] <- "B"
my_df01[ my_df01$K == 10 & my_df01$Group == 4, "Group"] <- "C"
my_df01[ my_df01$K == 10 & my_df01$Group == 5, "Group"] <- "E"
my_df01[ my_df01$K == 10 & my_df01$Group == 6, "Group"] <- "G"
my_df01[ my_df01$K == 10 & my_df01$Group == 7, "Group"] <- "D"
my_df01[ my_df01$K == 10 & my_df01$Group == 8, "Group"] <- "I"
my_df01[ my_df01$K == 10 & my_df01$Group == 9, "Group"] <- "F"
my_df01[ my_df01$K == 10 & my_df01$Group == 10, "Group"] <- "H"

# F B C D E A G H I J K
# F-A, F-H, C-F, C-I, K-G, B-I, D-J, E-I, I-B, G-B, C-G
my_df01[ my_df01$K == 11 & my_df01$Group == 1, "Group"] <- "H"
my_df01[ my_df01$K == 11 & my_df01$Group == 2, "Group"] <- "I"
my_df01[ my_df01$K == 11 & my_df01$Group == 3, "Group"] <- "F"
my_df01[ my_df01$K == 11 & my_df01$Group == 4, "Group"] <- "J"
my_df01[ my_df01$K == 11 & my_df01$Group == 5, "Group"] <- "C"
my_df01[ my_df01$K == 11 & my_df01$Group == 6, "Group"] <- "A"
my_df01[ my_df01$K == 11 & my_df01$Group == 7, "Group"] <- "K"
my_df01[ my_df01$K == 11 & my_df01$Group == 8, "Group"] <- "E"
my_df01[ my_df01$K == 11 & my_df01$Group == 9, "Group"] <- "G"
my_df01[ my_df01$K == 11 & my_df01$Group == 10, "Group"] <- "D"
my_df01[ my_df01$K == 11 & my_df01$Group == 11, "Group"] <- "B"


# B-J, I-G, E-K, H-L, C-G
my_df01[ my_df01$K == 12 & my_df01$Group == 1, "Group"] <- "A"
my_df01[ my_df01$K == 12 & my_df01$Group == 2, "Group"] <- "J"
my_df01[ my_df01$K == 12 & my_df01$Group == 3, "Group"] <- "G"
my_df01[ my_df01$K == 12 & my_df01$Group == 4, "Group"] <- "D"
my_df01[ my_df01$K == 12 & my_df01$Group == 5, "Group"] <- "K"
my_df01[ my_df01$K == 12 & my_df01$Group == 6, "Group"] <- "F"
my_df01[ my_df01$K == 12 & my_df01$Group == 7, "Group"] <- "I"
my_df01[ my_df01$K == 12 & my_df01$Group == 8, "Group"] <- "L"
my_df01[ my_df01$K == 12 & my_df01$Group == 9, "Group"] <- "C"
my_df01[ my_df01$K == 12 & my_df01$Group == 10, "Group"] <- "B"
my_df01[ my_df01$K == 12 & my_df01$Group == 11, "Group"] <- "E"
my_df01[ my_df01$K == 12 & my_df01$Group == 12, "Group"] <- "H"



p03 <- ggplot(my_df01, aes(x = Sample, y = Posterior, fill = Group))
p03 <- p03 + geom_bar(stat = "identity")
p03 <- p03 + facet_grid(K ~ pop, scales = "free_x", space = "free", 
                        labeller = labeller(K = grp2.labs))
p03 <- p03 + theme_bw()
p03 <- p03 + ylab("Posterior membership probability")
p03 <- p03 + theme(legend.position='none')
p03 <- p03 + scale_color_brewer(palette="Paired")
p03 <- p03 + scale_fill_brewer(palette="Paired")
p03 <- p03 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),panel.spacing.x=unit(0.1, "lines"))
#p03

png("kmeans.png", height=7, width=9, res=300, units="in")
ggarrange(ggarrange(p01,
                    p02,
                    ncol = 2, labels = c("A", "B")),
          p03,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)
dev.off()

########
# DAPC #
########

# optimize number of PCs to keep
dapc.x <- dapc(AllPops.gc, n.da=100, n.pca=50)
temp <- optim.a.score(dapc.x) #11 is the optimal number of PCs

# xval
xtab <- tab(AllPops.gc, NA.method="mean")
grp <- pop(AllPops.gc)

xval <- xvalDapc(xtab, grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval[2:6] # 20 PCs has the highest prediction and lowest error

# make the dapc
# set shapes - triangle for sexuals, circles for apos
my.pch <- c(17, 17, 17, 17, 17, 17, 21, 17, 17, 17, 17, 17, 17, 21, 21, 21, 21, 21, 17, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 17)
my.pch <-replace(my.pch,my.pch==21, 19)
my.pch.sub <- my.pch[-c(30,29,28,26)]
# for ms
setPop(AllPops.gc) <- ~ms
AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("S", "A"))
hookeri.dapc.ms <- dapc(AllPops.gc, grp=AllPops.gc$grp, n.pca=20, n.da=100)
scatter(hookeri.dapc.ms, grp = AllPops.gc$pop, cex = 2, legend = TRUE, clabel = T, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)

# all pops, but color by ms
setPop(AllPops.gc) <- ~pop
hookeri.dapc.msp <- dapc(AllPops.gc, grp=AllPops.gc$grp, n.pca=20, n.da=100)
scatter(hookeri.dapc.msp, grp = AllPops.gc$strata$ms, cex = 2, legend = TRUE, clabel = T,cstar=0, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch=c(17,19), col=cols.ms)



# for pops (all)#


setPop(AllPops.gc) <- ~pop
AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))
hookeri.dapc <- dapc(AllPops.gc, grp=AllPops.gc$grp, n.pca=20, n.da=100)
scatter(hookeri.dapc, grp = AllPops.gc$pop, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch=my.pch)

DAPC.all.df <- as.data.frame(hookeri.dapc$ind.coord)
DAPC.all.df$Group <- hookeri.dapc$grp
head(DAPC.all.df)

DAPC.all.df$Group <- as.character(DAPC.all.df$Group)
DAPC.all.df$Group <- factor(DAPC.all.df$Group, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))
DAPC.all.df <- separate(DAPC.all.df, "Group", sep= "-", c("pop", "ms"))

DAPC.all.df$group <- "blank"
DAPC.all.df$group[1:18] <- "CO-S"
DAPC.all.df$group[7:9] <- "Laramie-S"
DAPC.all.df$group[16:18] <- "YK-S"
DAPC.all.df$group[19:41] <- "Laramie-S"
DAPC.all.df$group[37] <- "L45-S"
DAPC.all.df$group[42:46] <- "ND-A"
DAPC.all.df$group[47:51] <- "MT.big-A"
DAPC.all.df$group[52:56] <- "BC-A"
DAPC.all.df$group[57:61] <- "MT.big-A"
DAPC.all.df$group[62:76] <- "MT.big-A"
DAPC.all.df$group[77:81] <- "Laramie.big-A"
DAPC.all.df$group[c(82,83,85)] <- "Laramie.small-A"
DAPC.all.df$group[c(84,86,87,88,89,90,91)] <- "Laramie.big-A"
DAPC.all.df$group[c(92,93,94,96)] <- "L39.big-A"
DAPC.all.df$group[95] <- "L39.small-A"
DAPC.all.df$group[97:100] <- "MT.big-A"
DAPC.all.df$group[101:104] <- "MT.big-A"
DAPC.all.df$group[105] <- "L62-A"
DAPC.all.df$group[106:109] <- "SK-A"
DAPC.all.df$group[110:114] <- "YK-A"

DAPC.all.df$group <- factor(DAPC.all.df$group, levels=c("CO-S", "Laramie-S", "L62-A", "Laramie.big-A", "Laramie.small-A", "L45-S", "L39.big-A", "L39.small-A", "MT.big-A", "ND-A", "SK-A", "BC-A", "YK-A", "YK-S"))

levels(DAPC.all.df$group)

DAPC.cols = c("gray27", "gray", "#CAB2D6", "#A6CEE3", "#FFFF99", "white", "#33A02C", "#B2DF8A", "#1F78B4", "#FF7F00", "#E31A1C", "#FB9A99", "#B15928", "#FDBF6F")

brewer.pal(12, "Paired")
# [1] "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C"
# [7] "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"


DAPC.gg <- ggplot(DAPC.all.df, aes(x = LD1, y = LD2))+ 
  geom_point(aes(fill=group, shape=ms), size=2.5, stroke=0.7)+ 
  theme_bw()+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values=DAPC.cols)+
  #geom_text(aes(label=pop), size=3)+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 24, 21, 21, 21, 24, 21, 21, 21, 21, 21, 21, 21, 24)), list(size=0.25)), shape=FALSE)+
  theme(legend.title=element_blank(), legend.text = element_text(size=5), legend.key.size = unit(0.3, 'cm'), legend.position=c(0.925, 0.2), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"))

png("dapc.png", height=6, width=8, res=300, units="in")
DAPC.gg
dev.off()


# tiff("dapc.all.tiff", height = 6, width = 8, res = 300, units="in")
# ggplot(DAPC.all.df, aes(x = LD1, y = LD2))+ 
#   geom_point(aes(fill=group, shape=ms), size=3, stroke=0.7)+ 
#   theme_bw()+
#   scale_shape_manual(values = c(21, 24))+
#   scale_fill_manual(values=DAPC.cols)+
#   #geom_text(aes(label=pop), size=3)+
#   guides(fill = guide_legend(override.aes=list(shape=c(24, 24, 21, 21, 21, 24, 21, 21, 21, 21, 21, 21, 21, 24))))
# dev.off()

DAPC.gg.sub <- ggplot(DAPC.all.df, aes(x = LD1, y = LD2))+ 
  geom_point(aes(fill=group, shape=ms), size=3, stroke=0.7)+ 
  theme_bw()+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values=DAPC.cols)+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 24, 21, 21, 21, 24, 21, 21, 21, 21, 21, 21, 21, 24))))+
  #geom_text(aes(label=pop), size=3)+
  ylim(-50,50)+
  xlim(-50,25)

# tiff("dapc.sub.tiff", height = 6, width = 8, res = 300, units="in")
# ggplot(DAPC.all.df, aes(x = LD1, y = LD2))+ 
#   geom_point(aes(fill=group, shape=ms), size=3, stroke=0.7)+ 
#   theme_bw()+
#   scale_shape_manual(values = c(21, 24))+
#   scale_fill_manual(values=DAPC.cols)+
#   guides(fill = guide_legend(override.aes=list(shape=c(24, 24, 21, 21, 21, 24, 21, 21, 21, 21, 21, 21, 21, 24))))+
#   #geom_text(aes(label=pop), size=3)+
#   ylim(-50,50)+
#   xlim(-50,25)
# dev.off()

DAPC.gg.sub2 <- ggplot(DAPC.all.df, aes(x = LD1, y = LD2))+ 
  geom_point(aes(fill=group, shape=ms), size=2)+ 
  theme_bw()+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values=DAPC.cols)+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 24, 21, 21, 21, 24, 21, 21, 21, 21, 21, 21, 21, 24))))+
  #geom_text(aes(label=pop), size=3)+
  ylim(0,25)+
  xlim(0,25)

#########
# Trees #
#########

# main tree with all individuals
hookeri.nj <- aboot(AllPops.gc, dist = provesti.dist, sample = 200, tree = "nj", cutoff = 50, quiet = TRUE)

#png("njtree_all.png", height=7, width=7, res=300, units="in")
plot.phylo(hookeri.nj, type="unrooted", cex=0.6, lab4ut = "axial", font=2, show.tip.label = FALSE, no.margin = TRUE)
tiplabels(pch=21, col="black", bg=DAPC.cols[DAPC.all.df$group])
add.scale.bar()
legend("bottomright", legend=levels(DAPC.all.df$group), cex=0.5, pch=21, col="black", pt.bg=DAPC.cols, pt.cex=1)
#dev.off()

DAPC.cols
# pop trees based on Nei's distance

Nei.mat <- dist.genpop(AllPops.gp, method=1)
Nei.tree <- nj(Nei.mat)

Nei.mat.dips <- dist.genpop(dips.gp, method=1)
Nei.tree.dips <- nj(Nei.mat.dips)

sexpop.cols <- c("black", "black", "black", "black", "black", "#FDBF6F", "gray49", "gray49", "gray49", "gray49", "gray49", "gray49", "gray49", "gray49")
# only diploids

png("neitree_dips.png", height=7, width=7, res=300, units="in")
plot.phylo(Nei.tree.dips, type="unrooted", tip.col=sexpop.cols, cex=0.6, lab4ut = "axial", font=2, show.tip.label = TRUE, no.margin = TRUE, cex=2)
#edgelabels(annot[annot>0], which(annot>0), frame="n")
add.scale.bar()
dev.off()

pop.cols<- c("black", "black", "black", "black", "black", "gray49", "#CAB2D6", "gray49", "gray49", "gray49", "gray49", "gray49", "gray49", "#A6CEE3", "#FF7F00", "#A6CEE3", "#33A02C", "#1F78B4", "gray49", "#1F78B4", "#1F78B4", "#1F78B4", "#1F78B4", "#1F78B4", "#1F78B4", "#FF7F00", "#FB9A99", "#E31A1C", "#B15928", "#FDBF6F")
# all pops
png("neitree.png", height=7, width=7, res=300, units="in")
plot.phylo(Nei.tree, type="unrooted", tip.col=pop.cols, cex=0.6, lab4ut = "axial", font=2, show.tip.label = TRUE, no.margin = TRUE)
#annot <- round(Nei.tree$edge.length,2)
#edgelabels(annot[annot>0], which(annot>0), frame="n")
add.scale.bar()
dev.off()

###########################
# Observed heterozygosity #
###########################

ind.Hobs <- data.frame(matrix(NA, nrow=114, ncol=3))
colnames(ind.Hobs) <- c("ind", "ploidy", "Hobs")
ind.Hobs$ms <- "NA"
ind.Hobs$pop <- "NA"
ind.Hobs$mlg <- "NA"
AllPops.gc[1]@mlg

for(i in 1:114) {
  ind.Hobs[i,5] <- as.character(AllPops.gc[i]@strata$pop)
  ind.Hobs[i,4] <- as.character(AllPops.gc[i]@strata$ms)
  #ind.Hobs[i,3] <- mean(summary(AllPops.gc[i])$Hobs, na.rm=TRUE)
  #ind.Hobs[i,2] <- ploidy(AllPops.gc[i])
  #ind.Hobs[i,1] <- indNames(AllPops.gc[i])
  print(i)
}

ind.Hobs$ms <- factor(ind.Hobs$ms, levels=c("S", "A"))
ind.Hobs$group <- DAPC.all.df$group

mean.Hobs <- ind.Hobs %>%
  group_by(ms) %>%
  summarise(mean=mean(Hobs), sd=sd(Hobs))

DAPC.all.df$group

gg.Hobs <- ggplot()+
  geom_boxplot(data=ind.Hobs, aes(y=Hobs, x=ms), alpha=0, width=0.2, outlier.shape=NA)+
  geom_jitter(data=ind.Hobs, aes(x=ms, y=Hobs, fill=group), alpha=1, shape=21, width=0.05, size=2)+
  scale_fill_manual(values=DAPC.cols, name="Mating system")+
  guides(fill=FALSE)+
  new_scale_fill()+
  geom_errorbar(data=mean.Hobs, aes(ymax=mean+sd, ymin=mean-sd, x=as.numeric(as.factor(ms))+0.4), width=0)+
  geom_point(data=mean.Hobs, aes(x=as.numeric(as.factor(ms))+0.4, y=mean, fill=ms), shape=21, size=3.5)+
  scale_fill_manual(values=c("white", "white"))+
  labs(x="Mating system", y="Mean observed heterozygosity")+
  theme_classic()+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))+
  theme(legend.title=element_blank(), legend.text = element_text(size=5), legend.key.size = unit(0.3, 'cm'), legend.position=c(0.925, 0.2), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"))



png("ind_het_ms.png", height=7, width=8, res=300, units="in")
gg.Hobs
dev.off()

##########
# ploidy #
##########

# mess with CNV
gt<-extract.gt(vcf)
hets<-is_het(gt)
is.na(vcf@gt[,-1][!hets]) <- TRUE

# try to clean up
ad <- extract.gt(vcf, element = "AD")
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)



freq1 <- ad1/(ad1+ad2)
freq2 <- ad2/(ad1+ad2)
myPeaks1 <- freq_peak(freq1, getPOS(vcf))
is.na(myPeaks1$peaks[myPeaks1$counts < 20]) <- TRUE
myPeaks2 <- freq_peak(freq2, getPOS(vcf), lhs = FALSE)
is.na(myPeaks2$peaks[myPeaks2$counts < 20]) <- TRUE
myPeaks1
peak_to_ploid(myPeaks1)
freq_peak_plot(pos = getPOS(vcf), ab1 = freq1, ab2 = freq2, fp1 = myPeaks1, fp2=myPeaks2)

# vcfR Knaus
# determine ploidy
knitr::kable(vcf@gt[c(1:2,11,30),1:4])
ad <- extract.gt(vcf, element = 'AD')
knitr::kable(ad[c(1:2,11,30),1:4])
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)

ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)

"indianred"
"tan2"
# diploid !!!!
hist(ad2[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# diploid
hist(ad2[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n", main="Diploid", xlab="")
hist(ad1[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))

tiff("ABplotDip.tiff", res=300, units="in", width=8.6, height=5.8)
hist(ad2[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n", main="Diploid", xlab="")
hist(ad1[,"B42-S_1"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
dev.off()

# triploid !!!!
hist(ad2[,"L17-A_2"], breaks = seq(0,1,by=0.02), col = "indianred", xaxt="n", main="Triploid")
hist(ad1[,"L17-A_2"], breaks = seq(0,1,by=0.02), col = "tan2", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

tiff("ABplotTrip.tiff", res=300, units="in", width=8.6, height=5.8)
hist(ad2[,"L17-A_1"], breaks = seq(0,1,by=0.02), col = "indianred", xaxt="n", main="Triploid")
hist(ad1[,"L17-A_1"], breaks = seq(0,1,by=0.02), col = "tan2", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))
dev.off()

# Dan's
hist(ad2[,"L62-A_5"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"L62-A_5"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# L39
hist(ad2[,"L39-A_4"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"L39-A_4"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# tetraploid???
tiff("ABplotTet.tiff", res=300, units="in", width=8.6, height=5.8)
hist(ad2[,"S03-A_4"], breaks = seq(0,1,by=0.02), col = "indianred", xaxt="n", main="Tetraploid", xlab="Allele balance")
hist(ad1[,"S03-A_4"], breaks = seq(0,1,by=0.02), col = "tan2", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))
dev.off()

# L16-A_2
hist(ad2[,"L16-A_5"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"L16-A_5"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# L45-A_1 (diploid in a triploid maze...)
hist(ad2[,"L45-S_1"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"L45-S_1"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# C59 (ind 3 looks diploid, the others look weird...maybe remove from analysis?)
hist(ad2[,"C59-S_2"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"C59-S_2"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))

# L11 S_1
hist(ad2[,"C87-A_1"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"C87-A_1"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))


