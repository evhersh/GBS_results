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
library(viridis)
library(RColorBrewer)

# Apos only

kmeans.apos <- find.clusters(apos.gc, max.n.clust=30) # looks like lowest BIC is 13 clusters...
apos.tab <- tab(apos.gc, NA.method="mean")
apos.xval <- xvalDapc(apos.tab, grp=kmeans.apos$grp)
apos.xval[2:6]

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
p01 # 12 clusters kind of looks right?

# plot 2 #
my_k00 <- 8

grp20_l <- vector(mode = "list", length = length(my_k00))
dapc00_l <- vector(mode = "list", length = length(my_k00))

for(i in 1:length(dapc00_l)){
  set.seed(70)
  grp20_l[[i]] <- find.clusters(apos.gc, n.pca = 200, n.clust = my_k00[i])
  dapc00_l[[i]] <- dapc(apos.gc, pop = grp20_l[[i]]$grp, n.pca = 10, n.da = my_k00[i])
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
  dapc01_l[[i]] <- dapc(apos.gc, pop = grp2_l[[i]]$grp, n.pca = 10, n.da = my_k01[i])
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

ggarrange(ggarrange(p01,
                    p02,
                    ncol = 2, labels = c("A", "B")),
          p03,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)



########### DAPC plots

my_k00 <- 9

grp20_l <- vector(mode = "list", length = length(my_k00))
dapc00_l <- vector(mode = "list", length = length(my_k00))

for(i in 1:length(dapc00_l)){
  set.seed(70)
  grp20_l[[i]] <- find.clusters(apos.gc, n.pca = 200, n.clust = my_k00[i])
  dapc00_l[[i]] <- dapc(apos.gc, pop = grp20_l[[i]]$grp, n.pca = 10, n.da = my_k00[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp2_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df00 <- as.data.frame(dapc00_l[[ length(dapc00_l) ]]$ind.coord)
my_df00$Group <- dapc00_l[[ length(dapc00_l) ]]$grp
head(my_df00)

my_df00$Group <- as.character(my_df00$Group)

subset(my_df01, Posterior>0.5 & K==9)
# group A = L06, L16, L17 (laramie big)
# group B = C27, MT big
# group C = L39-A_4, L62-A_5, L16
# group D = L39-A_1,2,3,5
# group E = C43
# group F = S03
# group G = C23
# group H = SM
# group I = C85-1, 2b, c87-5, l45-345
# k9
# 
my_df00[my_df00$Group == 1, "Group"] <- "B"
my_df00[my_df00$Group == 2, "Group"] <- "H"
my_df00[my_df00$Group == 3, "Group"] <- "C"
my_df00[my_df00$Group == 4, "Group"] <- "I"
my_df00[my_df00$Group == 5, "Group"] <- "D"
my_df00[my_df00$Group == 6, "Group"] <- "F"
my_df00[my_df00$Group == 7, "Group"] <- "G"
my_df00[my_df00$Group == 8, "Group"] <- "A"
my_df00[my_df00$Group == 9, "Group"] <- "E"


p04.1 <- ggplot(my_df00, aes(x = LD1, y = LD2, fill = Group))
p04.1 <- p04.1 + geom_point(size = 3, shape = 21)
p04.1 <- p04.1 + theme_bw()
p04.1 <- p04.1 + scale_color_brewer(palette="Paired")
p04.1 <- p04.1 + scale_fill_brewer(palette="Paired")
#p02 <- p02 + scale_color_manual(values=c(col_vector))
#p02 <- p02 + scale_fill_manual(values=c(paste(col_vector, "66", sep = "")))
p04.1



####### K=10 #########

my_k00 <- 10

grp20_l <- vector(mode = "list", length = length(my_k00))
dapc00_l <- vector(mode = "list", length = length(my_k00))

for(i in 1:length(dapc00_l)){
  set.seed(70)
  grp20_l[[i]] <- find.clusters(apos.gc, n.pca = 200, n.clust = my_k00[i])
  dapc00_l[[i]] <- dapc(apos.gc, pop = grp20_l[[i]]$grp, n.pca = 10, n.da = my_k00[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp2_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df00 <- as.data.frame(dapc00_l[[ length(dapc00_l) ]]$ind.coord)
my_df00$Group <- dapc00_l[[ length(dapc00_l) ]]$grp
head(my_df00)

my_df00$Group <- as.character(my_df00$Group)

subset(my_df01, Posterior>0.5 & K==10)
# group A = L06, L16, L17 (laramie big)
# group B = C27, MT big
# group C = L39-A_4, L62-A_5, L16
# group D = L39-A_1,2,3,5
# group E = C43
# group F = S03
# group G = C23
# group H = SM
# group I = MT big 2
# group J = MT big 3
# k10
# 
my_df00[my_df00$Group == 1, "Group"] <- "J"
my_df00[my_df00$Group == 2, "Group"] <- "A"
my_df00[my_df00$Group == 3, "Group"] <- "C"
my_df00[my_df00$Group == 4, "Group"] <- "I"
my_df00[my_df00$Group == 5, "Group"] <- "E"
my_df00[my_df00$Group == 6, "Group"] <- "F"
my_df00[my_df00$Group == 7, "Group"] <- "G"
my_df00[my_df00$Group == 8, "Group"] <- "B"
my_df00[my_df00$Group == 9, "Group"] <- "H"
my_df00[my_df00$Group == 10, "Group"] <- "D"


p04.2 <- ggplot(my_df00, aes(x = LD1, y = LD2, fill = Group))
p04.2 <- p04.2 + geom_point(size = 3, shape = 21)
p04.2 <- p04.2 + theme_bw()
p04.2 <- p04.2 + scale_color_brewer(palette="Paired")
p04.2 <- p04.2 + scale_fill_brewer(palette="Paired")
#p02 <- p02 + scale_color_manual(values=c(col_vector))
#p02 <- p02 + scale_fill_manual(values=c(paste(col_vector, "66", sep = "")))
p04.2


###### K11 ######

my_k00 <- 11

grp20_l <- vector(mode = "list", length = length(my_k00))
dapc00_l <- vector(mode = "list", length = length(my_k00))

for(i in 1:length(dapc00_l)){
  set.seed(70)
  grp20_l[[i]] <- find.clusters(apos.gc, n.pca = 200, n.clust = my_k00[i])
  dapc00_l[[i]] <- dapc(apos.gc, pop = grp20_l[[i]]$grp, n.pca = 10, n.da = my_k00[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp2_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df00 <- as.data.frame(dapc00_l[[ length(dapc00_l) ]]$ind.coord)
my_df00$Group <- dapc00_l[[ length(dapc00_l) ]]$grp
head(my_df00)

my_df00$Group <- as.character(my_df00$Group)

subset(my_df01, Posterior>0.1 & K==11)
# group A = L06, L16-3,5, L17 (laramie big)
# group B = C27, MT big
# group C = L39-A_4, L62-A_5
# group D = L39-A_1,2,3,5
# group E = C43
# group F = S03
# group G = C23
# group H = SM
# group I = MT big 2
# group J = MT big 3
# group K = L16 1,2,4
# k10
# 
my_df00[my_df00$Group == 1, "Group"] <- "G"
my_df00[my_df00$Group == 2, "Group"] <- "K"
my_df00[my_df00$Group == 3, "Group"] <- "E"
my_df00[my_df00$Group == 4, "Group"] <- "A"
my_df00[my_df00$Group == 5, "Group"] <- "H"
my_df00[my_df00$Group == 6, "Group"] <- "J"
my_df00[my_df00$Group == 7, "Group"] <- "F"
my_df00[my_df00$Group == 8, "Group"] <- "C"
my_df00[my_df00$Group == 9, "Group"] <- "D"
my_df00[my_df00$Group == 10, "Group"] <- "B"
my_df00[my_df00$Group == 11, "Group"] <- "I"

p04.3 <- ggplot(my_df00, aes(x = LD1, y = LD2, fill = Group))
p04.3 <- p04.3 + geom_point(size = 3, shape = 21)
p04.3 <- p04.3 + theme_bw()
p04.3 <- p04.3 + scale_color_brewer(palette="Paired")
p04.3 <- p04.3 + scale_fill_brewer(palette="Paired")
#p02 <- p02 + scale_color_manual(values=c(col_vector))
#p02 <- p02 + scale_fill_manual(values=c(paste(col_vector, "66", sep = "")))
p04.3


###### K12 #######


my_k00 <- 12

grp20_l <- vector(mode = "list", length = length(my_k00))
dapc00_l <- vector(mode = "list", length = length(my_k00))

for(i in 1:length(dapc00_l)){
  set.seed(70)
  grp20_l[[i]] <- find.clusters(apos.gc, n.pca = 200, n.clust = my_k00[i])
  dapc00_l[[i]] <- dapc(apos.gc, pop = grp20_l[[i]]$grp, n.pca = 10, n.da = my_k00[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp2_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df00 <- as.data.frame(dapc00_l[[ length(dapc00_l) ]]$ind.coord)
my_df00$Group <- dapc00_l[[ length(dapc00_l) ]]$grp
head(my_df00)

my_df00$Group <- as.character(my_df00$Group)

subset(my_df01, Posterior>0.1 & K==12)
# group A = L06, L16-3,5, L17 (laramie big)
# group B = C27, MT big
# group C = L39-A_4, L62-A_5
# group D = L39-A_1,2,3,5
# group E = C43
# group F = S03
# group G = C23
# group H = SM
# group I = MT big 2
# group J = L16-3,5
# group K = L16 1,2,4
# group L = C43-4
# k12
# 
my_df00[my_df00$Group == 1, "Group"] <- "I"
my_df00[my_df00$Group == 2, "Group"] <- "B"
my_df00[my_df00$Group == 3, "Group"] <- "C"
my_df00[my_df00$Group == 4, "Group"] <- "A"
my_df00[my_df00$Group == 5, "Group"] <- "F"
my_df00[my_df00$Group == 6, "Group"] <- "L"
my_df00[my_df00$Group == 7, "Group"] <- "E"
my_df00[my_df00$Group == 8, "Group"] <- "D"
my_df00[my_df00$Group == 9, "Group"] <- "H"
my_df00[my_df00$Group == 10, "Group"] <- "G"
my_df00[my_df00$Group == 11, "Group"] <- "K"
my_df00[my_df00$Group == 12, "Group"] <- "J"

p04.4 <- ggplot(my_df00, aes(x = LD1, y = LD2, fill = Group))
p04.4 <- p04.4 + geom_point(size = 3, shape = 21)
p04.4 <- p04.4 + theme_bw()
p04.4 <- p04.4 + scale_color_brewer(palette="Paired")
p04.4 <- p04.4 + scale_fill_brewer(palette="Paired")
#p02 <- p02 + scale_color_manual(values=c(col_vector))
#p02 <- p02 + scale_fill_manual(values=c(paste(col_vector, "66", sep = "")))
p04.4


png("kmeans-dapc.png", height=10, width=11, res=300, units="in")
ggarrange(p02,p04.1,p04.2,p04.3,p04.4, ncol=2,nrow=3, labels=c("K=8", "K=9", "K=10", "K=11", "K=12"), hjust = -3, vjust = 2.5, common.legend = TRUE, legend = "none")
dev.off()
