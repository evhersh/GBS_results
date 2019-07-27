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
p01 # 12 clusters kind of looks right?

# plot 2 #
my_k00 <- 8

grp20_l <- vector(mode = "list", length = length(my_k00))
dapc00_l <- vector(mode = "list", length = length(my_k00))

for(i in 1:length(dapc00_l)){
  set.seed(10)
  grp20_l[[i]] <- find.clusters(apos.gc, n.pca = 200, n.clust = my_k00[i])
  dapc00_l[[i]] <- dapc(apos.gc, pop = grp20_l[[i]]$grp, n.pca = 8, n.da = my_k00[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp2_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df00 <- as.data.frame(dapc00_l[[ length(dapc00_l) ]]$ind.coord)
my_df00$Group <- dapc00_l[[ length(dapc00_l) ]]$grp
head(my_df00)

p02 <- ggplot(my_df00, aes(x = LD1, y = LD2, color = Group, fill = Group))
p02 <- p02 + geom_point(size = 4, shape = 21)
p02 <- p02 + theme_bw()
p02 <- p02 + scale_color_brewer(palette="Paired")
p02 <- p02 + scale_fill_brewer(palette="Paired")
#p02 <- p02 + scale_color_manual(values=c(col_vector))
#p02 <- p02 + scale_fill_manual(values=c(paste(col_vector, "66", sep = "")))
#p02

# plot 3
my_k01 <- 8:10

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

#my_df001 <- my_df01
my_df01 <- my_df001

my_df01$Group <- as.character(my_df01$Group)
my_df01$Group <- as.character(my_df01$Group)

# G-B
my_df01[ my_df01$K == 8 & my_df01$Group == 1, "Group"] <- "A"
my_df01[ my_df01$K == 8 & my_df01$Group == 2, "Group"] <- "G"
my_df01[ my_df01$K == 8 & my_df01$Group == 3, "Group"] <- "C"
my_df01[ my_df01$K == 8 & my_df01$Group == 4, "Group"] <- "D"
my_df01[ my_df01$K == 8 & my_df01$Group == 5, "Group"] <- "E"
my_df01[ my_df01$K == 8 & my_df01$Group == 6, "Group"] <- "F"
my_df01[ my_df01$K == 8 & my_df01$Group == 7, "Group"] <- "B"
my_df01[ my_df01$K == 8 & my_df01$Group == 8, "Group"] <- "H"

# G B A F C D I H E
# A-C, C-I, C-E, I-G, G-B
my_df01[ my_df01$K == 9 & my_df01$Group == 1, "Group"] <- "B"
my_df01[ my_df01$K == 9 & my_df01$Group == 2, "Group"] <- "G"
my_df01[ my_df01$K == 9 & my_df01$Group == 3, "Group"] <- "A"
my_df01[ my_df01$K == 9 & my_df01$Group == 4, "Group"] <- "F"
my_df01[ my_df01$K == 9 & my_df01$Group == 5, "Group"] <- "C"
my_df01[ my_df01$K == 9 & my_df01$Group == 6, "Group"] <- "D"
my_df01[ my_df01$K == 9 & my_df01$Group == 7, "Group"] <- "I"
my_df01[ my_df01$K == 9 & my_df01$Group == 8, "Group"] <- "H"
my_df01[ my_df01$K == 9 & my_df01$Group == 9, "Group"] <- "E"

# J A G B E C D I F H
# B-A, D-G, G-B, J-H, G-C, I-J, C-J, C-F, G-B
my_df01[ my_df01$K == 10 & my_df01$Group == 1, "Group"] <- "J"
my_df01[ my_df01$K == 10 & my_df01$Group == 2, "Group"] <- "A"
my_df01[ my_df01$K == 10 & my_df01$Group == 3, "Group"] <- "B"
my_df01[ my_df01$K == 10 & my_df01$Group == 4, "Group"] <- "G"
my_df01[ my_df01$K == 10 & my_df01$Group == 5, "Group"] <- "E"
my_df01[ my_df01$K == 10 & my_df01$Group == 6, "Group"] <- "C"
my_df01[ my_df01$K == 10 & my_df01$Group == 7, "Group"] <- "D"
my_df01[ my_df01$K == 10 & my_df01$Group == 8, "Group"] <- "I"
my_df01[ my_df01$K == 10 & my_df01$Group == 9, "Group"] <- "F"
my_df01[ my_df01$K == 10 & my_df01$Group == 10, "Group"] <- "H"

# F B C D E A G H I J K
# F-A, F-H, C-F, C-I, K-G, B-I, D-J, E-I, I-B
# my_df01[ my_df01$K == 11 & my_df01$Group == 1, "Group"] <- "H"
# my_df01[ my_df01$K == 11 & my_df01$Group == 2, "Group"] <- "E"
# my_df01[ my_df01$K == 11 & my_df01$Group == 3, "Group"] <- "F"
# my_df01[ my_df01$K == 11 & my_df01$Group == 4, "Group"] <- "J"
# my_df01[ my_df01$K == 11 & my_df01$Group == 5, "Group"] <- "B"
# my_df01[ my_df01$K == 11 & my_df01$Group == 6, "Group"] <- "A"
# my_df01[ my_df01$K == 11 & my_df01$Group == 7, "Group"] <- "K"
# my_df01[ my_df01$K == 11 & my_df01$Group == 8, "Group"] <- "I"
# my_df01[ my_df01$K == 11 & my_df01$Group == 9, "Group"] <- "C"
# my_df01[ my_df01$K == 11 & my_df01$Group == 10, "Group"] <- "D"
# my_df01[ my_df01$K == 11 & my_df01$Group == 11, "Group"] <- "G"
# 
# 
# 
# my_df01[ my_df01$K == 12 & my_df01$Group == 1, "Group"] <- "A"
# my_df01[ my_df01$K == 12 & my_df01$Group == 2, "Group"] <- "B"
# my_df01[ my_df01$K == 12 & my_df01$Group == 3, "Group"] <- "C"
# my_df01[ my_df01$K == 12 & my_df01$Group == 4, "Group"] <- "D"
# my_df01[ my_df01$K == 12 & my_df01$Group == 5, "Group"] <- "E"
# my_df01[ my_df01$K == 12 & my_df01$Group == 6, "Group"] <- "F"
# my_df01[ my_df01$K == 12 & my_df01$Group == 7, "Group"] <- "G"
# my_df01[ my_df01$K == 12 & my_df01$Group == 8, "Group"] <- "H"
# my_df01[ my_df01$K == 12 & my_df01$Group == 9, "Group"] <- "I"
# my_df01[ my_df01$K == 12 & my_df01$Group == 10, "Group"] <- "J"
# my_df01[ my_df01$K == 12 & my_df01$Group == 11, "Group"] <- "K"
# my_df01[ my_df01$K == 12 & my_df01$Group == 12, "Group"] <- "L"



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




