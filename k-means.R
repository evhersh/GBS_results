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

# All pops

grp <- find.clusters(AllPops.gc, max.n.clust=30) # looks like lowest BIC is 13 clusters...
names(grp)
table(pop(AllPops.gc), grp$grp)
table.value(table(pop(AllPops.gc), grp$grp), col.lab=paste("inf", 1:13),
            row.lab=paste("ori", 1:13))

## from vcfR -- run k-means multiple times
# library(vcfR)
# vcf <- read.vcfR("prubi_gbs.vcf.gz")
# pop.data <- read.table("population_data.gbs.txt", sep = "\t", header = TRUE)
# all(colnames(vcf@gt)[-1] == pop.data$AccessID)
# ## [1] TRUE
# gl_rubi <- vcfR2genlight(vcf)

# library(adegenet)
maxK <- 30
myMat <- matrix(nrow=50, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(AllPops.gc, n.pca = 200, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

my_df2 <- melt(myMat)
colnames(my_df2)[1:3] <- c("Group", "K", "BIC")
my_df2$K <- as.factor(my_df2$K)
head(my_df2)

p1 <- ggplot(my_df2, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1 # 10 clusters kind of looks right

# plot 2 #
my_k0 <- 9

grp0_l <- vector(mode = "list", length = length(my_k0))
dapc0_l <- vector(mode = "list", length = length(my_k0))

for(i in 1:length(dapc0_l)){
  set.seed(10)
  grp0_l[[i]] <- find.clusters(AllPops.gc, n.pca = 200, n.clust = my_k0[i])
  dapc0_l[[i]] <- dapc(AllPops.gc, pop = grp0_l[[i]]$grp, n.pca = 20, n.da = my_k0[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df0 <- as.data.frame(dapc0_l[[ length(dapc0_l) ]]$ind.coord)
my_df0$Group <- dapc0_l[[ length(dapc0_l) ]]$grp
head(my_df0)

p2 <- ggplot(my_df0, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(col_vector))
p2 <- p2 + scale_fill_manual(values=c(paste(col_vector, "66", sep = "")))
p2

# plot 3
my_k <- 9:12

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(10)
  grp_l[[i]] <- find.clusters(AllPops.gc, n.pca = 200, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(AllPops.gc, pop = grp_l[[i]]$grp, n.pca = 20, n.da = my_k[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Sample <- rownames(tmp)
tmp <- melt(tmp, id = c("Sample", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$pop <- mystrata$pop[match(tmp$Sample, mystrata$id)]
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Sample <- rownames(tmp)
  tmp <- melt(tmp, id = c("Sample", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$pop <- mystrata$pop[match(tmp$Sample, mystrata$id)]
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

my_df$pop <- factor(my_df$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))

# my_df0$Group <- as.character(my_df0$Group)
# my_df$Group <- as.character(my_df$Group)
# 
# my_df0[ my_df0$K == 10 & my_df0$Group == 1, "Group"] <- "A"
# my_df0[ my_df0$K == 10 & my_df0$Group == 2, "Group"] <- "B"
# my_df0[ my_df0$K == 10 & my_df0$Group == 3, "Group"] <- "C"
# my_df0[ my_df0$K == 10 & my_df0$Group == 4, "Group"] <- "D"
# my_df0[ my_df0$K == 10 & my_df0$Group == 5, "Group"] <- "E"
# my_df0[ my_df0$K == 10 & my_df0$Group == 6, "Group"] <- "F"
# my_df0[ my_df0$K == 10 & my_df0$Group == 7, "Group"] <- "G"
# my_df0[ my_df0$K == 10 & my_df0$Group == 8, "Group"] <- "H"
# my_df0[ my_df0$K == 10 & my_df0$Group == 9, "Group"] <- "I"
# my_df0[ my_df0$K == 10 & my_df0$Group == 10, "Group"] <- "J"
# 
# my_df[ my_df$K == 10 & my_df$Group == 1, "Group"] <- "A"
# my_df[ my_df$K == 10 & my_df$Group == 2, "Group"] <- "B"
# my_df[ my_df$K == 10 & my_df$Group == 3, "Group"] <- "C"
# my_df[ my_df$K == 10 & my_df$Group == 4, "Group"] <- "D"
# my_df[ my_df$K == 10 & my_df$Group == 5, "Group"] <- "E"
# my_df[ my_df$K == 10 & my_df$Group == 6, "Group"] <- "F"
# my_df[ my_df$K == 10 & my_df$Group == 7, "Group"] <- "G"
# my_df[ my_df$K == 10 & my_df$Group == 8, "Group"] <- "H"
# my_df[ my_df$K == 10 & my_df$Group == 9, "Group"] <- "I"
# my_df[ my_df$K == 10 & my_df$Group == 10, "Group"] <- "J"
# 
# my_df[ my_df$K == 11 & my_df$Group == 1, "Group"] <- "A"
# my_df[ my_df$K == 11 & my_df$Group == 2, "Group"] <- "B"
# my_df[ my_df$K == 11 & my_df$Group == 3, "Group"] <- "C"
# my_df[ my_df$K == 11 & my_df$Group == 4, "Group"] <- "D"
# my_df[ my_df$K == 11 & my_df$Group == 5, "Group"] <- "E"
# my_df[ my_df$K == 11 & my_df$Group == 6, "Group"] <- "F"
# my_df[ my_df$K == 11 & my_df$Group == 7, "Group"] <- "G"
# my_df[ my_df$K == 11 & my_df$Group == 8, "Group"] <- "H"
# my_df[ my_df$K == 11 & my_df$Group == 9, "Group"] <- "I"
# my_df[ my_df$K == 11 & my_df$Group == 10, "Group"] <- "J"
# my_df[ my_df$K == 11 & my_df$Group == 11, "Group"] <- "K"
# 
# my_df[ my_df$K == 12 & my_df$Group == 1, "Group"] <- "A"
# my_df[ my_df$K == 12 & my_df$Group == 2, "Group"] <- "B"
# my_df[ my_df$K == 12 & my_df$Group == 3, "Group"] <- "C"
# my_df[ my_df$K == 12 & my_df$Group == 4, "Group"] <- "D"
# my_df[ my_df$K == 12 & my_df$Group == 5, "Group"] <- "E"
# my_df[ my_df$K == 12 & my_df$Group == 6, "Group"] <- "F"
# my_df[ my_df$K == 12 & my_df$Group == 7, "Group"] <- "G"
# my_df[ my_df$K == 12 & my_df$Group == 8, "Group"] <- "H"
# my_df[ my_df$K == 12 & my_df$Group == 9, "Group"] <- "I"
# my_df[ my_df$K == 12 & my_df$Group == 10, "Group"] <- "J"
# my_df[ my_df$K == 12 & my_df$Group == 11, "Group"] <- "K"
# my_df[ my_df$K == 12 & my_df$Group == 12, "Group"] <- "L"
# 
# my_df[ my_df$K == 13 & my_df$Group == 1, "Group"] <- "A"
# my_df[ my_df$K == 13 & my_df$Group == 2, "Group"] <- "B"
# my_df[ my_df$K == 13 & my_df$Group == 3, "Group"] <- "C"
# my_df[ my_df$K == 13 & my_df$Group == 4, "Group"] <- "D"
# my_df[ my_df$K == 13 & my_df$Group == 5, "Group"] <- "E"
# my_df[ my_df$K == 13 & my_df$Group == 6, "Group"] <- "F"
# my_df[ my_df$K == 13 & my_df$Group == 7, "Group"] <- "G"
# my_df[ my_df$K == 13 & my_df$Group == 8, "Group"] <- "H"
# my_df[ my_df$K == 13 & my_df$Group == 9, "Group"] <- "I"
# my_df[ my_df$K == 13 & my_df$Group == 10, "Group"] <- "J"
# my_df[ my_df$K == 13 & my_df$Group == 11, "Group"] <- "K"
# my_df[ my_df$K == 13 & my_df$Group == 12, "Group"] <- "L"
# my_df[ my_df$K == 13 & my_df$Group == 13, "Group"] <- "M"


p3 <- ggplot(my_df, aes(x = Sample, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ pop, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
#p3 <- p3 + scale_color_brewer(palette="Paired")
p3 <- p3 + scale_fill_manual(values=c(col_vector))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),panel.spacing.x=unit(0.1, "lines"))
p3

ggarrange(ggarrange(p1,
                    p2,
                    ncol = 2, labels = c("A", "B")),
          p3,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)


# Apos only

