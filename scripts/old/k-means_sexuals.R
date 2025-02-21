# k-means sexuals

# All pops

# grp3 <- find.clusters(dips.gc, max.n.clust=30) # looks like lowest BIC is 13 clusters...
kmeans.dips <- find.clusters(dips.gc)
dips.tab <- tab(dips.gc, NA.method="mean")
dips.xval <- xvalDapc(dips.tab, grp=kmeans.dips$grp)
dips.xval[2:6]

dips.dapc <- dapc(dips.gc, pop=kmeans.dips$grp)
# names(grp3)
# table(pop(dips.gc), grp3$grp3)
# table.value(table(pop(dips.gc), grp3$grp3), col.lab=paste("inf", 1:13),
#             row.lab=paste("ori", 1:13))

## from vcfR -- run k-means multiple times
# library(vcfR)
# vcf <- read.vcfR("prubi_gbs.vcf.gz")
# pop.data <- read.table("population_data.gbs.txt", sep = "\t", header = TRUE)
# all(colnames(vcf@gt)[-1] == pop.data$AccessID)
# ## [1] TRUE
# gl_rubi <- vcfR2genlight(vcf)

# library(adegenet)
maxK3 <- 10
myMat3 <- matrix(nrow=50, ncol=maxK3)
colnames(myMat3) <- 1:ncol(myMat3)
for(i in 1:nrow(myMat3)){
  grp3 <- find.clusters(dips.gc, n.pca = 200, choose.n.clust = FALSE,  max.n.clust = maxK3)
  myMat3[i,] <- grp3$Kstat
}

my_df3 <- melt(myMat3)
colnames(my_df3)[1:3] <- c("Group", "K", "BIC")
my_df3$K <- as.factor(my_df3$K)
head(my_df3)

p001 <- ggplot(my_df3, aes(x = K, y = BIC))
p001 <- p001 + geom_boxplot()
p001 <- p001 + theme_bw()
p001 <- p001 + xlab("Number of groups (K)")
p001 # X clusters kind of looks right

# plot 2 #
my_k0000 <- 4

grp30_l <- vector(mode = "list", length = length(my_k0000))
dapc00_l <- vector(mode = "list", length = length(my_k0000))

for(i in 1:length(dapc00_l)){
  set.seed(10)
  grp30_l[[i]] <- find.clusters(dips.gc, n.pca = 200, n.clust = my_k0000[i])
  dapc00_l[[i]] <- dapc(dips.gc, pop = grp30_l[[i]]$grp, n.pca = 5, n.da = my_k0000[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp3_l[[i]]$grp3, n.pca = 3, n.da = 2)
}

my_df00 <- as.data.frame(dapc00_l[[ length(dapc00_l) ]]$ind.coord)
my_df00$Group <- dapc00_l[[ length(dapc00_l) ]]$grp
head(my_df00)

p002 <- ggplot(my_df00, aes(x = LD1, y = LD2, fill = Group))
p002 <- p002 + geom_point(size = 4, shape = 21, color="black")
p002 <- p002 + theme_bw()
#p002 <- p002 + scale_color_manual(values=c(col_vector))
p002 <- p002 + scale_fill_manual(values=cols)
p002

# plot 3
my_k00 <- 2:4

grp003_l <- vector(mode = "list", length = length(my_k00))
dapc00_l <- vector(mode = "list", length = length(my_k00))

for(i in 1:length(dapc_l)){
  set.seed(10)
  grp003_l[[i]] <- find.clusters(dips.gc, n.pca = 200, n.clust = my_k00[i])
  dapc00_l[[i]] <- dapc(dips.gc, pop = grp003_l[[i]]$grp, n.pca = 5, n.da = my_k00[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp3_l[[i]]$grp3, n.pca = 3, n.da = 2)
}

my_df000 <- as.data.frame(dapc00_l[[ length(dapc00_l) ]]$ind.coord)
my_df000$Group <- dapc00_l[[ length(dapc00_l) ]]$grp
head(my_df000)

tmp00 <- as.data.frame(dapc00_l[[1]]$posterior)
tmp00$K <- my_k00[1]
tmp00$Sample <- rownames(tmp00)
tmp00 <- melt(tmp00, id = c("Sample", "K"))
names(tmp00)[3:4] <- c("Group", "Posterior")
tmp00$pop <- mystrata$pop[match(tmp00$Sample, mystrata$id)]
my_df000 <- tmp00

for(i in 2:length(dapc00_l)){
  tmp00 <- as.data.frame(dapc00_l[[i]]$posterior)
  tmp00$K <- my_k00[i]
  tmp00$Sample <- rownames(tmp00)
  tmp00 <- melt(tmp00, id = c("Sample", "K"))
  names(tmp00)[3:4] <- c("Group", "Posterior")
  tmp00$pop <- mystrata$pop[match(tmp00$Sample, mystrata$id)]
  my_df000 <- rbind(my_df000, tmp00)
}

grp3.labs <- paste("K =", my_k00)
names(grp3.labs) <- my_k00

my_df000$pop <- factor(my_df000$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S","L45-S", "C59-S"))

# my_df00$Group <- as.character(my_df00$Group)
my_df000$Group <- as.character(my_df000$Group)

#black  #lightgrey #brown #darkgrey
my_df000[ my_df000$K == 2 & my_df000$Group == 1, "Group"] <- "C"
my_df000[ my_df000$K == 2 & my_df000$Group == 2, "Group"] <- "A"

#black  #lightgrey #brown #darkgrey
my_df000[ my_df000$K == 3 & my_df000$Group == 1, "Group"] <- "A"
my_df000[ my_df000$K == 3 & my_df000$Group == 2, "Group"] <- "B"
my_df000[ my_df000$K == 3 & my_df000$Group == 3, "Group"] <- "C"

#black  #lightgrey #brown #darkgrey
my_df000[ my_df000$K == 4 & my_df000$Group == 1, "Group"] <- "D"
my_df000[ my_df000$K == 4 & my_df000$Group == 2, "Group"] <- "C"
my_df000[ my_df000$K == 4 & my_df000$Group == 3, "Group"] <- "B"
my_df000[ my_df000$K == 4 & my_df000$Group == 4, "Group"] <- "A"


p003 <- ggplot(my_df000, aes(x = Sample, y = Posterior, fill = Group))
p003 <- p003 + geom_bar(stat = "identity")
p003 <- p003 + facet_grid(K ~ pop, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp3.labs))
p003 <- p003 + theme_bw()
p003 <- p003 + ylab("Posterior membership probability")
p003 <- p003 + theme(legend.position='none')
#p003 <- p003 + scale_color_brewer(palette="Paired")
#black  #lightgrey #peach #darkgrey
p003 <- p003 + scale_fill_manual(values=c("black", "azure3", "#FDBF6F", "azure4"))
p003 <- p003 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),panel.spacing.x=unit(0.1, "lines"))
p003

ggarrange(ggarrange(p001,
                    p2,
                    ncol = 2, labels = c("A", "B")),
          p3,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)
