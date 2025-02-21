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
library(ggpattern)
library(tidyverse)

##### KMEANS #####

kmeans.all <- find.clusters(AllPops.gc)
all.tab <- tab(AllPops.gc, NA.method="mean")
all.xval <- xvalDapc(all.tab, grp=kmeans.all$grp)
all.xval[2:6]

dapc.all <- dapc(AllPops.gc, n.pca=50, n.da=20)

optim.a.score(dapc.all)

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

##### plot 1 #####

p1 <- ggplot(my_df2, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1 # 4-6 clusters kind of looks right

##### plot 2 #####
my_k0 <- 4

grp0_l <- vector(mode = "list", length = length(my_k0))
dapc0_l <- vector(mode = "list", length = length(my_k0))

for(i in 1:length(dapc0_l)){
  set.seed(10)
  grp0_l[[i]] <- find.clusters(AllPops.gc, n.pca = 200, n.clust = my_k0[i])
  dapc0_l[[i]] <- dapc(AllPops.gc, pop = grp0_l[[i]]$grp, n.pca = 10, n.da = my_k0[i])
  #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
}

my_df0 <- as.data.frame(dapc0_l[[ length(dapc0_l) ]]$ind.coord)
my_df0$Group <- dapc0_l[[ length(dapc0_l) ]]$grp
head(my_df0)

p2 <- ggplot(my_df0, aes(x = LD1, y = LD2, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21, color="black")
p2 <- p2 + theme_bw()
#p2 <- p2 + scale_color_manual(values=c(col_vector))
p2 <- p2 + scale_fill_manual(values=cols)
p2

##### plot 3 #####
my_k <- 4:6

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(10)
  grp_l[[i]] <- find.clusters(AllPops.gc, n.pca = 200, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(AllPops.gc, pop = grp_l[[i]]$grp, n.pca = 10, n.da = my_k[i])
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
#tmp$pop <- AllPops.gc@other$group[match(tmp$Sample, indNames(AllPops.gc))]
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

my_df$pop <- factor(my_df$pop, levels=c("CO1-S", "CO2-S", "CO3-S", 
                                        "CO4-S", "CO5-S", "WY1-S", 
                                        "WY1-A", "WY2-S", "WY3-S", 
                                        "WY4-S", "WY6-S", "WY7-S", 
                                        "WY8-S", "WY9-A", "WY10-A", 
                                        "WY11-A", "WY12-A", "WY13-A",
                                        "WY14-S", "WY14-A", "WY15-A", 
                                        "MT1-A", "MT2-A", "MT3-A", 
                                        "MT4-A", "ND-A", "BC-A", 
                                        "SK-A", "YK1-A", "YK2-S"))

# my_df$pop <- factor(my_df$pop, levels=c(
#   "CO-S",
#   "Laramie-S",
#   "L62-A",
#   "Laramie.big-A",
#   "Laramie.small-A",
#   "L45-S",
#   "L39.big-A",
#   "L39.small-A",
#   "MT.big-A",
#   "ND-A",
#   "SK-A",
#   "BC-A",
#   "YK-A",
#   "YK-S"
# ))

my_df <- my_df %>%
  mutate(ms = case_when(str_detect(pop, "-S") ~ "S", str_detect(pop, "-A") ~ "A"))

# my_df0$Group <- as.character(my_df0$Group)
my_df$Group <- as.character(my_df$Group)
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
# #black  #grey   #blue       #brown    # yellow    #orange
my_df[ my_df$K == 4 & my_df$Group == 1, "Group"] <- "D"
my_df[ my_df$K == 4 & my_df$Group == 2, "Group"] <- "A"
my_df[ my_df$K == 4 & my_df$Group == 3, "Group"] <- "C"
my_df[ my_df$K == 4 & my_df$Group == 4, "Group"] <- "B"
# my_df[ my_df$K == 10 & my_df$Group == 5, "Group"] <- "E"
# my_df[ my_df$K == 10 & my_df$Group == 6, "Group"] <- "F"
# my_df[ my_df$K == 10 & my_df$Group == 7, "Group"] <- "G"
# my_df[ my_df$K == 10 & my_df$Group == 8, "Group"] <- "H"
# my_df[ my_df$K == 10 & my_df$Group == 9, "Group"] <- "I"
# my_df[ my_df$K == 10 & my_df$Group == 10, "Group"] <- "J"
# #black  #grey   #blue       #brown    # yellow    #orange
my_df[ my_df$K == 5 & my_df$Group == 1, "Group"] <- "C"
my_df[ my_df$K == 5 & my_df$Group == 2, "Group"] <- "D"
my_df[ my_df$K == 5 & my_df$Group == 3, "Group"] <- "A"
my_df[ my_df$K == 5 & my_df$Group == 4, "Group"] <- "B"
my_df[ my_df$K == 5 & my_df$Group == 5, "Group"] <- "E"
# #black  #grey   #blue       #brown    # yellow    #orange
my_df[ my_df$K == 6 & my_df$Group == 1, "Group"] <- "A"
my_df[ my_df$K == 6 & my_df$Group == 2, "Group"] <- "C"
my_df[ my_df$K == 6 & my_df$Group == 3, "Group"] <- "B"
my_df[ my_df$K == 6 & my_df$Group == 4, "Group"] <- "F"
my_df[ my_df$K == 6 & my_df$Group == 5, "Group"] <- "D"
my_df[ my_df$K == 6 & my_df$Group == 6, "Group"] <- "E"

p3 <- ggplot(my_df, aes(x = Sample, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ pop, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
#p3 <- p3 + scale_color_brewer(palette="Paired")
                                        #black  #grey   #blue       #peach    # yellow    #orange
# p3 <- p3 + scale_fill_manual(values=c("#37474f",
#                                       "#607d8b",
#                                       "#79554b",
#                                       "#4e342e",
#                                       "#cfd8dc",
#                                       "#bcaaa4"))

p3 <- p3 + scale_fill_manual(values=c("#4e342e",
                                      "#79554b",
                                      "#607d8b",
                                      "#37474f",
                                      "#bcaaa4",
                                      "#cfd8dc"))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),panel.spacing.x=unit(0.1, "lines"),
                 strip.text.x.top = element_text(angle=90))
p3


png("kmeans-all.png", height=3, width=11, res=300, units="in")
p3
dev.off()

# ggarrange(ggarrange(p1,
#                     p2,
#                     ncol = 2, labels = c("A", "B")),
#           p3,
#           nrow = 2,
#           labels = c("", "C"),
#           heights = c(1, 2)
# )



# Apos only
# Ensure enough colors for the unique groups
unique_groups <- length(unique(my_df$Group))  # Count unique values in Group
colors <- RColorBrewer::brewer.pal(max(3, unique_groups), "Dark2")  # Pick a palette


# First Plot (p3)
p3 <- ggplot(my_df, aes(x = Sample, y = Posterior, fill = Group)) +
  geom_bar(stat = "identity", color = "white") +
  facet_grid(K ~ pop, scales = "free_x", space = "free", labeller = labeller(K = grp.labs)) +
  theme_bw() +
  ylab("Posterior membership probability") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    panel.spacing.x = unit(0.1, "lines"),
    strip.text.x.top = element_text(angle = 90)
  ) +
  scale_fill_viridis(discrete = TRUE)  # Use viridis color scale for fill

p3

# Combine with Hobs plot
theme_adjustment <- theme(
  plot.margin = margin(0, 0, 0, 0),  # Adjust margins for both plots
  panel.spacing = unit(0, "lines")    # No space between panels
)

p3 <- p3 + theme_adjustment
gg.Hobs.bar <- gg.Hobs.bar + theme_adjustment

# Stack with minimal white space
kh.plot <- gg.Hobs.bar / plot_spacer() / p3 + plot_layout(ncol = 1, heights = c(1,-0.11, 1), guides="collect")

png("./figures/kmeans_hobs.png", height=7, width=10, res=300, units="in")
kh.plot
dev.off()
