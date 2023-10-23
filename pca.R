#########
# PCA ###
#########

# replace NA's
sum(is.na(AllPops.gc$tab))
AllPops.gc.scaled <- scaleGen(AllPops.gc, NA.method="mean")

hookeri.pca.all <- dudi.pca(AllPops.gc.scaled, cent=FALSE, scale=FALSE, scannf=FALSE, nf=8)

pca.all.df <- as.data.frame(hookeri.pca.all$li)
pca.all.df$pop <- pop(AllPops.gc)
pca.all.df$ms <- AllPops.gc$strata$ms

pca.all.df$group <- "blank"
pca.all.df$group[1:18] <- "CO-S"
pca.all.df$group[7:9] <- "Laramie-S"
pca.all.df$group[16:18] <- "YK-S"
pca.all.df$group[19:41] <- "Laramie-S"
pca.all.df$group[37] <- "L45-S"
pca.all.df$group[42:46] <- "ND-A"
pca.all.df$group[47:51] <- "MT.big-A"
pca.all.df$group[52:56] <- "BC-A"
pca.all.df$group[57:61] <- "MT.big-A"
pca.all.df$group[62:76] <- "MT.big-A"
pca.all.df$group[77:81] <- "Laramie.big-A"
pca.all.df$group[c(82,83,85)] <- "Laramie.small-A"
pca.all.df$group[c(84,86,87,88,89,90,91)] <- "Laramie.big-A"
pca.all.df$group[c(92,93,94,96)] <- "L39.big-A"
pca.all.df$group[95] <- "L39.small-A"
pca.all.df$group[97:100] <- "MT.big-A"
pca.all.df$group[101:104] <- "MT.big-A"
pca.all.df$group[105] <- "L62-A"
pca.all.df$group[106:109] <- "SK-A"
pca.all.df$group[110:114] <- "YK-A"

pca.all.df$group <- factor(pca.all.df$group, levels=c("CO-S", "Laramie-S", "L62-A", "Laramie.big-A", "Laramie.small-A", "L45-S", "L39.big-A", "L39.small-A", "MT.big-A", "ND-A", "SK-A", "BC-A", "YK-A", "YK-S"))

pca.gg.12 <- ggplot(pca.all.df, aes(x = Axis1, y = Axis2))+ 
  geom_point(aes(fill=group, shape=ms), size=3, stroke=0.7)+ 
  theme_bw()+
  scale_shape_manual(values = c(24,21))+
  scale_fill_manual(values=group.cols)+
  #geom_text(aes(label=pop), size=3)+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 24, 21, 21, 21, 24, 21, 21, 21, 21, 21, 21, 21, 24))))

pca.gg.13 <- ggplot(pca.all.df, aes(x = Axis1, y = Axis3))+ 
  geom_point(aes(fill=group, shape=ms), size=3, stroke=0.7)+ 
  theme_bw()+
  scale_shape_manual(values = c(24,21))+
  scale_fill_manual(values=group.cols)+
  #geom_text(aes(label=pop), size=3)+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 24, 21, 21, 21, 24, 21, 21, 21, 21, 21, 21, 21, 24))))

pca.gg.23 <- ggplot(pca.all.df, aes(x = Axis2, y = Axis3))+ 
  geom_point(aes(fill=group, shape=ms), size=3, stroke=0.7)+ 
  theme_bw()+
  scale_shape_manual(values = c(24,21))+
  scale_fill_manual(values=group.cols)+
  #geom_text(aes(label=pop), size=3)+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 24, 21, 21, 21, 24, 21, 21, 21, 21, 21, 21, 21, 24))))
