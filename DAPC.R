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


##########################
###### for pops (all)#####
##########################

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
  geom_point(aes(fill=group, shape=ms), size=2, stroke=0.7)+ 
  theme_bw()+
  scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values=DAPC.cols)+
  #geom_text(aes(label=pop), size=3)+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 24, 21, 21, 21, 24, 21, 21, 21, 21, 21, 21, 21, 24))))

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



# sub a few pops
noYK.gc <- popsub(AllPops.gc, blacklist=c("C59-S", "SM-A", "C23-A", "S03-A"))
hookeri.dapc2 <- dapc(noYK.gc, grp=noYK.gc$pop, n.pca=20, n.da=100)
scatter(hookeri.dapc2, grp = noYK.gc$pop, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch=my.pch.sub)

# apo only pops
setPop(apos.gc) <- ~pop
apos.gc$pop <- factor(apos.gc$pop, levels=c("L62-A", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A"))
apo.dapc <- dapc(apos.gc, grp=apos.gc$grp, n.pca=10, n.da=100)
scatter(apo.dapc, grp = apos.gc$pop, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)

# sex only pops
setPop(dips.gc) <- ~pop
dips.gc$pop <- factor(dips.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S","L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L45-S", "C59-S"))
dips.dapc <- dapc(dips.gc, grp=dips.gc$grp)
scatter(dips.dapc, grp = dips.gc$pop, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)

# OneMLG
pch.OneMLG <- c(17, 17, 17, 17, 17, 17, 21, 17, 17, 17, 17, 17, 17, 21, 21, 21, 21, 17, 21, 21, 21, 21, 17)
pch.OneMLG <-replace(my.pch,my.pch==21, 19)
setPop(OneMLG.gc) <- ~pop
OneMLG.gc$pop <- factor(OneMLG.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))
OneMLG.dapc <- dapc(OneMLG.gc, grp=OneMLG.gc$grp)
scatter(OneMLG.dapc, grp = OneMLG.gc$pop, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch=pch.OneMLG)
