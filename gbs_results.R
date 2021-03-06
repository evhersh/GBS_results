## @knitr load

# packages
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
library(rmarkdown)
library(ggpubr)


# Data import 

# strata (sample, pop, ms)
mystrata <- read.csv("~/Google Drive/GitHub/Hookeri-GBS/popmap_all.csv")

# vcf to vcfR
# vcf <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/final.filtered.snps.vcf")
# vcf.dips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.dips.vcf")
# vcf.trips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.trips.vcf")
# vcf.tets <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.tets.vcf")
# 
# # vcfR to genind to genclone
# dips.gi <- vcfR2genind(vcf.dips, sep = "/", ploidy=2)
# dips.gc <- as.genclone(dips.gi)
# sampleorder <- match(indNames(dips.gc), mystrata$id)
# strata(dips.gc) <- mystrata[sampleorder,]
# setPop(dips.gc) <- ~pop
# 
# trips.gi <- vcfR2genind(vcf.trips, sep = "/", ploidy=3)
# trips.gc <- as.genclone(trips.gi)
# sampleorder <- match(indNames(trips.gc), mystrata$id)
# strata(trips.gc) <- mystrata[sampleorder,]
# setPop(trips.gc) <- ~pop
# 
# tets.gi <- vcfR2genind(vcf.tets, sep = "/", ploidy=4)
# tets.gc <- as.genclone(tets.gi)
# sampleorder <- match(indNames(tets.gc), mystrata$id)
# strata(tets.gc) <- mystrata[sampleorder,]
# setPop(tets.gc) <- ~pop
# 
# # combine genclones
# dipsNtripsNtets.gc <- repool(dips.gc,trips.gc,tets.gc)
# #dipsNtripsNtets.gc$pop <- factor(dipsNtripsNtets.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))
# 
# AllPops.gc <-as.genclone(dipsNtripsNtets.gc)
# AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))

# just re-save in 'poppr.R' if I re-call the SNPs or anything...
load("AllPops.gc.RData")
mll(AllPops.gc)

set.seed(420)

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))

#cols <- brewer.pal(n = nPop(AllPops.gc), name = "Dark2")



## @knitr MLGs
load("AllPops.gc.RData")
apo.list <- c("C23-A", "C27-A", "C43-A", "C85-A", "C86-A", "C87-A", "C88-A", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A", "L45-A", "L62-A", "S03-A", "SM-A")
sex.list <- c("B42-S", "B46-S", "B49-S", "B53-S", "B60-S", "C59-S", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L45-S", "L62-S")
# calculate raw euclidian distance
ddist2 <- prevosti.dist(AllPops.gc)

# assign MLG's using raw euclidian distance from dist() [above]
fstats <- filter_stats(AllPops.gc, distance=ddist2, plot=TRUE, stats="THRESHOLD")
cutoff_predictor(fstats$farthest)
# looks like this gives the same clone mlg assignments as my IBS stuff
mlg.filter(AllPops.gc, distance=ddist2) <- 0.1

mlgtab <- mlg.table(AllPops.gc)
#table.value(mlgtab)

apo.table <- mlg.table(AllPops.gc, sublist=apo.list, plot=FALSE)

apo.table.melt <- melt(apo.table)

# ggplot(subset(apo.table.melt, value>0), aes(x=Var2, y=Var1))+
#   geom_point(aes(size=value), shape=21, colour="black", fill="cornflowerblue")+
#   labs(x="population", y="MLG", size="# of individuals", title="Apo MLGs by Population")+
#   scale_size_area(max_size=15)+
#   theme_bw()


## @knitr DAPC.ms

my.pch <- c(17, 17, 17, 17, 17, 17, 21, 17, 17, 17, 17, 17, 17, 21, 21, 21, 21, 21, 17, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 17)
my.pch <-replace(my.pch,my.pch==21, 19)
my.pch.sub <- my.pch[-c(30,29,28,26)]
# for ms
setPop(AllPops.gc) <- ~ms
AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("S", "A"))
hookeri.dapc.ms <- dapc(AllPops.gc, grp=AllPops.gc$grp, n.pca=20, n.da=100)
scatter(hookeri.dapc.ms, grp = AllPops.gc$pop, cex = 2, legend = TRUE, clabel = T, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, col=c("red", "blue"))

## @knitr DAPC.popms

# all pops, but color by ms
setPop(AllPops.gc) <- ~pop
hookeri.dapc.msp <- dapc(AllPops.gc, grp=AllPops.gc$grp, n.pca=20, n.da=100)
scatter(hookeri.dapc.msp, grp = AllPops.gc$strata$ms, cex = 2, legend = TRUE, clabel = T,cstar=0, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch=c(17,19), col=c("red", "blue"))

## @knitr DAPC.popall

# for pops (all) 
setPop(AllPops.gc) <- ~pop
AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))
hookeri.dapc <- dapc(AllPops.gc, grp=AllPops.gc$grp, n.pca=20, n.da=100)
scatter(hookeri.dapc, grp = AllPops.gc$pop, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch=my.pch)

## @knitr DAPC.popsub

# sub a few pops
noYK.gc <- popsub(AllPops.gc, blacklist=c("C59-S", "SM-A", "C23-A", "S03-A"))
hookeri.dapc2 <- dapc(noYK.gc, grp=noYK.gc$pop, n.pca=20, n.da=100)
scatter(hookeri.dapc2, grp = noYK.gc$pop, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch=my.pch.sub)

## @knitr DAPC.member

summary(hookeri.dapc)

assignplot(hookeri.dapc)

compoplot(hookeri.dapc,posi="bottomright", lab="", ncol=1, xlab="individuals")

## @knitr compoplot

# structure-style plot
dapc.results <- as.data.frame(hookeri.dapc$posterior)
dapc.results$pop <- pop(AllPops.gc)
dapc.results$indNames <- rownames(dapc.results)

library(reshape2)
dapc.results <- melt(dapc.results)

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

# Plot posterior assignments from DAPC (how is this different from Structure?)
p4 <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p4 <- p4 + geom_bar(stat='identity') 
p4 <- p4 + scale_fill_manual(values = col_vector) 
p4 <- p4 + facet_grid(~Original_Pop, scales = "free")
p4 <- p4 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p4


## @knitr kmeans

# maxK <- 30
# myMat <- matrix(nrow=50, ncol=maxK)
# colnames(myMat) <- 1:ncol(myMat)
# for(i in 1:nrow(myMat)){
#   set.seed(9)
#   grp <- find.clusters(AllPops.gc, n.pca = 200, choose.n.clust = FALSE,  max.n.clust = maxK)
#   myMat[i,] <- grp$Kstat
# }
# 
# my_df2 <- melt(myMat)
# colnames(my_df2)[1:3] <- c("Group", "K", "BIC")
# my_df2$K <- as.factor(my_df2$K)
# #head(my_df2)
# 
# p1 <- ggplot(my_df2, aes(x = K, y = BIC))
# p1 <- p1 + geom_boxplot()
# p1 <- p1 + theme_bw()
# p1 <- p1 + xlab("Number of groups (K)")
# #p1 # 10 clusters kind of looks right
# 
# # plot 2 #
# my_k0 <- 10
# 
# grp0_l <- vector(mode = "list", length = length(my_k0))
# dapc0_l <- vector(mode = "list", length = length(my_k0))
# 
# for(i in 1:length(dapc0_l)){
#   set.seed(9)
#   grp0_l[[i]] <- find.clusters(AllPops.gc, n.pca = 200, n.clust = my_k0[i])
#   dapc0_l[[i]] <- dapc(AllPops.gc, pop = grp0_l[[i]]$grp, n.pca = 20, n.da = my_k0[i])
#   #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
# }
# 
# my_df0 <- as.data.frame(dapc0_l[[ length(dapc0_l) ]]$ind.coord)
# my_df0$Group <- dapc0_l[[ length(dapc0_l) ]]$grp
# #head(my_df0)
# 
# p2 <- ggplot(my_df0, aes(x = LD1, y = LD2, color = Group, fill = Group))
# p2 <- p2 + geom_point(size = 4, shape = 21)
# p2 <- p2 + theme_bw()
# p2 <- p2 + scale_color_manual(values=c(col_vector))
# p2 <- p2 + scale_fill_manual(values=c(paste(col_vector, "66", sep = "")))
# #p2
# 
# # plot 3
# my_k <- 10:13
# 
# grp_l <- vector(mode = "list", length = length(my_k))
# dapc_l <- vector(mode = "list", length = length(my_k))
# 
# for(i in 1:length(dapc_l)){
#   set.seed(9)
#   grp_l[[i]] <- find.clusters(AllPops.gc, n.pca = 200, n.clust = my_k[i])
#   dapc_l[[i]] <- dapc(AllPops.gc, pop = grp_l[[i]]$grp, n.pca = 20, n.da = my_k[i])
#   #  dapc_l[[i]] <- dapc(gl_rubi, pop = grp_l[[i]]$grp, n.pca = 3, n.da = 2)
# }
# 
# my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
# my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
# #head(my_df)
# 
# tmp <- as.data.frame(dapc_l[[1]]$posterior)
# tmp$K <- my_k[1]
# tmp$Sample <- rownames(tmp)
# tmp <- melt(tmp, id = c("Sample", "K"))
# names(tmp)[3:4] <- c("Group", "Posterior")
# tmp$pop <- mystrata$pop[match(tmp$Sample, mystrata$id)]
# my_df <- tmp
# 
# for(i in 2:length(dapc_l)){
#   tmp <- as.data.frame(dapc_l[[i]]$posterior)
#   tmp$K <- my_k[i]
#   tmp$Sample <- rownames(tmp)
#   tmp <- melt(tmp, id = c("Sample", "K"))
#   names(tmp)[3:4] <- c("Group", "Posterior")
#   tmp$pop <- mystrata$pop[match(tmp$Sample, mystrata$id)]
#   my_df <- rbind(my_df, tmp)
# }
# 
# grp.labs <- paste("K =", my_k)
# names(grp.labs) <- my_k
# 
# my_df$pop <- factor(my_df$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))
# 
# p3 <- ggplot(my_df, aes(x = Sample, y = Posterior, fill = Group))
# p3 <- p3 + geom_bar(stat = "identity")
# p3 <- p3 + facet_grid(K ~ pop, scales = "free_x", space = "free", 
#                       labeller = labeller(K = grp.labs))
# p3 <- p3 + theme_bw()
# p3 <- p3 + ylab("Posterior membership probability")
# p3 <- p3 + theme(legend.position='none')
# #p3 <- p3 + scale_color_brewer(palette="Paired")
# p3 <- p3 + scale_fill_manual(values=c(col_vector))
# p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
# #p3
# 
# ggarrange(ggarrange(p1,
#                     p2,
#                     ncol = 2, labels = c("A", "B")),
#           p3,
#           nrow = 2,
#           labels = c("", "C"),
#           heights = c(1, 2)
# )



## @knitr MSN
load("AllPops.gc.RData")
setPop(AllPops.gc, ~ms/pop)
msn <- poppr.msn(AllPops.gc, ddist2, showplot = FALSE)

# inds="none" to remove names
my.pch <- c(17, 17, 17, 17, 17, 17, 21, 17, 17, 17, 17, 17, 17, 21, 21, 21, 21, 21, 17, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 17)
my.pch <-replace(my.pch,my.pch==21, 19)
my.cols.ms <- replace(my.pch,my.pch==19, "blue")
my.cols.ms <- replace(my.cols.ms,my.cols.ms==17, "red")
replace(my.pch,my.pch==21, 19)

# inds="none" to remove names
plot_poppr_msn(AllPops.gc, msn, inds="none", palette=my.cols.ms, pop.leg = FALSE, size.leg = FALSE)





## @knitr AMOVA
hookeri.amova <- poppr.amova(AllPops.gc, ~ms/pop, within=FALSE, cutoff = 0.1)
hookeri.amova

## @knitr AMOVA.cc
hookeri.amova.cc <- poppr.amova(AllPops.gc, ~ms/pop, within=FALSE, cutoff = 0.1, clonecorrect = TRUE)
hookeri.amova.cc



## @knitr Poppr
# hookeri.poppr.ms <- poppr(AllPops.gc, sample=999, clonecorrect = TRUE, strata=~ms/pop/id)
# hookeri.poppr.ms

load("hookeri.poppr.pop.RData")
hookeri.poppr.pop




