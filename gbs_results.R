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


# Data import 

# strata (sample, pop, ms)
mystrata <- read.csv("~/Google Drive/GitHub/Hookeri-GBS/popmap_all.csv")

# vcf to vcfR
vcf <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/final.filtered.snps.vcf")
vcf.dips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.dips.vcf")
vcf.trips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.trips.vcf")
vcf.tets <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.tets.vcf")

# vcfR to genind to genclone
dips.gi <- vcfR2genind(vcf.dips, sep = "/", ploidy=2)
dips.gc <- as.genclone(dips.gi)
sampleorder <- match(indNames(dips.gc), mystrata$id)
strata(dips.gc) <- mystrata[sampleorder,]
setPop(dips.gc) <- ~pop

trips.gi <- vcfR2genind(vcf.trips, sep = "/", ploidy=3)
trips.gc <- as.genclone(trips.gi)
sampleorder <- match(indNames(trips.gc), mystrata$id)
strata(trips.gc) <- mystrata[sampleorder,]
setPop(trips.gc) <- ~pop

tets.gi <- vcfR2genind(vcf.tets, sep = "/", ploidy=4)
tets.gc <- as.genclone(tets.gi)
sampleorder <- match(indNames(tets.gc), mystrata$id)
strata(tets.gc) <- mystrata[sampleorder,]
setPop(tets.gc) <- ~pop

# combine genclones
dipsNtripsNtets.gc <- repool(dips.gc,trips.gc,tets.gc)
#dipsNtripsNtets.gc$pop <- factor(dipsNtripsNtets.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))

AllPops.gc <-as.genclone(dipsNtripsNtets.gc)
AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))

mll(AllPops.gc)






## @knitr MLGs

# calculate raw euclidian distance
dist <- dist(AllPops.gc)

# assign MLG's using raw euclidian distance from dist() [above]
fstats <- filter_stats(AllPops.gc, distance=dist, plot=TRUE)

# looks like this gives the same clone mlg assignments as my IBS stuff
mlg.filter(AllPops.gc, distance=dist) <- 100

mlgtab <- mlg.table(AllPops.gc)




## @knitr DAPC.ms

my.pch <- c(17, 17, 17, 17, 17, 17, 21, 17, 17, 17, 17, 17, 17, 21, 21, 21, 21, 21, 17, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 17)
my.pch <-replace(my.pch,my.pch==21, 19)
my.pch.sub <- my.pch[-c(30,29,28,26)]
# for ms
setPop(AllPops.gc) <- ~ms
AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("S", "A"))
hookeri.dapc.ms <- dapc(AllPops.gc, grp=AllPops.gc$grp, n.pca=20, n.da=100)
scatter(hookeri.dapc.ms, grp = AllPops.gc$pop, cex = 2, legend = TRUE, clabel = T, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)

## @knitr DAPC.popms

# all pops, but color by ms
setPop(AllPops.gc) <- ~pop
hookeri.dapc.msp <- dapc(AllPops.gc, grp=AllPops.gc$grp, n.pca=20, n.da=100)
scatter(hookeri.dapc.msp, grp = AllPops.gc$strata$ms, cex = 2, legend = TRUE, clabel = T, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch=c(17,19))

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

# structure-style plot
dapc.results <- as.data.frame(hookeri.dapc$posterior)
dapc.results$pop <- pop(AllPops.gc)
dapc.results$indNames <- rownames(dapc.results)

library(reshape2)
dapc.results <- melt(dapc.results)

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

# Plot posterior assignments from DAPC (how is this different from Structure?)
p2 <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p2 <- p2 + geom_bar(stat='identity') 
p2 <- p2 + scale_fill_manual(values = col_vector) 
p2 <- p2 + facet_grid(~Original_Pop, scales = "free")
p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p2






## @knitr MSN

setPop(AllPops.gc, ~ms/pop)
msn <- poppr.msn(AllPops.gc, dist, showplot = FALSE)

# inds="none" to remove names
my.cols.ms <- replace(my.pch,my.pch==19, "red")
my.cols.ms <- replace(my.cols.ms,my.cols.ms==17, "blue")
replace(my.pch,my.pch==21, 19)

# inds="none" to remove names
plot_poppr_msn(AllPops.gc, msn, inds="none", palette=my.cols.ms)
