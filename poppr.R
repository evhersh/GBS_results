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

###############
# Data import #
###############

# strata (sample, pop, ms)
mystrata <- read.csv("~/Google Drive/GitHub/Hookeri-GBS/popmap_all.csv")

# vcf to vcfR
vcf <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/final.filtered.snps.vcf")
vcf.dips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.dips.vcf")
vcf.trips <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.trips.vcf")
vcf.tets <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.tets.vcf")

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

mll(AllPops.gc)
# save(AllPops.gc, file="AllPops.gc.RData")
# 
# x.mat <- as.matrix(AllPops.gc)
# x.mat[x.mat == 0] <- "1/1"
# x.mat[x.mat == 1] <- "1/2"
# x.mat[x.mat == 2] <- "2/2"
# 
# genind2genalex(AllPops.gc, "hookeri_genalex.txt", sep="\t", sequence=TRUE, overwrite = TRUE)

hook.fstat <- genind2hierfstat(AllPops.gc)
hookeri.AF <-genind2df(AllPops.gc)
hookeri.AF <-cbind(ind = rownames(hookeri.AF), hookeri.AF)
hookeri.AF[4]

rownames(hookeri.AF) <- NULL
write.table(hookeri.AF, file="Hookeri_AF.txt", sep="\t", row.names=TRUE, quote=FALSE)

# make more colors

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

cols <- brewer.pal(n = nPop(AllPops.gc), name = "Paired")

###################
# data conversion #
###################

genind2structure <- function(obj, file="", pops=FALSE){
  if(!"genind" %in% class(obj)){
    warning("Function was designed for genind objects.")
  }
  
  # get the max ploidy of the dataset
  pl <- max(obj@ploidy)
  # get the number of individuals
  S <- adegenet::nInd(obj)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind=rep(indNames(obj), each=pl))
  # column of pop ids to write
  if(pops){
    popnums <- 1:adegenet::nPop(obj)
    names(popnums) <- as.character(unique(adegenet::pop(obj)))
    popcol <- rep(popnums[as.character(adegenet::pop(obj))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- adegenet::locNames(obj) 
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=adegenet::nLoc(obj),
                           dimnames=list(NULL,loci)))
  
  # begin going through loci
  for(L in loci){
    thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), 
                              dimnames(obj@tab)[[2]]), 
                        drop = FALSE] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  
  # export table
  write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
}

# genind2structure(AllPops.gc, file="AllPops.structure", pops=TRUE)

############
# Clone ID #
############

# calculate raw euclidian distance
dist <- dist(AllPops.gc)

# assign MLG's using raw euclidian distance from dist() [above]
fstats <- filter_stats(AllPops.gc, distance=dist, plot=TRUE)

# looks like this gives the same clone mlg assignments as my IBS stuff
mlg.filter(AllPops.gc, distance=dist) <- 100


mlg.table(AllPops.gc)

#mll(AllPops.gc)

# indNames(AllPops.gc)

#########
# Poppr #
#########
set.seed(420)
load("hookeri.poppr.RData")
#hookeri.poppr <- poppr(AllPops.gc, sample=999, clonecorrect = TRUE, strata=~ms/pop/id)
#hookeri.poppr.pop <- poppr(AllPops.gc, sample=999, clonecorrect = TRUE, strata=~pop/id)


###################
# private alleles #
###################
my.ploidy <- c(17, 17, 17, 17, 17, 17, 21, 17, 17, 17, 17, 17, 17, 21, 21, 21, 21, 21, 17, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 17)
my.ploidy <-replace(my.ploidy,my.ploidy==21, 3)
my.ploidy <-replace(my.ploidy,my.ploidy==17, 2)
my.ploidy[29] <- 4
my.ploidy
## tabcount is a matrix pop x alleles, counting alleles per pop
tabcount <- apply(tab(AllPops.gc), 2, tapply, pop(AllPops.gc), sum, na.rm=FALSE)
AllPops.gp <- new("genpop", tabcount, type="codom", ploidy=my.ploidy)
summary(AllPops.gp)
popNames(All)

pal <- private_alleles(AllPops.gp, level="population", report="data.frame")



ggplot(pal) + geom_tile(aes(x = population, y = allele, fill = count))
ggplot(pal) + geom_boxplot(aes(x=population, y= count))

#########
# AMOVA #
#########
hookeri.amova <- poppr.amova(AllPops.gc, ~ms/pop, within=FALSE, cutoff = 0.1)
hookeri.amova.cc <- poppr.amova(AllPops.gc, ~ms/pop, within=FALSE, cutoff = 0.1, clonecorrect = TRUE)
hookeri.amova.cc.pop <- poppr.amova(AllPops.gc, ~pop, within=FALSE, cutoff = 0.1, clonecorrect = TRUE)

# pop_combinations <- combn(popNames(AllPops.gc), 2)
# amova_list.cc <- apply(pop_combinations, MARGIN = 2, function(i) poppr.amova(AllPops.gc[pop = i], ~pop, within=FALSE, cutoff=0.1, clonecorrect = TRUE))

#mlg.table(AllPops.gc)

#mll.custom(AllPops.gc) <- c(1,2,3, 4,5,6, 7,8,9, 10,11,12, 13,14,15, 16,17,18, 19,20,21, 22,23,24, 25,26,27, 28,29,30, 31,32,33, 34,35,36, 37,38,39,40, #B42 - L62 
                     #C23           #C27            #C43              #C85           #C86         #C87             #C88           #L06          #L16            #L17             #L39           #L41        #L45 change?    #L62    #S03           #SM
#                41,41,41,41,41, 42,42,42,42,42, 43,43,43,43,43, 42,42,42,42,42, 42,42,42,42,42, 42,42,42,42,42, 42,42,42,42,42, 44,44,44,44,44, 45,45,44,45,44, 44,44,44,44,44, 46,46,46,46,46, 42,42,42,42, 47,42,42,42,42, 47, 48,48,48,48, 49,49,49,49,49)

sampleorder <- match(indNames(AllPops.gc), mystrata$id)


#######
# Fst #
#######
my.ploidy <- c(17, 17, 17, 17, 17, 17, 21, 17, 17, 17, 17, 17, 17, 21, 21, 21, 21, 21, 17, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 17)
my.ploidy <-replace(my.ploidy,my.ploidy==21, 3)
my.ploidy <-replace(my.ploidy,my.ploidy==17, 2)
my.ploidy[29] <- 4

matFST <- pairwise.fst(AllPops.gc, res.type="matrix")

Nei.mat <- dist.genpop(AllPops.gp, method=1)
Nei.tree <- nj(Nei.mat)

plot(Nei.tree, type="unr", tip.col=my.cols.ms, font=2)
annot <- round(Nei.tree$edge.length,2)
edgelabels(annot[annot>0], which(annot>0), frame="n")
add.scale.bar()
##################
##### DAPC #######
##################



  ###########
  # K-means #
  ###########

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
  set.seed(9)
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
my_k0 <- 10

grp0_l <- vector(mode = "list", length = length(my_k0))
dapc0_l <- vector(mode = "list", length = length(my_k0))

for(i in 1:length(dapc0_l)){
  set.seed(9)
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
my_k <- 10:13

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
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

my_df[ my_df$K == 2 & my_df$Group == 2, "Group"]

p3 <- ggplot(my_df, aes(x = Sample, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ pop, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
#p3 <- p3 + scale_color_brewer(palette="Paired")
p3 <- p3 + scale_fill_manual(values=c(col_vector))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3

ggarrange(ggarrange(p1,
                    p2,
                    ncol = 2, labels = c("A", "B")),
          p3,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)


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
scatter(hookeri.dapc.msp, grp = AllPops.gc$strata$ms, cex = 2, legend = TRUE, clabel = T, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch=c(17,19))

# for pops (all)
setPop(AllPops.gc) <- ~pop
AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))
hookeri.dapc <- dapc(AllPops.gc, grp=AllPops.gc$grp, n.pca=20, n.da=100)
scatter(hookeri.dapc, grp = AllPops.gc$pop, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch=my.pch)

# sub a few pops
noYK.gc <- popsub(AllPops.gc, blacklist=c("C59-S", "SM-A", "C23-A", "S03-A"))
hookeri.dapc2 <- dapc(noYK.gc, grp=noYK.gc$pop, n.pca=20, n.da=100)
scatter(hookeri.dapc2, grp = noYK.gc$pop, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75, pch=my.pch.sub)

# loadings
contib <- loadingplot(hookeri.dapc$var.contr, axis=1, thres=0.0003, lab.jitter = 1)

# membership probabilities
summary(hookeri.dapc)

assignplot(hookeri.dapc)

compoplot(hookeri.dapc,posi="bottomright", lab="", ncol=1, xlab="individuals")

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


# create figure
# tiff("DAPC.tiff", res=300, units="in", width=8.6, height=5.8)
# scatter(hookeri.dapc2, grp = noYK.gc$pop, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)
# dev.off()

#seasun(length(levels(mdapc$grp)))[mdapc$grp]
#col = col_vector
#grp = dipsNtripsNtets.gc@strata$region




#########
# Trees #
#########

# aboot(AllPops.gc, dist = dist(), sample = 200, tree = "nj", cutoff = 50, quiet = FALSE)
# 
# theTree <- dist %>%
#   nj() %>%    # calculate neighbor-joining tree
#   ladderize() # organize branches by clade
# plot(theTree)
# add.scale.bar(length = 0.05)

# for inds
hookeri.nj <- aboot(AllPops.gc, dist = provesti.dist, sample = 200, tree = "nj", cutoff = 50, quiet = TRUE)

# for pops
hookpop.nj <- aboot(AllPops.gp, dist = provesti.dist, sample = 200, tree = "nj", cutoff = 50, quiet = TRUE)

plot.phylo(hookpop.nj, cex=0.8, tip.color = my.cols.ms)
nodelabels(hookpop.nj$node.label, adj = c(1.5, -0.7), frame = "n", cex = 0.8,
           font = 3, xpd = TRUE)
axisPhylo(3)

# AllPops.gc %>%
#   genind2genpop(pop = ~pop) %>%
#   aboot(cutoff = 50, quiet = FALSE, sample = 1000, distance = nei.dist)




#######
# MSN #
#######
setPop(AllPops.gc, ~ms/pop)
msn <- poppr.msn(AllPops.gc, dist, showplot = FALSE)

# inds="none" to remove names
my.cols.ms <- replace(my.pch,my.pch==19, "red")
my.cols.ms <- replace(my.cols.ms,my.cols.ms==17, "blue")
replace(my.pch,my.pch==21, 19)

# inds="none" to remove names
plot_poppr_msn(AllPops.gc, msn, inds="none", palette=my.cols.ms)
