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

# ploidy order of samples in vcf
vcfploidy <- c(2,2,2, 2,2,2, 2,2,2, 2,2,2, 2,2,2, 3,3,3,3,3, 3,3,3,3,3, 3,3,3,3,3, 2,2,2, 3,3,3,3,3, 3,3,3,3,3, 3,3,3,3,3, 3,3,3,3,3, 2,2,2, 3,3,3,3,3, 2,2,2, 2,2,2, 2,2,2, 2,2,2, 2,2,2, 3,3,3,3,3, 3,3,3,3,3, 3,3,3,3,3, 3,3,3,3, 3,3,3,3,2, 3,2,2,2,2, 3,3,3,3, 4,4,4,4,4)



# vcf to vcfR
MULTI.vcf <- read.vcfR("MULTI.final.filtered.snps.vcf")
#MULTI.vcf.dips <- read.vcfR("MULTI.filtered.dips.vcf")
#MULTI.vcf.trips <- read.vcfR("MULTI.filtered.trips.vcf")
#MULTI.vcf.tets <- read.vcfR("MULTI.filtered.tets.vcf")
OneMLG.vcf <- read.vcfR("OneMLG.vcf")

MULTI.gi <- vcfR2genind(MULTI.vcf, sep = "/", return.alleles = TRUE, ploidy=vcfploidy)
sampleorder <- match(indNames(MULTI.gi), mystrata$id)
strata(MULTI.gi) <- mystrata[sampleorder,]
setPop(MULTI.gi) <- ~pop
MULTI.gc <- as.genclone(MULTI.gi)

# set population levels (order from south to north)
save(MULTI.gc, file="MULTI.gc.RData")

apo.list <- c("C23-A", "C27-A", "C43-A", "C85-A", "C86-A", "C87-A", "C88-A", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A", "L45-A", "L62-A", "S03-A", "SM-A")
sex.list <- c("B42-S", "B46-S", "B49-S", "B53-S", "B60-S", "C59-S", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L45-S", "L62-S")

MULTI.dips.gc <- popsub(MULTI.gc, sublist=sex.list, drop=FALSE)
MULTI.apos.gc <- popsub(MULTI.gc, sublist=apo.list, drop=FALSE)

# # vcfR to genind to genclone. Keep return.alleles=TRUE change from 0/1 to A/T 
# dips.gi.ALL <- vcfR2genind(vcf.dips.ALL, sep = "/", ploidy=2, return.alleles = TRUE)
# #dips.gi <- vcfR2genind(vcf.dips, sep = "/", ploidy=2)
# dips.gc.ALL <- as.genclone(dips.gi.ALL)
# sampleorder <- match(indNames(dips.gc.ALL), mystrata$id)
# strata(dips.gc.ALL) <- mystrata[sampleorder,]
# setPop(dips.gc.ALL) <- ~pop
# save(dips.gc.ALL, file="dips.gc.ALL.RData")
# 
# trips.gi.ALL <- vcfR2genind(vcf.trips.ALL, sep = "/", ploidy=3, return.alleles = TRUE)
# #trips.gi <- vcfR2genind(vcf.trips, sep = "/", ploidy=3)
# trips.gc.ALL <- as.genclone(trips.gi.ALL)
# sampleorder <- match(indNames(trips.gc.ALL), mystrata$id)
# strata(trips.gc.ALL) <- mystrata[sampleorder,]
# setPop(trips.gc.ALL) <- ~pop
# save(trips.gc.ALL, file="trips.gc.ALL.RData")
# 
# tets.gi.ALL <- vcfR2genind(vcf.tets.ALL, sep = "/", ploidy=4, return.alleles = TRUE)
# #tets.gi <- vcfR2genind(vcf.tets, sep = "/", ploidy=4)
# tets.gc.ALL <- as.genclone(tets.gi.ALL)
# sampleorder <- match(indNames(tets.gc.ALL), mystrata$id)
# strata(tets.gc.ALL) <- mystrata[sampleorder,]
# setPop(tets.gc.ALL) <- ~pop
# save(tets.gc.ALL, file="tets.gc.ALL.RData")
# 
# # combine genclones
# dipsNtripsNtets.gc.ALL <- repool(dips.gc.ALL,trips.gc.ALL,tets.gc.ALL)
# #dipsNtripsNtets.gc$pop <- factor(dipsNtripsNtets.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))
# AllPops.gc.ALL <-as.genclone(dipsNtripsNtets.gc.ALL)

# set population levels (order from south to north)
MULTI.gc$pop <- factor(AllPops.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))
save(MULTI.gc, file="MULTI.gc.RData")

# gives info on ploidy per locus (confirm that poppr is reading alleles properly)
ploidy.infotable <- info_table(AllPops.gc, type="ploidy")
ploidy.infotable[60:90,30:60]

# looks at multi-locus lineages
mll(MULTI.gc)

apo.list <- c("C23-A", "C27-A", "C43-A", "C85-A", "C86-A", "C87-A", "C88-A", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A", "L45-A", "L62-A", "S03-A", "SM-A")
sex.list <- c("B42-S", "B46-S", "B49-S", "B53-S", "B60-S", "C59-S", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L45-S", "L62-S")

MULTI.dips.gc <- popsub(MULTI.gc, sublist=sex.list, drop=FALSE)
MULTI.apos.gc <- popsub(MULTI.gc, sublist=apo.list, drop=FALSE)

#OnePerMLGapos.gc <- apos.gc[apo.keep.list]

#OneMLG.gi <- repool(dips.gc, OnePerMLGapos.gc)
#OneMLG.gc <- as.genclone(OneMLG.gi)
#OneMLG.inds <- indNames(OneMLG.gc)

MULTI.OneMLG.gc <- MULTI.gc[OneMLG.inds]

save(MULTI.OneMLG.gc, file="MULTI.OneMLG.gc.RData")

load("OneMLG.gc.ALL.RData")
####################
# make more colors #
####################

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

cols <- brewer.pal(n = nPop(AllPops.gc.ALL), name = "Paired")

###################
# data conversion #
###################

# use vcfkeepsamples from vcflib to extract the inds from OneMLG.inds

# convert for polyrelatedness

my_gt <- extract.gt(OneMLG.vcf, convertNA = FALSE)
my_gt <- t(my_gt)

my_gt[] <- gsub("[/|]", "", my_gt) # remove separators
my_gt[] <- gsub("5", "6", my_gt) # 5 to 6
my_gt[] <- gsub("4", "5", my_gt) # 4 to 5
my_gt[] <- gsub("3", "4", my_gt) # 3 to 4
my_gt[] <- gsub("2", "3", my_gt) # 2 to 3
my_gt[] <- gsub("1", "2", my_gt) # change 1 to 2
my_gt[] <- gsub("0", "1", my_gt) # change 0 to 1
my_gt[] <- gsub("\\.", "0", my_gt) # change . to 0
my_gt[1:30, 1:30]

my_gt <-cbind(ind = rownames(my_gt), my_gt)
#rownames(my_gt) <- NULL
write.table(my_gt, file="in.txt", sep="\t", row.names=TRUE, quote=FALSE)

relat_mat <- read.csv("relatedness_mat.csv", header=FALSE)
relat_mat <- as.matrix(relat_mat)
ncol(relat_mat)
nrow(relat_mat)

inds <- rownames(my_gt)

colnames(relat_mat) <- inds
rownames(relat_mat) <- inds

relat_mat <- as.matrix(relat_mat)
