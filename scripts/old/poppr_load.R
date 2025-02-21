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


head(is.biallelic(vcf.dips))


# vcfR to genind to genclone. Keep return.alleles=TRUE change from 0/1 to A/T 
dips.gi <- vcfR2genind(vcf.dips, sep = "/", ploidy=2, return.alleles = TRUE)
#dips.gi <- vcfR2genind(vcf.dips, sep = "/", ploidy=2)
dips.gc <- as.genclone(dips.gi)
sampleorder <- match(indNames(dips.gc), mystrata$id)
sampleorder <- match(indNames(dips.gi), mystrata$id)
strata(dips.gc) <- mystrata[sampleorder,]
setPop(dips.gc) <- ~pop
strata(dips.gi) <- mystrata[sampleorder,]
setPop(dips.gi) <- ~pop
save(dips.gc, file="dips.gc.RData")
save(dips.gi, file="dips.gi.RData")


trips.gi <- vcfR2genind(vcf.trips, sep = "/", ploidy=3, return.alleles = TRUE)
#trips.gi <- vcfR2genind(vcf.trips, sep = "/", ploidy=3)
trips.gc <- as.genclone(trips.gi)
sampleorder <- match(indNames(trips.gc), mystrata$id)
strata(trips.gc) <- mystrata[sampleorder,]
setPop(trips.gc) <- ~pop
save(trips.gc, file="trips.gc.RData")

tets.gi <- vcfR2genind(vcf.tets, sep = "/", ploidy=4, return.alleles = TRUE)
#tets.gi <- vcfR2genind(vcf.tets, sep = "/", ploidy=4)
tets.gc <- as.genclone(tets.gi)
sampleorder <- match(indNames(tets.gc), mystrata$id)
strata(tets.gc) <- mystrata[sampleorder,]
setPop(tets.gc) <- ~pop
save(tets.gc, file="tets.gc.RData")

# combine genclones
dipsNtripsNtets.gc <- repool(dips.gc,trips.gc,tets.gc)
#dipsNtripsNtets.gc$pop <- factor(dipsNtripsNtets.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))
AllPops.gc <-as.genclone(dipsNtripsNtets.gc)

# set population levels (order from south to north)
AllPops.gc$pop <- factor(AllPops.gc$pop, levels=c("B53-S", "B60-S", "B42-S", "B46-S", "B49-S", "L62-S", "L62-A", "L05-S", "L08-S", "L10-S", "L11-S", "L12-S", "L13-S", "L06-A", "L16-A", "L17-A", "L39-A", "L41-A","L45-S", "L45-A", "C87-A", "C86-A", "C88-A", "C85-A", "C27-A", "C23-A", "C43-A", "S03-A", "SM-A", "C59-S"))
save(AllPops.gc, file="AllPops.gc.RData")

# gives info on ploidy per locus (confirm that poppr is reading alleles properly)
ploidy.infotable <- info_table(AllPops.gc, type="ploidy")
ploidy.infotable[60:90,30:60]

# looks at multi-locus lineages
mll(AllPops.gc)

dips.gc <- popsub(AllPops.gc, sublist=sex.list, drop=FALSE)
apos.gc <- popsub(AllPops.gc, sublist=apo.list, drop=FALSE)

#OnePerMLGapos.gc <- apos.gc[apo.keep.list]

#OneMLG.gi <- repool(dips.gc, OnePerMLGapos.gc)
#OneMLG.gc <- as.genclone(OneMLG.gi)
#OneMLG.inds <- indNames(OneMLG.gc)

OneMLG.gc <- AllPops.gc[OneMLG.inds]

save(OneMLG.gc, file="OneMLG.gc.RData")

load("OneMLG.gc.RData")
####################
# make more colors #
####################

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

cols <- brewer.pal(n = nPop(AllPops.gc), name = "Paired")

###################
# data conversion #
###################

# use vcfkeepsamples from vcflib to extract the inds from OneMLG.inds

# convert for polyrelatedness

my_gt <- extract.gt(vcf, convertNA = FALSE)
my_gt <- t(my_gt)

my_gt[] <- gsub("[/|]", "", my_gt) # remove separators
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
