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

# # read data into vcfR
# vcf.unsplit <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/unsplit.apos.vcf")
# 
# # convert to genind/genclone
# unsplit.gi <- vcfR2genind(vcf.unsplit, sep = "/", ploidy=3, return.alleles = TRUE)
# unsplit.gc <- as.genclone(unsplit.gi)
# save(unsplit.gc, file="unsplit.gc.RData")
load("unsplit.gc.RData")

# calculate distance
distmat.unsplit <- as.matrix(dist(unsplit.gc))
labels(distmat.unsplit)

# C23
rep.dists <- distmat.unsplit["C23-A_1A", "C23-A_1B"]
rep.dists <- c(rep.dists, distmat.unsplit["C23-A_2A", "C23-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["C23-A_3A", "C23-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["C23-A_4A", "C23-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["C23-A_5A", "C23-A_5B"])

# C27
rep.dists <- c(rep.dists, distmat.unsplit["C27-A_1A", "C27-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["C27-A_2A", "C27-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["C27-A_3A", "C27-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["C27-A_4A", "C27-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["C27-A_5A", "C27-A_5B"])

# C43
rep.dists <- c(rep.dists, distmat.unsplit["C43-A_1A", "C43-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["C43-A_2A", "C43-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["C43-A_3A", "C43-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["C43-A_4A", "C43-A_4B"]) # 
rep.dists <- c(rep.dists, distmat.unsplit["C43-A_5A", "C43-A_5B"])

# C85
rep.dists <- c(rep.dists, distmat.unsplit["C85-A_1A", "C85-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["C85-A_2A", "C85-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["C85-A_3A", "C85-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["C85-A_4A", "C85-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["C85-A_5A", "C85-A_5B"])

# C86
rep.dists <- c(rep.dists, distmat.unsplit["C86-A_1A", "C86-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["C86-A_2A", "C86-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["C86-A_3A", "C86-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["C86-A_4A", "C86-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["C86-A_5A", "C86-A_5B"])

# C87
rep.dists <- c(rep.dists, distmat.unsplit["C87-A_1A", "C87-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["C87-A_2A", "C87-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["C87-A_3A", "C87-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["C87-A_4A", "C87-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["C87-A_5A", "C87-A_5B"])

# C88
rep.dists <- c(rep.dists, distmat.unsplit["C88-A_1A", "C88-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["C88-A_2A", "C88-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["C88-A_3A", "C88-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["C88-A_4A", "C88-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["C88-A_5A", "C88-A_5B"])

# L06
rep.dists <- c(rep.dists, distmat.unsplit["L06-A_1A", "L06-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["L06-A_2A", "L06-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["L06-A_3A", "L06-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["L06-A_4A", "L06-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["L06-A_5A", "L06-A_5B"])

# L16
rep.dists <- c(rep.dists, distmat.unsplit["L16-A_1A", "L16-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["L16-A_2A", "L16-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["L16-A_3A", "L16-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["L16-A_4A", "L16-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["L16-A_5A", "L16-A_5B"])

# L17
rep.dists <- c(rep.dists, distmat.unsplit["L17-A_1A", "L17-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["L17-A_2A", "L17-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["L17-A_3A", "L17-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["L17-A_4A", "L17-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["L17-A_5A", "L17-A_5B"])

# L39
rep.dists <- c(rep.dists, distmat.unsplit["L39-A_1A", "L39-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["L39-A_2A", "L39-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["L39-A_3A", "L39-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["L39-A_4A", "L39-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["L39-A_5A", "L39-A_5B"])

# L41
rep.dists <- c(rep.dists, distmat.unsplit["L41-A_1A", "L41-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["L41-A_2A", "L41-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["L41-A_3A", "L41-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["L41-A_4A", "L41-A_4B"])

# L45
rep.dists <- c(rep.dists, distmat.unsplit["L45-A_1A", "L45-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["L45-A_2A", "L45-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["L45-A_3A", "L45-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["L45-A_4A", "L45-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["L45-A_5A", "L45-A_5B"])

# L62
rep.dists <- c(rep.dists, distmat.unsplit["L62-A_5A", "L62-A_5B"])

# S03
rep.dists <- c(rep.dists, distmat.unsplit["S03-A_1A", "S03-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["S03-A_2A", "S03-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["S03-A_3A", "S03-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["S03-A_4A", "S03-A_4B"])

mean(rep.dists[-c(14,17,26)]) # 54.67
hist(rep.dists)
