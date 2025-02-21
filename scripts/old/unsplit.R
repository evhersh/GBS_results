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
library(gmodels)
library(tidyverse)

# # read data into vcfR


vcf.unsplit1 <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.unsplit.trips.vcf")
vcf.unsplit2 <- read.vcfR("~/Google Drive/GitHub/Hookeri-GBS/Data/filtered.unsplit.tets.vcf")

# 
# # convert to genind/genclone
unsplit1.gi <- vcfR2genind(vcf.unsplit1, sep = "/", ploidy=3, return.alleles = TRUE)
unsplit2.gi <- vcfR2genind(vcf.unsplit2, sep = "/", ploidy=4, return.alleles = TRUE)

unsplit3.gi <- repool(unsplit1.gi, unsplit2.gi)
unsplit.gc <- as.genclone(unsplit3.gi)
save(unsplit.gc, file="unsplit.gc.RData")
load("unsplit.gc.RData")

# calculate distance
#distmat.unsplit <- as.matrix(dist(unsplit.gc))
distmat.unsplit <- as.matrix(prevosti.dist(unsplit.gc))
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
rep.dists <- c(rep.dists, distmat.unsplit["C43-A_4A", "C43-A_4B"]) # dist of 0.088
rep.dists <- c(rep.dists, distmat.unsplit["C43-A_5A", "C43-A_5B"])

# C85
rep.dists <- c(rep.dists, distmat.unsplit["C85-A_1A", "C85-A_1B"])
#rep.dists <- c(rep.dists, distmat.unsplit["C85-A_2A", "C85-A_2B"])
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
rep.dists <- c(rep.dists, distmat.unsplit["C87-A_1A", "C87-A_1B"]) # dist of 0.102
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
rep.dists <- c(rep.dists, distmat.unsplit["L16-A_1A", "L16-A_1B"]) # dist of 0.09
rep.dists <- c(rep.dists, distmat.unsplit["L16-A_2A", "L16-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["L16-A_3A", "L16-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["L16-A_4A", "L16-A_4B"]) # dist of 0.075
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

# SM
rep.dists <- c(rep.dists, distmat.unsplit["SM-A_1A", "SM-A_1B"])
rep.dists <- c(rep.dists, distmat.unsplit["SM-A_2A", "SM-A_2B"])
rep.dists <- c(rep.dists, distmat.unsplit["SM-A_3A", "SM-A_3B"])
rep.dists <- c(rep.dists, distmat.unsplit["SM-A_4A", "SM-A_4B"])
rep.dists <- c(rep.dists, distmat.unsplit["SM-A_5A", "SM-A_5B"])

mean(rep.dists[-c(14,17,26)]) # 54.67
hist(rep.dists)

error.df <- as.data.frame(rep.dists)
error.df$ind <- 1:length(rep.dists)

error.mean <- error.df %>%
  summarize(mean=mean(rep.dists), sd=sd(rep.dists), n=n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)


ggplot(error.df, aes(y=rep.dists, x=ind))+
  geom_jitter()


