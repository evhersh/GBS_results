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

# Load data
vcf.MA <- read.vcfR("MULTI.final.filtered.snps.vcf")
is.biallelic(vcf.MA)

gt_MA <- extract.gt(vcf.MA, convertNA = TRUE)
gt_MA[880:900, 1:30]

gt_MA_transposed <- t(gt_MA)

write.table(gt_MA_transposed, file="snp_matrix_transposed.txt", sep="\t", row.names=TRUE, quote=FALSE)
