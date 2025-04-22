library(gbs2ploidy)
library(vcfR)

vcf <- read.vcfR("./data/final.filtered.snps.vcf", verbose = TRUE)

# Suppose you have AD as a character matrix like "12,8"
ad <- extract.gt(vcf, element = "AD", as.numeric = FALSE)

# Split AD into ref and alt matrices
ad_split <- strsplit(ad, ",")
refmat <- matrix(NA, nrow = nrow(ad), ncol = ncol(ad))
altmat <- matrix(NA, nrow = nrow(ad), ncol = ncol(ad))
rownames(refmat) <- rownames(altmat) <- rownames(ad)
colnames(refmat) <- colnames(altmat) <- colnames(ad)

for (i in 1:length(ad_split)) {
  val <- ad_split[[i]]
  if (length(val) >= 2) {
    refmat[i] <- as.numeric(val[1])
    altmat[i] <- as.numeric(val[2])
  }
}

# Transpose matrices to have samples as rows and loci as columns
cov1 <- t(refmat)
cov2 <- t(altmat)

# Estimate allelic proportions
props <- estprops(cov1 = cov1, cov2 = cov2, mcmc.steps = 1000, mcmc.burnin = 500, mcmc.thin = 10)

# Calculate heterozygosity
het <- apply(!is.na(cov1 + cov2), 1, mean)

# Calculate depth
depth <- apply(cov1 + cov2, 1, mean, na.rm = TRUE)

# Infer ploidy without a training set
ploidy_results <- estploidy(alphas = props, het = het, depth = depth, train = FALSE, nclasses = 2, ids = rownames(cov1), pcs = 1:2)

# View posterior probabilities
head(ploidy_results$pp)