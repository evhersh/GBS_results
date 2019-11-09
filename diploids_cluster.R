library("poppr")
library("hierfstat")
library(pegas)
library(adegenet)

load("dips.gc.RData")

##### HW test and F-statistics

dips.gp <- genind2genpop(dips.gc)

dips.hwt <- hw.test(dips.gc)
dips.hwt

fstat.dips <- fstat(dips.gc)

Gtest.dips <- gstat.randtest(dips.gc,nsim=99)
Gtest.dips
plot(Gtest.dips)

matFst.dips <- pairwise.fst(dips.gc)

##### inbreeding
seppop.dips <- seppop(dips.gc)
seppop.dips

temp <- inbreeding(seppop.dips, N=100)

Fbar <- sapply(temp, mean)

save.image(file='diploid_stats.RData')

q()
