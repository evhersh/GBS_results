library("poppr")
library("hierfstat")
library(pegas)
library(adegenet)
install.packages("genetics")
library(genetics)

load("dips.gc.RData")

##### HW test and F-statistics

dips.gp <- genind2genpop(dips.gc)

dips.summary <- summary(dips.gc)

to_keep <- propTyped(dips.gc, by = "loc") > .99

x <- dips.gc[loc = to_keep]

dips.hwt <- hw.test(dips.gc)
#dips.hwt

fstat.dips <- fstat(dips.gc)
fstat.dips

Gtest.dips <- gstat.randtest(x,nsim=99)
Gtest.dips
plot(Gtest.dips)

matFst.dips <- pairwise.fst(dips.gc)
matFst.dips
##### inbreeding

temp <- inbreeding(dips.gc, N=100)

Fbar <- sapply(temp, mean)

save.image(file='diploid_stats.RData')

##### AMOVA

dips.amova <- poppr.amova(dips.gc, ~pop)
dips.amova

dips.amova.sig <- randtest(dips.amova, nrepet = 99)
dips.amova.sig


q()


