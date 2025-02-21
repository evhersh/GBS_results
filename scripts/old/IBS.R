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
require(reshape2)
require(Rmisc)
library(SeqArray)

# SNPRelate
vcf.unsplit <-"~/Google Drive/GitHub/Hookeri-GBS/Data/nomissing.unsplit.apos.vcf"
vcf.fn <-"~/Google Drive/GitHub/Hookeri-GBS/Data/final.filtered.snps.vcf"
snpgdsVCF2GDS(vcf.fn, "test1.gds")
snpgdsSummary("test1.gds")
genofile <- snpgdsOpen("test1.gds")
#pop_code <- scan("~/Google Drive/GitHub/Hookeri-GBS/pop_code_unsplit.txt", what=character())
pop_code <- scan("~/Google Drive/GitHub/Hookeri-GBS/pop_code.txt", what=character())
table(pop_code)
head(pop_code)
samp_id <- read.gdsn(index.gdsn(genofile, "sample.id"))

pca <- snpgdsPCA(genofile, snp.id=genofile$snp.id, num.thread=2, autosome.only = FALSE)

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
head(cbind(sample.id, pop_code))

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))

# SNPRelate IBS
ibs <- snpgdsIBS(genofile, num.thread=2, autosome.only = FALSE)
distance.ibs <- as.dist(ibs$ibs)
sample.id

B42S.ibs <- as.dist(ibs$ibs[1:3,1:3]) #  83-84%
B46S.ibs <- as.dist(ibs$ibs[4:6,4:6]) #  79-84%
B49S.ibs <- as.dist(ibs$ibs[7:9,7:9]) #  83-85%
B53S.ibs <- as.dist(ibs$ibs[10:12,10:12]) #  83%
B60S.ibs <- as.dist(ibs$ibs[13:15,13:15]) #  85-86%
L05S.ibs <- as.dist(ibs$ibs[54:56,54:56]) #  79-86%
L08S.ibs <- as.dist(ibs$ibs[62:64,62:64]) #  80%
L10S.ibs <- as.dist(ibs$ibs[65:67,65:67]) #  81%
L11S.ibs <- as.dist(ibs$ibs[68:70,68:70]) #  80%
L12S.ibs <- as.dist(ibs$ibs[71:73,71:73]) #  80%
L13S.ibs <- as.dist(ibs$ibs[74:76,74:76]) #  80%
L62AS.ibs <- as.dist(ibs$ibs[101:105,101:105]) # 80-81%
C59S.ibs <- as.dist(ibs$ibs[31:33,31:33]) # 92-94% ****

wInPopIBS.sex <- as.vector(c(B42S.ibs, B46S.ibs, B49S.ibs,B53S.ibs,B60S.ibs,L05S.ibs,L08S.ibs,L10S.ibs,L11S.ibs,L12S.ibs,L13S.ibs,C59S.ibs))
wInPopIBS.sex$ms <- rep("sex",nrow(wInPopIBS.sex))
colnames(wInPopIBS.sex) <- c("ibs","ms")

wInPopIBS.apo$ms <- rep("apo",nrow(wInPopIBS.apo))
colnames(wInPopIBS.apo) <- c("ibs","ms")


#aes(y=..density..),

# gg.wInPopIBS.sex <- ggplot(data=wInPopIBS.sex, aes(x=ibs))+
#   geom_histogram(bins=25, alpha=0.75, fill="lightblue", color="black")+
#   theme_classic()+
#   xlim(.75,1.00)
# 
# gg.wInPopIBS.apo <- ggplot(data=wInPopIBS.apo, aes(x=ibs))+
#   geom_histogram(bins=30, alpha=0.75, fill="red", color="black")+
#   theme_classic()+
#   xlim(.75,1.00)

L06A.ibs <- as.dist(ibs$ibs[57:61,57:61]) #  97% (all same)
L16A.ibs <- as.dist(ibs$ibs[77:81,77:81]) #  84-96% (different)****
L17A.ibs <- as.dist(ibs$ibs[82:86,82:86]) #  97% (all same)
L39A.ibs <- as.dist(ibs$ibs[87:91,87:91]) #  87-97% (different)****
L41A.ibs <- as.dist(ibs$ibs[92:95,92:95]) #  96-97% (all same)
L45A.ibs <- as.dist(ibs$ibs[96:100,96:100]) # 74--95% (different)****
C23A.ibs <- as.dist(ibs$ibs[16:20,16:20]) #  96-97% (all same)
C27A.ibs <- as.dist(ibs$ibs[21:25,21:25]) #  97% (all same)
C43A.ibs <- as.dist(ibs$ibs[26:30,26:30]) #  94-97% (more or less the same?)
C85A.ibs <- as.dist(ibs$ibs[34:38,34:38]) #  96-97% (all same)
C86A.ibs <- as.dist(ibs$ibs[39:43,39:43]) #  97% (all same)
C87A.ibs <- as.dist(ibs$ibs[44:48,44:48]) #  96-97% (all same)
C88A.ibs <- as.dist(ibs$ibs[49:53,49:53]) #  97% (all same)
S03A.ibs <- as.dist(ibs$ibs[106:109,106:109]) #  97% (all same)
SMA.ibs <- as.dist(ibs$ibs[110:114,110:114]) #   96-97% (all same)

wInPopIBS.apo <- as.vector(c(L06A.ibs,L16A.ibs,L17A.ibs,L39A.ibs,L41A.ibs,L45A.ibs,C23A.ibs,C27A.ibs,C43A.ibs,C85A.ibs,C86A.ibs,C87A.ibs,C88A.ibs,S03A.ibs,SMA.ibs))
hist(wInPopIBS.apo)

hist(distance.ibs)

wInPopIBS <- rbind(wInPopIBS.sex,wInPopIBS.apo)


gg.wInPopIBS <- ggplot(data=wInPopIBS, aes(x=ibs, fill=ms))+
  geom_histogram(aes(y=..density..), bins=50, alpha=0.6, color="black", position="identity", size=.75)+
  theme_classic()+
  scale_x_continuous(breaks=seq(0.70, 1, 0.025))+
  scale_fill_manual(values=c("red", "blue"), name="Mating System", labels=c("Apomictic", "Sexual"))+
  geom_density(alpha=.3)+
  geom_vline(aes(xintercept=0.9681932), color="black", linetype="dashed")+
  labs(x="Indentity by state (IBS)", y="Density")

tiff("wInPopIBS.tiff", res=300, units="in", width=8.6, height=5.8)
gg.wInPopIBS
dev.off()

#### pop comparisons ####

## Laramie group ##

# L06:L16
ibs$ibs[57:61,77:81]

# L16:L17
ibs$ibs[77:81,82:86]

# L06:L17
ibs$ibs[57:61,82:86]

# L45:C27
ibs$ibs[96:100,21:25]

### for unsplit apos ###
vcf.unsplit <-"~/Google Drive/GitHub/Hookeri-GBS/Data/nomissing.unsplit.apos.vcf"
snpgdsVCF2GDS(vcf.unsplit, "test2.gds")
snpgdsSummary("test2.gds")
genofile2 <- snpgdsOpen("test2.gds")
pop_code2 <- scan("~/Google Drive/GitHub/Hookeri-GBS/pop_code_unsplit.txt", what=character())
table(pop_code2)
head(pop_code2)
samp_id2 <- read.gdsn(index.gdsn(genofile2, "sample.id"))
head(cbind(samp_id2, pop_code2))

ibs2 <- snpgdsIBS(genofile2, num.thread=2, autosome.only = FALSE)

even_index <- seq(2, 138, 2)
odd_index <- seq(1,137, 2)

sample.pairs <- samp_id2[seq(1,138,2)]

datalist=list()
for(i in c(1:69)){
  dat <- ibs2$ibs[odd_index[i],even_index[i]]
  datalist[[i]] <- dat
}

replicate.ibs <- as.numeric(datalist)
replicate.ibs <- replicate.ibs[-17]
hist(replicate.ibs)

replicate.ibs <- data.frame(replicate.ibs)
mean(replicate.ibs$replicate.ibs)

replicate.ibs$rep <- seq.int(nrow(replicate.ibs))
replicate.ibs$grp <- rep("replicate.pairs",nrow(replicate.ibs))

repSE <- summarySE(replicate.ibs, measurevar="replicate.ibs", groupvars="grp")

# replicate ibs (error rates)
pd = position_dodge(0.1)
gg.replicateIBS <- ggplot(data=replicate.ibs, aes(x=grp, y=replicate.ibs))+
  geom_jitter(width=0.05, alpha=.5, shape=21, size=2, colour="black", fill="darkred")+
  geom_point(data=repSE, aes(y=replicate.ibs), size=3)+
  geom_errorbar(data=repSE, aes(ymin=replicate.ibs-ci, ymax=replicate.ibs+ci), width = 0.1, position=pd, size=.75)+
  theme_classic()+
  ylim(values=c(.9, 1))+
  labs(x="Replicate pairs (apomictic)", y="Identity by State (IBS)")+
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())

#geom_hline(aes(yintercept=0.9681932), linetype="dashed"

tiff("replicateIBS.tiff", res=300, units="in", width=5.8, height=8.6)
gg.replicateIBS
dev.off()

# big_data = do.call(rbind, datalist)
# big_data <- cbind(big_data,sample.pairs)
# colnames(big_data) <- c("ibs","sample.pairs")
# big_data <- as.data.frame(big_data)
# big_data$ibs <- as.numeric(big_data$ibs)
# small.ibs <- big_data %>%
#   filter(ibs < 0.95)
# mean(big_data)



# pop.idx <- order(pop_code)
# image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(16))
# 
# loc <- cmdscale(1 - ibs$ibs, k = 2)
# x <- loc[, 1]; y <- loc[, 2]
# race <- as.factor(pop_code)
# 
# plot(x, y, col=race, xlab = "", ylab = "",
#      main = "Multidimensional Scaling Analysis (IBS)")
# legend("topleft", legend=levels(race), text.col=1:nlevels(race))
# 
# set.seed(100)
# ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2, autosome.only = FALSE))
# rv <- snpgdsCutTree(ibs.hc)
# plot(rv$dendrogram, leaflab="none", main="HapMap Phase II")
# 
# table(rv$samp.group)
# 
# rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))
# 
# plot(rv2$dendrogram, leaflab="none", main="HapMap Phase II")
# legend("topright", legend=levels(race), col=1:nlevels(race), pch=19, ncol=4)