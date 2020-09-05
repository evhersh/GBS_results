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
library(dplyr)
library(plotrix)



###



#AllPops.inds <- as.list(rep(NA, 114))
ind.Hobs <- data.frame(matrix(NA, nrow=114, ncol=3))
colnames(ind.Hobs) <- c("ind", "ploidy", "Hobs")
ind.Hobs$ms <- "NA"
ind.Hobs$pop <- "NA"
ind.Hobs$mlg <- "NA"
AllPops.gc[1]@mlg

for(i in 1:114) {
  ind.Hobs[i,5] <- as.character(AllPops.gc[i]@strata$pop)
  ind.Hobs[i,4] <- as.character(AllPops.gc[i]@strata$ms)
  #ind.Hobs[i,3] <- mean(summary(AllPops.gc[i])$Hobs, na.rm=TRUE)
  #ind.Hobs[i,2] <- ploidy(AllPops.gc[i])
  #ind.Hobs[i,1] <- indNames(AllPops.gc[i])
  print(i)
}

ind.Hobs$ms <- factor(ind.Hobs$ms, levels=c("S", "A"))

mean.Hobs <- ind.Hobs %>%
  group_by(ms) %>%
  summarise(mean=mean(Hobs), sd=sd(Hobs))

gg.Hobs <-ggplot()+
  geom_boxplot(data=ind.Hobs, aes(y=Hobs, x=ms), alpha=0, width=0.3, outlier.shape=NA)+
  geom_jitter(data=ind.Hobs, aes(x=ms, y=Hobs, fill=ms), alpha=1, shape=21, width=0.05, size=2)+
  geom_point(data=mean.Hobs, aes(x=as.numeric(as.factor(ms))+0.4, y=mean, fill=ms), shape=21, size=3.5)+
  geom_errorbar(data=mean.Hobs, aes(ymax=mean+sd, ymin=mean-sd, x=as.numeric(as.factor(ms))+0.4), width=0)+
  scale_fill_manual(values=c("red3","slateblue4"), name="Mating system", labels=c("Sexual", "Apomictic"))+
  labs(x="Mating system", y="Mean observed heterozygosity")+
  theme_classic()+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))

png("ind_het_ms.png", height=7, width=8, res=300, units="in")
gg.Hobs
dev.off()

kruskal.test(data=ind.Hobs, Hobs~ms)
