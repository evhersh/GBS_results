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

# boxplot for dips (should do one box for all loci?)

mean(summary(dips.gc)$Hobs) # mean = 0.214936
mean(summary(trips.gc)$Hobs) # mean = 0.4256852
mean(tets.summary$Hobs) # mean = 0.4174576

dips.summary <- summary(dips.gc)
trips.summary <- summary(trips.gc)
tets.summary <- summary(tets.gc)

summary(apos.gc)

dips.Hobs <- data.frame(group="diploid", value=dips.summary$Hobs)
trips.Hobs <- data.frame(group="triploid", value=trips.summary$Hobs)
tets.Hobs <- data.frame(group="tetraploid", value=tets.summary$Hobs)

all.Hobs <- rbind(dips.Hobs, trips.Hobs)

ggplot(all.Hobs, aes(x=group, y=value, fill=group))+
  geom_boxplot(color="black")+
  scale_fill_manual(values=c("coral3", "cornflowerblue"), name="Mating System", labels=c("Sexual", "Apomictic"))+
  stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=18, size=3,show_guide = FALSE) +
  theme_bw()+
  labs(x="Mating System", y="Observed heterozygosity")
  
apos.subYK.gc <- setPop(apos.subYK.gc, formula=~mlg)
# Hobs per pop - dips
n.pop.trips <- seppop(apos.subYK.gc)

mlg110.sum <- summary(n.pop.trips$`110`)
mlg70.sum <- summary(n.pop.trips$`70`)
mlg55.sum <- summary(n.pop.trips$`55`)
mlg103.sum <- summary(n.pop.trips$`103`)
mlg91.sum <- summary(n.pop.trips$`91`)
mlg112.sum <- summary(n.pop.trips$`112`)
mlg105.sum <- summary(n.pop.trips$`105`)
mlg104.sum <- summary(n.pop.trips$`104`)
mlg5.sum <- summary(n.pop.trips$`5`)

mlg110.hobs <- data.frame(group="110", value=mlg110.sum$Hobs)
mlg70.hobs <- data.frame(group="70", value=mlg70.sum$Hobs)
mlg55.hobs <- data.frame(group="55", value=mlg55.sum$Hobs)
mlg103.hobs <- data.frame(group="103", value=mlg103.sum$Hobs)
mlg91.hobs <- data.frame(group="91", value=mlg91.sum$Hobs)
mlg112.hobs <- data.frame(group="112", value=mlg112.sum$Hobs)
mlg105.hobs <- data.frame(group="105", value=mlg105.sum$Hobs)
mlg104.hobs <- data.frame(group="104", value=mlg70.sum$Hobs)
mlg5.hobs <- data.frame(group="5", value=mlg70.sum$Hobs)

mlg.Hobs <- rbind(mlg110.hobs, mlg70.hobs, mlg55.hobs, mlg103.hobs, mlg91.hobs, mlg112.hobs, mlg105.hobs, mlg104.hobs, mlg5.hobs)

ggplot(mlg.Hobs, aes(x=group, y=value, color=group))+
  geom_boxplot()+
  theme_bw()

#pop.hobs.trips <- do.call("c", lapply(n.pop.trips, function(x) summary(x)$Hobs)) 



#mean.hobs[is.nan(mean.hobs)] <- NA  # probably not necessary

barplot(pop.hobs.dips)
