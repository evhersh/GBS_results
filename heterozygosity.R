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
  
# apos.subYK.gc <- setPop(apos.subYK.gc, formula=~mlg)
# # Hobs per pop - dips
# n.pop.trips <- seppop(apos.subYK.gc)
# 
# mlg110.sum <- summary(n.pop.trips$`110`)
# mlg70.sum <- summary(n.pop.trips$`70`)
# mlg55.sum <- summary(n.pop.trips$`55`)
# mlg103.sum <- summary(n.pop.trips$`103`)
# mlg91.sum <- summary(n.pop.trips$`91`)
# mlg112.sum <- summary(n.pop.trips$`112`)
# mlg105.sum <- summary(n.pop.trips$`105`)
# mlg104.sum <- summary(n.pop.trips$`104`)
# mlg5.sum <- summary(n.pop.trips$`5`)
# 
# mlg110.hobs <- data.frame(group="110", value=mlg110.sum$Hobs)
# mlg70.hobs <- data.frame(group="70", value=mlg70.sum$Hobs)
# mlg55.hobs <- data.frame(group="55", value=mlg55.sum$Hobs)
# mlg103.hobs <- data.frame(group="103", value=mlg103.sum$Hobs)
# mlg91.hobs <- data.frame(group="91", value=mlg91.sum$Hobs)
# mlg112.hobs <- data.frame(group="112", value=mlg112.sum$Hobs)
# mlg105.hobs <- data.frame(group="105", value=mlg105.sum$Hobs)
# mlg104.hobs <- data.frame(group="104", value=mlg70.sum$Hobs)
# mlg5.hobs <- data.frame(group="5", value=mlg70.sum$Hobs)
# 
# mlg.Hobs <- rbind(mlg110.hobs, mlg70.hobs, mlg55.hobs, mlg103.hobs, mlg91.hobs, mlg112.hobs, mlg105.hobs, mlg104.hobs, mlg5.hobs)
# 
# ggplot(mlg.Hobs, aes(x=group, y=value, color=group))+
#   geom_boxplot()+
#   theme_bw()

#pop.hobs.trips <- do.call("c", lapply(n.pop.trips, function(x) summary(x)$Hobs)) 



#mean.hobs[is.nan(mean.hobs)] <- NA  # probably not necessary

barplot(pop.hobs.dips)

# MLG heterozygosity

mlg57.sum <- summary(mlg57.gc)
mlg103.sum <- summary(mlg103.gc)
mlg7.sum <- summary(mlg7.gc)
mlg12.sum <- summary(mlg12.gc)
mlg52.sum <- summary(mlg52.gc)
mlg91.sum <- summary(mlg91.gc)
mlg104.sum <- summary(mlg104.gc)
mlg110.sum <- summary(mlg110.gc)
mlg114.sum <- summary(mlg114.gc)
mlg105.sum <- summary(mlg105.gc)

mlg57.hobs <- data.frame(group="MT_big", value=mlg57.sum$Hobs) # mean = 0.4422016
mlg103.hobs <- data.frame(group="Lrmy_big", value=mlg103.sum$Hobs) # 0.4281033
mlg7.hobs <- data.frame(group="SK", value=mlg7.sum$Hobs) # 0.3864135
mlg12.hobs <- data.frame(group="YK", value=mlg12.sum$Hobs) # 0.4174576
mlg52.hobs <- data.frame(group="BC", value=mlg52.sum$Hobs) # 0.5080658
mlg91.hobs <- data.frame(group="Lrmy_small", value=mlg91.sum$Hobs) # 0.4476231
mlg104.hobs <- data.frame(group="Dan's", value=mlg104.sum$Hobs) # 0.3443815
mlg110.hobs <- data.frame(group="ND", value=mlg110.sum$Hobs)
mlg114.hobs <- data.frame(group="WY_4", value=mlg114.sum$Hobs)
mlg105.hobs <- data.frame(group="WY_1", value=mlg105.sum$Hobs)

mlg.Hobs <- rbind(mlg57.hobs,mlg103.hobs,mlg7.hobs,mlg12.hobs,mlg52.hobs,mlg91.hobs,mlg104.hobs,mlg110.hobs,mlg114.hobs,mlg105.hobs)

mlg.Hobs$group <- factor(mlg.Hobs$group, levels=c("Dan's", "Lrmy_small", "Lrmy_big", "WY_1", "WY_4", "MT_big", "ND", "SK", "BC", "YK"))


ggplot(mlg.Hobs, aes(x=group, y=value, fill=group))+
  stat_summary(fun.y=mean, geom="point", na.rm=TRUE, size=4, pch=21, colour="black")+
  stat_summary(fun.data=mean_se, geom="errorbar", na.rm=TRUE)
