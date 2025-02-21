##### Observed heterozygosity #####

ind.Hobs <- data.frame(matrix(NA, nrow=114, ncol=3))
colnames(ind.Hobs) <- c("ind", "ploidy", "Hobs")
ind.Hobs$ms <- "NA"
ind.Hobs$pop <- "NA"
ind.Hobs$mlg <- "NA"
AllPops.gc[1]@mlg

for(i in 1:114) {
  ind.Hobs[i,5] <- as.character(AllPops.gc[i]@strata$pop)
  ind.Hobs[i,4] <- as.character(AllPops.gc[i]@strata$ms)
  ind.Hobs[i,3] <- mean(summary(AllPops.gc[i])$Hobs, na.rm=TRUE)
  ind.Hobs[i,2] <- ploidy(AllPops.gc[i])
  ind.Hobs[i,1] <- indNames(AllPops.gc[i])
  print(i)
}

ind.Hobs$ms <- factor(ind.Hobs$ms, levels=c("S", "A"))
ind.Hobs$group <- AllPops.gc@other$group

# Second Plot (p4)
gg.Hobs.bar <- ggplot(ind.Hobs, aes(x = ind, y = Hobs)) +
  geom_bar(stat = "identity", color = "white", fill = "grey40") +  # Ensure bars are grey
  facet_grid(. ~ pop, scales = "free_x", space = "free") +
  #scale_y_reverse() +  # Flip the y-axis without changing values
  theme_bw() +
  ylab("Observed Heterozygosity (Hobs)") +  # Keep y-axis label
  xlab(NULL) +  # Remove x-axis title
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    panel.spacing.x = unit(0.1, "lines"),
    strip.text.x = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet background
    panel.border = element_blank(),  # Remove black facet borders
    panel.grid = element_blank()  # Remove all grid lines
  )

gg.Hobs.bar

mean.Hobs <- ind.Hobs %>%
  group_by(ms) %>%
  summarise(mean=mean(Hobs), sd=sd(Hobs))

gg.Hobs.boxplot <- ggplot()+
  geom_boxplot(data=ind.Hobs, aes(y=Hobs, x=ms), alpha=0, width=0.2, outlier.shape=NA)+
  geom_jitter(data=ind.Hobs, aes(x=ms, y=Hobs, fill=group), alpha=1, shape=21, width=0.05, size=2.5)+
  scale_fill_manual(values=group.cols3, name="Mating system")+
  #guides(fill="none")+
  new_scale_fill()+
  geom_errorbar(data=mean.Hobs, aes(ymax=mean+sd, ymin=mean-sd, x=as.numeric(as.factor(ms))+0.4), width=0)+
  geom_point(data=mean.Hobs, aes(x=as.numeric(as.factor(ms))+0.4, y=mean, fill=ms), shape=21, size=3.5, show.legend = FALSE)+
  scale_fill_manual(values=c("white", "white"))+
  labs(x="Mating system", y="Mean observed heterozygosity")+
  theme_bw()+
  scale_x_discrete(labels=c("S" = "Sexual", "A" = "Apomictic"))+
  theme(legend.title=element_blank(), 
        legend.text = element_text(size=10), 
        legend.key.size = unit(0.3, 'cm'),
        legend.position=c(0.9, 0.2), 
        legend.background = element_blank(), 
        legend.box.background = element_rect(colour = "black"))


png("ind_het_ms.png", height=7, width=8, res=300, units="in")
gg.Hobs
dev.off()