##### NJ TREE #####
hookeri.nj1000 <- aboot(AllPops.gc, dist = provesti.dist, sample = 1000, tree = "nj", cutoff = 50, quiet = TRUE)

saveRDS(hookeri.nj1000, here("data", "hookeri_nj1000.rds"))

hookeri.nj1000$node.label2 <- as.numeric(hookeri.nj1000$node.label)
hookeri.nj1000$node.label2[hookeri.nj1000$node.label<90] <- 1
hookeri.nj1000$node.label2[hookeri.nj1000$node.label>90] <- 2
hookeri.nj1000$node.label2[is.na(hookeri.nj1000$node.label2)] <- 1

tree.pch <- group.cols3[AllPops.gc@other$group]
tree.pch[1:41] <- 24
tree.pch[42:114] <- 21
tree.pch <- as.numeric(tree.pch)

tree.legend.pch <- c(24, 24, 21, 21, 21, 24, 21, 21, 21, 21, 21, 21, 21, 24)

png("njtree_all_v2.png", height=7, width=7, res=300, units="in")
plot.phylo(hookeri.nj1000, type="unrooted", cex=0.6, lab4ut = "axial", font=2, show.tip.label = FALSE, no.margin = TRUE)
tiplabels(pch=tree.pch, col="black", bg=group.cols3[AllPops.gc@other$group], cex=1.3)
add.scale.bar()
legend("bottomright", legend=levels(AllPops.gc@other$group), cex=0.8, pch=tree.legend.pch, col="black", pt.bg=group.cols3, pt.cex=1.7, inset=0.01)
dev.off()
