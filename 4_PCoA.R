compositional <- microbiome::transform(cricket, transform = "compositional")

bray <- ordinate(compositional, method = "PCoA", distance = "bray")

B <- plot_ordination(physeq = compositional,
                ordination = bray,
                type = "samples",
                color = "type") + 
  theme_bw() + scale_color_manual(values=c("#7FB3D5","#E9C46A","#264653", "#F4A261","#2A9D8F")) +
  geom_point(size=3) + ggtitle("B. Bray-Curtis distance matrix by type") +
  theme(legend.position="bottom")
B + stat_ellipse(type = "t", linetype = 2)


S <- plot_ordination(physeq = compositional,
                     ordination = bray,
                     type = "samples",
                     shape = "HighORLow", color = "HighORLow") + 
  theme_bw() + scale_color_manual(values=c("#7FB3D5","#E9C46A"))+
  scale_shape_manual(values = c(19, 1)) +
  geom_point(size=3) + ggtitle("A. Bray-Curtis distance matrix by sperm viability")+
  theme(legend.position="bottom")
S + stat_ellipse(type = "t", linetype = 2)

grid.arrange(S,B, ncol=2)


rm(B)
rm(S)
rm(bray)
rm(compositional)
