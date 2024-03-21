library(plotrix)
library(microbiome)
library(gridExtra)

# # check read depth & select rarefaction level
sort(sample_sums(cricket))

# rarefy
cricket_rare<- rarefy_even_depth(cricket, sample.size = 2133, rngseed = 1)
#21 samples and 245 ASVs were removed

# extract metadata from subsetted file
divMeta <- meta(cricket_rare)

# add diversity index data to metadata file
divMeta <- cbind(divMeta, estimate_richness(cricket_rare))
write.csv(divMeta, "adiv.csv")


##BoxPlot###

##by tissue type
#reorder levels of type
divMeta$type <- factor(divMeta$type, levels=c("digestive_tract","accessory_gland","seminal_vesicle","testes","spermatophore"))


cricket.obs <- ggplot(divMeta, aes(x=type, y=Observed, fill=type)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black")+
  scale_fill_manual(values=c("#E9C46A", "#7FB3D5", "#264653", "#2A9D8F", "#F4A261"))+ lims(y=c(0,160))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+ theme(axis.title.x = element_blank())+theme(axis.text.x = element_blank())+
  ggtitle("A. Observed ASVs")
print(cricket.obs)

cricket.sim <- ggplot(divMeta, aes(x=type, y=Simpson, fill=type)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black")+
  scale_fill_manual(values=c("#E9C46A", "#7FB3D5", "#264653", "#2A9D8F", "#F4A261"))+ lims(y=c(0,1))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+theme(axis.title.x = element_blank())+theme(axis.text.x = element_blank())+
  ggtitle("C. Simpson")
print(cricket.sim)

cricket.shan <- ggplot(divMeta, aes(x=type, y=Shannon, fill=type)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black")+
  scale_fill_manual(values=c("#E9C46A", "#7FB3D5", "#264653", "#2A9D8F", "#F4A261"))+ lims(y=c(0,4.5))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+theme(axis.title.x = element_blank())+theme(axis.text.x = element_blank())+
  ggtitle("D. Shannon")
print(cricket.shan)

cricket.Chao1 <- ggplot(divMeta, aes(x=type, y=Chao1, fill=type)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black")+
  scale_fill_manual(values=c("#E9C46A", "#7FB3D5", "#264653", "#2A9D8F", "#F4A261"))+ lims(y=c(0,175))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+theme(axis.title.x = element_blank())+theme(axis.text.x = element_blank())+
  ggtitle("B. Chao1")
print(cricket.Chao1)


grid.arrange(cricket.obs,cricket.Chao1,cricket.sim,cricket.shan, nrow=2, ncol=2)

##by spermviability
#reorder levels of type
divMeta$HighORLow <- factor(divMeta$HighORLow, levels=c("Low","High"))


cricket.obs <- ggplot(divMeta, aes(x=HighORLow, y=Observed, fill=HighORLow)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black")+
  lims(y=c(0,180))+theme(axis.title.y = element_blank())+
  scale_fill_manual(values=c("#E9C46A", "#7FB3D5"))+
  theme(axis.text.y = element_text(size = 14))+ theme(axis.title.x = element_blank())+theme(axis.text.x = element_blank())+
  ggtitle("A. Observed ASVs")
print(cricket.obs)

cricket.sim <- ggplot(divMeta, aes(x=HighORLow, y=Simpson, fill=HighORLow)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black")+
  scale_fill_manual(values=c("#E9C46A", "#7FB3D5"))+ lims(y=c(0,1))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+theme(axis.title.x = element_blank())+theme(axis.text.x = element_blank())+
  ggtitle("C. Simpson")
print(cricket.sim)

cricket.shan <- ggplot(divMeta, aes(x=HighORLow, y=Shannon, fill=HighORLow)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black")+
  scale_fill_manual(values=c("#E9C46A", "#7FB3D5"))+ lims(y=c(0,4.5))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+theme(axis.title.x = element_blank())+theme(axis.text.x = element_blank())+
  ggtitle("D. Shannon")
print(cricket.shan)

cricket.Chao1 <- ggplot(divMeta, aes(x=HighORLow, y=Chao1, fill=HighORLow)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=4, color="black")+
  scale_fill_manual(values=c("#E9C46A", "#7FB3D5"))+ lims(y=c(0,200))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+theme(axis.title.x = element_blank())+theme(axis.text.x = element_blank())+
  ggtitle("B. Chao1")
print(cricket.Chao1)


grid.arrange(cricket.obs,cricket.Chao1,cricket.sim,cricket.shan, nrow=2, ncol=2)


rm(cricket_rare)
rm(cricket.obs)
rm(cricket.shan)
rm(cricket.Chao1)
rm(cricket.sim)

