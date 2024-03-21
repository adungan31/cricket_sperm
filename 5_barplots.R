library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("microbiome")

##########################################
##                Bar plots             ##
##########################################

#Wolbachia ASVs by Tissue Type

merge <- merge_samples(cricket, "treat")
comp <- transform_sample_counts(merge, function(x) 100 * x/sum(x))
Wol <- subset_taxa(comp, Genus=="Wolbachia")

#reorder levels of type
sample_data(Wol)$order <- factor(sample_names(Wol))
sample_data(Wol)$order <- factor(sample_data(Wol)$order, levels=c("D_Low","D_High","AG_Low","AG_High","SV_Low","SV_High","T_Low","T_High","S_Low","S_High"))
meta <- meta(Wol)

p = plot_bar(Wol, fill = "OTU", x="order")+theme(axis.title.x = element_blank())
print(p)

 ##Genus level##
# To represent at the Genus level (instead of ASV level)
Genus <- tax_glom(merge,taxrank = "Genus")


# Transform counts in relative abundance and select most abundant families
Genus <- transform_sample_counts(Genus, function(x) 100 * x/sum(x))
genus <- psmelt(Genus)
genus$Genus <- as.character(genus$Genus)

#rename Genera with <4% abundance
genus$Genus[genus$Abundance < 4]<- "< 4% Abundance"

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$Genus)))

#write.csv(genus, "genus.csv", row.names=FALSE)
genus <- read.csv("genus.csv")


#reorder levels of type
genus$Tissue <- factor(genus$Tissue, levels=c("D","AG","SV","T","S"))
genus$Sperm <- factor(genus$Sperm, levels = c("Low","High"))


A <- ggplot(genus, aes(x = Tissue, y = Tissue_Abun, fill = reorder(Genus, Tissue_Abun))) +  
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Genera \n") +  
  scale_fill_manual(values = c( "#7FB3D5","#E9C46A","#000000", "#2A9D8F", "#264653", "#F4A261","#D1E7F1","#CCCCCC", "#7ECCA4", "#1F78B4","#B2DF8A","#555555","#33A02C","#FB9A99", "#6A3D9A","#FFFFFF", "#B7B5E4"))
print(A)

B <- ggplot(genus, aes(x = Sperm, y = Sperm_Abun, fill = reorder(Genus, Sperm_Abun))) +  
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Genera (%) \n") +  
  scale_fill_manual(values = c( "#7FB3D5","#E9C46A","#000000", "#2A9D8F", "#264653", "#F4A261","#D1E7F1","#CCCCCC", "#7ECCA4", "#1F78B4","#B2DF8A","#555555","#33A02C","#FB9A99", "#6A3D9A","#FFFFFF", "#B7B5E4"))
print(B)

##Phylum level##
# To represent at the Genus level (instead of ASV level)
P <- tax_glom(merge,taxrank = "Phylum")


# Transform counts in relative abundance and select most abundant families
P <- transform_sample_counts(P, function(x) 100 * x/sum(x))
phylum <- psmelt(P)
phylum$Phylum <- as.character(phylum$Phylum)

#rename Genera with <1% abundance
phylum$Phylum[phylum$Abundance < 1]<- "< 1% Abundance"

#How many levels in Phylum
HowMany <- length(levels(as.factor(phylum$Phylum)))

#write.csv(phylum, "phylum.csv", row.names=FALSE)
phylum <- read.csv("phylum.csv")


#reorder levels of type
phylum$Tissue <- factor(phylum$Tissue, levels=c("D","AG","SV","T","S"))
phylum$Sperm <- factor(phylum$Sperm, levels = c("Low","High"))


C <- ggplot(phylum, aes(x = Tissue, y = Tissue_Abun, fill = reorder(Phylum, Tissue_Abun))) +  
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Phyla \n") +  
  scale_fill_manual(values = c( "#7FB3D5","#E9C46A","#000000", "#2A9D8F", "#264653", "#F4A261","#D1E7F1","#CCCCCC"))
print(C)

D <- ggplot(phylum, aes(x = Sperm, y = Sperm_Abun, fill = reorder(Phylum, Sperm_Abun))) +  
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_blank()) +
  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Phyla (%) \n") +  
  scale_fill_manual(values = c( "#7FB3D5","#E9C46A","#000000", "#2A9D8F", "#264653", "#F4A261","#D1E7F1","#CCCCCC"))


print(D)

rm(comp)
rm(A)
rm(B)
rm(C)
rm(D)
rm(genus)
rm(Genus)
rm(merge)
rm(meta)
rm(P)
rm(p)
rm(phylum)
rm(Wol)
rm(HowMany)
