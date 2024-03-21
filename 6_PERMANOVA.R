##############################
##        PERMANOVA         ##
##############################

library(vegan)
library(pairwiseAdonis)

##A PERMANOVA was run on the data set with a Bray-Curtis distance matrix

compositional <- microbiome::transform(cricket, transform = "compositional")

#Generate Unifrac distance matrix
bray_dist_matrix <- phyloseq::distance(compositional, method = "bray")

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis

dispr.bray <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(compositional)$type)
dispr.bray
plot(dispr.bray)
anova(dispr.bray)#p=4.77e-09
#reject the assumption of homogeneity of dispersion by type

dispr.bray <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(compositional)$HighORLow)
dispr.bray
plot(dispr.bray)
anova(dispr.bray)#p=0.2708
#fail to reject the assumption of homogeneity of dispersion by sperm viability

#Because the bray matrix met the assumption of homogeneity for sperm viability only, we continue with that distance matrix but 
#add a permutation test for multivariate dispersion

#ADONIS test
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(compositional)$type*phyloseq::sample_data(compositional)$HighORLow, strata = phyloseq::sample_data(compositional)$maleID, perm=9999, permutest="dispersion") 

#Number of permutations: 9999

#vegan::adonis2(formula = bray_dist_matrix ~ phyloseq::sample_data(compositional)$type * phyloseq::sample_data(compositional)$HighORLow, permutations = 9999, strata = phyloseq::sample_data(compositional)$maleID, permutest = "dispersion")
#                                                                                          Df SumOfSqs      R2      F Pr(>F)    
#phyloseq::sample_data(compositional)$type                                                  4    7.606 0.17535 7.0585 0.0001 ***
#phyloseq::sample_data(compositional)$HighORLow                                             1    0.384 0.00884 1.4238 0.0006 ***
#phyloseq::sample_data(compositional)$type:phyloseq::sample_data(compositional)$HighORLow   4    0.904 0.02085 0.8391 0.5930    
#Residual                                                                                 128   34.481 0.79496                  
#Total                                                                                    137   43.375 1.00000                                        

###Running a pairwise adonis

#pull out metadata from phyloseq object
meta <- meta(compositional)

results <- pairwise.adonis(bray_dist_matrix, factors = meta$type, perm = 9999, p.adjust.m = "holm")
#                                 pairs Df   SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#1    spermatophore vs accessory_gland  1 0.009354203  4.801096 0.08603945   0.001      0.010   *
#2    spermatophore vs digestive_tract  1 0.035503736 22.921889 0.31008256   0.001      0.010   *
#3    spermatophore vs seminal_vesicle  1 0.004867090  2.461629 0.04692247   0.006      0.024   .
#4             spermatophore vs testes  1 0.007325000  3.766601 0.06635241   0.001      0.010   *
#5  accessory_gland vs digestive_tract  1 0.026547484 19.738442 0.26768184   0.001      0.010   *
#6  accessory_gland vs seminal_vesicle  1 0.003091504  1.771479 0.03234309   0.069      0.207    
#7           accessory_gland vs testes  1 0.001883092  1.090468 0.01910071   0.350      0.402    
#8  digestive_tract vs seminal_vesicle  1 0.031020745 22.796271 0.30075716   0.001      0.010   *
#9           digestive_tract vs testes  1 0.022557534 16.548894 0.22810678   0.001      0.010   *
#10          seminal_vesicle vs testes  1 0.002282172  1.304805 0.02317395   0.201      0.402    


rm(results)
rm(compositional)
rm(dispr.bray)
rm(bray_dist_matrix)
rm(meta)

