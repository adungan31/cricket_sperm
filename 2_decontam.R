library(decontam)
library(microbiome)

# identify contaminants first from PCR negatives
consList <- isContaminant(seqtab = phy, neg = "PCR_control", method = "prevalence")

# pull out the names of contaminants
cons <- rownames(consList)[consList$contaminant=="TRUE"]
cons <- as.character(cons) #6 contaminants identified from PCR_controls

# - - - - - - - - - - - - - - - - - - - - - - - - - #

# to get info on the contaminants, uncomment the following code to
# run it ON THE FILE WITH THE CONTAMINANT ASVs IN-SITU,
# then combine the consPer.csv and taxonomy.csv file data

# subset the non neg-control anemone samples
#vvv <- subset_samples(bacteria, PCR_control == "FALSE")
# merge the samples
#yyy <- merge_samples(vvv, "PCR_control", fun = sum)
# transform counts to percentages
#yyy <- transform_sample_counts(yyy, function(x) 100 * x/sum(x))
# extract the cons percentage data
#zzz <- prune_taxa(x = yyy, taxa = cons)
# write otu table to dataframe
#xxx <- data.frame(t(zzz@otu_table))
# write xxx to csv
#write.csv(x = xxx, row.names = TRUE, file = "consPCRcontrol.csv")
# subset the contaminant ASVs
#bacteriaCons <- prune_taxa(bacteria, taxa = cons)
# write the contaminants to a file for reference
#contaminants <-bacteriaCons@tax_table
#contaxa <- contaminants@.Data
#write.csv(contaxa, "contaxaPCR_control.csv")

# 6 contaminant ASVs in PCR negatives
# total contamination in the samples = 1.04%


# - - - - - - - - - - - - - - - - - - - - - - - - - #

# remove the contaminants from the main bacteria phyloseq file
bacteria <- remove_taxa(phy, taxa = cons)


#Remove PCR_controls
bacteria <- subset_samples(bacteria, PCR_control == "FALSE") 


# identify contaminants from extract blanks
consList2 <- isContaminant(seqtab = bacteria, neg = "Extract_control", method = "prevalence")

# pull out the names of contaminants
cons2 <- rownames(consList2)[consList2$contaminant=="TRUE"]
cons2 <- as.character(cons2) #6 contaminants identified from extraction blanks

# - - - - - - - - - - - - - - - - - - - - - - - - - #

# to get info on the contaminants, uncomment the following code to
# run it ON THE FILE WITH THE CONTAMINANT ASVs IN-SITU,
# then combine the consPer.csv and taxonomy.csv file data

# subset the non neg-control anemone samples
#vvv2 <- subset_samples(bacteria, Extract_control == "FALSE")
# merge the samples
#yyy2 <- merge_samples(vvv2, group = "Extract_control", fun = sum)
# transform counts to percentages
#yyy2 <- transform_sample_counts(yyy2, function(x) 100 * x/sum(x))
# extract the cons percentage data
#zzz2 <- prune_taxa(x = yyy2, taxa = cons2)
# write otu table to dataframe
#xxx2 <- data.frame(t(zzz2@otu_table))
# write xxx to csv
#write.csv(x = xxx2, row.names = TRUE, file = "consExtraction.csv")
# subset the contaminant ASVs
#bacteriaCons2 <- prune_taxa(bacteria, taxa = cons2)
# write the contaminants to a file for reference
#contaminants2 <-bacteriaCons2@tax_table
#contaxa2 <- contaminants2@.Data
#write.csv(contaxa2, "contaxaExtraction.csv")


# 6 contaminant ASVs in Extraction blanks
# total contamination in the samples = 0.02%

# remove the contaminants from the main bacteria phyloseq file
bacteria <- remove_taxa(bacteria, taxa = cons2)

# Remove extraction controls from bacteria 
bacteria <- subset_samples(bacteria, Extract_control=="FALSE")

#subset mock community sample and remove them from  bacteria phyloseq object
cricket <- remove_samples("ZymoMC.F18" , bacteria)
cricket <- prune_taxa((taxa_sums(cricket) > 0), cricket)# 1892 ASVs in 145 samples

#view reads by sample
sort(sample_sums(cricket)) 
#12 samples with < 1000 reads
#C128.F17 
#C128SV.F17   
#C19D.F12    
#C86.F17 
#C28SV.F12  
#C19AG.F11  
#C82AG.F13  
#C73AG.F14  
#C75AG.F13  
#C73SV.F13 
#C127T.F13  
#C75SV.F13

#remove samples with less than 400 reads
cricket <- prune_samples((sample_sums(cricket) > 400), cricket)
cricket <- prune_taxa((taxa_sums(cricket) > 0), cricket)# 1887 ASVs in 138 samples
phy_tree(cricket) <- root(phy_tree(cricket),sample(taxa_names(cricket),1), resolve.root = TRUE)
# - - - - - - - - - - - - - - - - - - - - - - - - - #

rm(consList)
rm(phy)
rm(bacteria)
rm(cons)
rm(consList2)
rm(cons2)
