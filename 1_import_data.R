#set working directory

# load useful function
KillZeroRCs <- function(x) {
  x[ which( rowSums(x) != 0) , ] -> x
  x[ , which( colSums(x) != 0) ] -> x
  return(x)
}

# read in OTU table
ASV <- read.table(
  "asv.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# remove zero-sum row/columns
ASV <- KillZeroRCs(ASV)

# read in tax table
# fill logical indicates that rows have unequal lengths due to blank fields
TAXA <- read.table(
  "taxonomy.txt",
  sep = "\t",
  fill = TRUE,
  row.names = 1
)

# add levels of taxonomy to tax table
colnames(TAXA) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")

# read in METAdata
META <- read.table(
  "metadataR.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

META$treat <- paste(META$type_ID, META$HighORLow, sep = "_")

# convert OTU and tax tables to matrices for phyloseq
ASV <- as.matrix(ASV)
TAXA <- as.matrix(TAXA)

# load libraries
library(ape)
library(phyloseq)
library(microbiome)

# convert tree file to phyloseq object
TREE <- read_tree("tree.nwk")

# combine OTU, tax and METAdata files into a phyloseq object
phy <- phyloseq(
  otu_table(ASV,taxa_are_rows = T),
  tax_table(TAXA),
  sample_data(META)
)

# merge phyloseq objects
phy <- merge_phyloseq(phy, TREE)

# Delete the matrices as they are not used in further analyses
rm(ASV)
rm(TAXA)
rm(META)
rm(TREE)

# remove zero sum ASV
phy <- prune_taxa((taxa_sums(phy) > 0), phy) #1928 ASVs in 156 samples

#adjust tree for ASVs that have been removed
phy_tree(phy) <- root(phy_tree(phy),sample(taxa_names(phy),1), resolve.root = TRUE)

rm(KillZeroRCs)

