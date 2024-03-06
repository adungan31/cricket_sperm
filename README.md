# cricket_sperm
QIIME2 and R code for "Fecal bacterial communities of the platypus (Ornithorhynchus anatinus) reflect captivity status – implications for conservation and management"

All raw data from the Illumina MiSeq run are available under NCBI BioProject ID PRJNA971672.

Files provided
Platypus_metarbarcoding.pdf: This outlines the QIIME2 code for processing all raw data files, generating the precursor files for data analysis in R. Users unfamiliar with QIIME2 can follow along with a tutorial I designed in partnership with Melbourne Bioinformatics

Data files used for analysis in R:

a. metadata.txt

b. taxonomy.txt

c. asv.txt

d. tree.nwk

R code files:

a. 1_import_data.R

b. 2_decontam.R

c. 3_diversity_metrics.R

d. 3_diversity_metrics_stats.R

e. 4_PCoA.R

f. 5_barplots.R

g. 6_PERMANOVA.R

h. 10_most_abundant_ASVs.R

i. 14_microbiomeMarker.R

j. 16_UpSetR.R

Should you have any questions or find a mistake, please contact Ashley Dungan at adungan31@gmail.com
You are welcome to use this code for any of your research - please reference this GitHub page if you do.
