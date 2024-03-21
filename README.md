# cricket_sperm
QIIME2 and R code for "Microbial biomarkers as indicators of sperm quality in an insect"

All raw data from the Illumina MiSeq run are available under NCBI BioProject ID PRJNA971357.

Files provided
Cricket_metabarcoding_QIIME2_code.pdf: This outlines the QIIME2 code for processing all raw data files, generating the precursor files for data analysis in R. Users unfamiliar with QIIME2 can follow along with a [tutorial I designed in partnership with Melbourne Bioinformatics](https://www.melbournebioinformatics.org.au/tutorials/tutorials/qiime2/qiime2/). 

Data files used for analysis in R:

a. metadataR.txt

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


Should you have any questions or find a mistake, please contact Ashley Dungan at adungan31@gmail.com
You are welcome to use this code for any of your research - please reference this GitHub page if you do.
