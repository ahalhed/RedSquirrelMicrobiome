# I was working with this R file with the files on SHARCNET, but
# I decided to download all to my local device to become more
# familiar with the workings of QIIME2 artifacts in R.
```{bash}
module load miniconda3
conda activate qiime2-2019.10
# so there must apparently be zeroes somewhere in the table
# take them out to make qiime2R happy
qiime feature-table filter-samples \
  --i-table filtered-table.qza \
  --p-min-features 1 \
  --o-filtered-table filtered-table-no0.qza

# allocation for testing (64G needed for qza_to_phyloseq)
salloc --time=1:0:0 --mem-per-cpu 64G --account=def-cottenie

# get/activate R within allocation
module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.0
R
```
# set working directory
setwd("home/ahalhed/red-squirrel-w2020/R-env")
# load workspace image
load(file = "RedSquirrel.RData")
# attach required packages (install them on login node)
library(tidyverse)
# devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(phyloseq)

# bring QIIME2 artifacts into R (needs BIG memory - 64G worked)
# data into phyloseq
ps <- qza_to_phyloseq(features = "../filtered-table-no0.qza",
                      tree = "../trees/rooted_tree.qza",
                      taxonomy = "../taxonomy/GG-taxonomy.qza",
                      metadata = "../input/RS_meta.tsv")
# distance matrices
dm_phylo_bray <- read_qza("../core-metrics-phylogenetic/bray_curtis_distance_matrix.qza")
dm_phylo_jaccard <- read_qza("../core-metrics-phylogenetic/jaccard_distance_matrix.qza")
dm_phylo_wu <- read_qza("../core-metrics-phylogenetic/weighted_unifrac_distance_matrix.qza")
dm_phylo_uu <- read_qza("../core-metrics-phylogenetic/unweighted_unifrac_distance_matrix.qza")

# explore the phyloseq object
ntaxa(ps) # 541922
nsamples(ps) # 909
sample_names(ps)[1:5]
#[1] "B2003.11254.F.1.KL.2010" "B2571.10374.F.3.KL.2010"
#[3] "B2197.10318.F.4.KL.2010" "B2067.11272.F.1.KL.2010"
#[5] "P3169.10686.F.1.KL.2008"
rank_names(ps) #"Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
sample_variables(ps)
# [1] "Grid"           "Location.X"     "Location.Y"     "Sex"           
# [5] "Age"            "Month"          "Season"         "Year"          
# [9] "Squirrel.ID"    "SireID"         "DamID"          "CollectionDate"
#[13] "FoodSupplement" "BirthYear"      "Location"       "Date"    
taxa_names(ps)[1:10] #this isn't interesting to look at
# access different components of the phyloseq object
otu_table(ps)
tax_table(ps)
phy_tree(ps) # Phylogenetic tree with 541922 tips and 541921 internal nodes.


# multivariate ordination for the full data set
ord_ps_NMDS_bray <- ordinate(ps, "NMDS", dm_phylo_bray$data)
ord_ps_RDA_bray <- ordinate(ps, "RDA", dm_phylo_bray$data)
ord_ps_DCA_bray <- ordinate(ps, "DCA", dm_phylo_bray$data)

# plot the ordinations
# NDMS bray aes(color = sample_data(resampled.ps)$Squirrel.ID)
pdf("./plots/NDMSbray.pdf")
plot_ordination(ps, ord_ps_NMDS_bray) + 
  geom_point(size = 2, 
	aes(color = factor(sample_data(ps)$Grid))) +
  stat_ellipse(type = "norm", aes(group = factor(sample_data(ps)$Grid))) +
  labs(color = "Grid") + 
  ggtitle("NDMS Ordination Plot") + 
  theme_classic()
dev.off()
# RDA bray
pdf("./plots/RDAbray.pdf")
plot_ordination(ps, ord_ps_RDA_bray) + 
  geom_point(size = 2, 
	aes(color = factor(sample_data(ps)$Grid))) +
  stat_ellipse(type = "norm", aes(group = factor(sample_data(ps)$Grid)))+
  labs(color = "Grid") + 
  ggtitle("RDA Ordination Plot") + theme_classic()
dev.off()
# DCA bray
pdf("./plots/DCAbray.pdf")
plot_ordination(ps, ord_ps_DCA_bray) + 
  geom_point(size = 2, 
	aes(color = factor(sample_data(ps)$Grid))) +
  stat_ellipse(type = "norm", aes(group = factor(sample_data(ps)$Grid))) +
  labs(color = "Grid") + 
  ggtitle("DCA Ordination Plot") + theme_classic()
dev.off()

# subset for only resampled squirrels
resampled.ps <- prune_samples(sample_data(ps)$Squirrel.ID %in% 
	c('10050', '10061', '10082', '10086', '10112', '10119', '10121', '10125', '10126', '10130', '10131', '10133', '10137', '10138', '10140', '10145', '10153', '10161', '10165', '10168', '10176', '10177', '10181', '10182', '10183', '10184', '10200', '10212', '10213', '10228', '10239', '10260', '10262', '10264', '10265', '10273', '10274', '10275', '10278', '10285', '10304', '10305', '10307', '10318', '10320', '10322', '10324', '10334', '10338', '10340', '10341', '10342', '10355', '10369', '10370', '10374', '10375', '10376', '10377', '10389', '10390', '10391', '10393', '10395', '10408', '10410', '10411', '10416', '10461', '10517', '10569', '10571', '10588', '10627', '10636', '10640', '10642', '10645', '10646', '10658', '10678', '10679', '10686', '10699', '10701', '10713', '10724', '10725', '10736', '10746', '10747', '10748', '10752', '10754', '10759', '10760', '10761', '10782', '10784', '10818', '10828', '10939', '10960', '11150', '11164', '11171', '11176', '11177', '11181', '11201', '11202', '11203', '11210', '11211', '11254', '11256', '11266', '11272', '11280', '11282', '11319', '11324', '11344', '11488', '11489', '11882', '13659', '13660', '13663', '13680', '13682', '13688', '13692', '8036', '8082', '8095', '8297', '8349', '8385', '8386', '8404', '8407', '8497', '8505', '8510', '8513', '8529', '8557', '8565', '8572', '8576', '8579', '8582', '8588', '8592', '8619', '8636', '8643', '8656', '8657', '8659', '8662', '8680', '8690', '8691', '8700', '8703', '8720', '8726', '8727', '8728', '8729', '8730', '8732', '8735', '8745', '8755', '8771', '8815'), 
	ps)
# make sure the sampling worked - which it has
nsamples(resampled.ps) # 725


shannon <- read_qza("../core-metrics-phylogenetic/shannon_vector.qza")
head(shannon$data)

evenness <- read_qza("../core-metrics-phylogenetic/evenness_vector.qza")

# save environment to file
save.image(file = "RedSquirrel.RData")

# to change RAM allocation for R
# https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos
#cd ~
#touch .Renviron
#open .Renviron
# in .Renviron, type in this as the first line (no comment symbol)
#R_MAX_VSIZE=100Gb 