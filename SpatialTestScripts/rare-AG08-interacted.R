#---
#title: "PCNM for Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "02/13/2020"
#---

print("Set up (working directory, theme, and packages)")
# set working directory
setwd("/home/ahalhed/red-squirrel-w2020/R-env")

# attach required packages
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(vegan)

# set theme for plots
theme_set(theme_bw())

print("Initiate functions for analysis")
# subset the XY's by grid and year
XY_year <- function(metadata, grid, year) {
  df1 <- subset(metadata, Grid == grid, 
                select = c("SampleID", "Location X", "Location Y", "Year"))
  df2 <- subset(df1, Year == year, 
                select = c("SampleID", "Location X", "Location Y"))
  df3 <- column_to_rownames(remove_rownames(df2), var = "SampleID")
  return(df3)
}

# get the data
print("Read in the Data")
print("Building phyloseq object")
ps <- qza_to_phyloseq(features = "/home/ahalhed/red-squirrel-w2020/filtered-table.qza",
                      tree = "/home/ahalhed/red-squirrel-w2020/trees/rooted_tree.qza",
                      taxonomy = "/home/ahalhed/red-squirrel-w2020/taxonomy/GG-taxonomy.qza",
                      metadata = "/home/ahalhed/red-squirrel-w2020/input/RS_meta.tsv")

print("Read in the metadata")
rs_q2_metadata <- read.table("/home/ahalhed/red-squirrel-w2020/input/RS_meta.tsv", sep="\t")
colnames(rs_q2_metadata) <- c("SampleID", "Grid", "Location X", "Location Y", "Sex", "Age", "Month", "Season", "Year", "Squirrel.ID", "SireID", "DamID", "CollectionDate", "FoodSupplement", "BirthYear", "Location", "Date")

print("Getting the XY data")
XY_sub <- XY_year(rs_q2_metadata, "AG", 2008)
print("Computing Euclidean Distances")
eDIST <- dist(XY_sub)
# get OTUs (aka a community matrix)
print("Finding core microbiome")
print("Accessing full OTU table as a community object")
OTU_full <- otu_table(ps) %>% as.matrix %>% 
  as.data.frame %>% 
  t %>% as.data.frame 
# transform the counts to relative abundance
print("Transforming the OTU counts to relative abundance")
RA_full <- transform_sample_counts(ps, function(x) x/sum(x)) %>% 
  otu_table %>% as.matrix %>% as.data.frame %>% 
  t %>% as.data.frame
# 95% threshold
# minimum number of samples that an OTU must occur in 
# to be considered as part of the "core" microbiome
print("Calculating 95% threshold")
th <- 0.95*nrow(RA_full)
th
# Access the number of non-zero occurrences of each OTU
print("Accessing the number of non-zero occurrences of each OTU")
nonzero <- function(x) sum(x!=0)
nz <- plyr::numcolwise(nonzero)(RA_full)
# average relative abundance across all 909 samples
print("Calculating average relative abundance")
RA <- RA_full %>% colSums() %>% as.data.frame
RA_001 <- RA %>%
  rownames_to_column(var = "OTU") %>% 
  rename("ColSum" = ".") %>% 
  mutate(RelativeAbundance = .$ColSum/nrow(OTU_full)) %>%
  # greater than 0.001 (or 0.1%)
  filter(RelativeAbundance > 0.001)
# extract the OTUs that fall above the 95% threshold
# based on the number of non-zero occurrences of each OTU
print("Extract the OTUs that occur in >95% of samples")
OTU95 <- t(nz) %>% as.data.frame %>% 
  rownames_to_column(var = "OTU") %>% 
  rename("Occurrences" = "V1") %>%
  subset(., Occurrences > th)
# combine to ensure that those in 95% of the samples
# have a total relative abundance of 0.1%
print("Find OTUs that occur in >95% of samples and have an average relative abundance >0.1%")
cOTU <- inner_join(OTU95, RA_001)
print("How many OTUs meet these criteria?")
nrow(cOTU)
print("Subset the OTU table to find rare OTUs")
OTU_rare <- OTU_full %>% 
  select(-one_of(cOTU$OTU))
# subset the samples from the rare microbiome
print("Build the rare community object (OTU table) for grid/year")
comm_obj <- OTU_rare %>% 
  subset(., rownames(.) %in% rownames(XY_sub)) %>%
  .[ rowSums(.)>0, ]
# sample ID's are rownames
# get the metadata subset
print("Extract metadata for grid/year")
meta_sub <- rs_q2_metadata %>% 
  subset(., SampleID %in% rownames(XY_sub)) %>%
  select_if(~ !any(is.na(.))) %>% 
  # age and birth year are collinear
  select(Sex, Age, Month, Season, CollectionDate, BirthYear)
# Remove objects we're done with
print("Removing phyloseq obejct and full metadata data frame")
rm(rs_q2_metadata, ps)

# unweighted PCNM
print("Unweighted PCNM")
UWpcnm <- pcnm(eDIST)
# spatial variables
print("Spatial variables")
pcnm_df <- as.data.frame(scores(UWpcnm))

# picking up where the batch job got killed
step.space <- rda(formula = decostand(comm_obj, "hel") ~
                    PCNM55 + PCNM56 + PCNM22 + PCNM27 + PCNM74 + PCNM31 + PCNM73, 
                  data = pcnm_df)
step.space
# Call: rda(formula = decostand(comm_obj, "hel") ~ PCNM55 + PCNM56 +
# PCNM22 + PCNM27 + PCNM74 + PCNM31 + PCNM73, data = pcnm_df)

# Inertia Proportion Rank
# Total         0.71540    1.00000     
# Constrained   0.03157    0.04413    7
# Unconstrained 0.68384    0.95587  210
# Inertia is variance 

# Eigenvalues for constrained axes:
#   RDA1     RDA2     RDA3     RDA4     RDA5     RDA6     RDA7 
# 0.010326 0.005376 0.003844 0.003605 0.003060 0.002924 0.002433 

# Eigenvalues for unconstrained axes:
#   PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8 
# 0.06126 0.02458 0.01939 0.01663 0.01345 0.01039 0.00968 0.00920 
# (Showing 8 of 210 unconstrained eigenvalues)



anova(step.space)
# this is a summary of the selection process
step.space$anova
# save plot
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/core_AG2008_step_space.pdf")
plot(step.space)
dev.off()

# Partition Bray-Curtis dissimilarities
print("Partition Bray-Curtis dissimilarities")
varpart(vegdist(comm_obj), ~ ., scores(UWpcnm), data = meta_sub)

# variation decomposition with parsimonious variables
# if environmental or spatial variables don't have significant
# axes, the below will fail
print("Variation decomposition with parsimonious variables")
mod.pars <- varpart(comm_obj, ~ ., 
                    pcnm_df[, names(step.space$terminfo$ordered)], 
                    data = meta_sub[, names(step.env$terminfo$ordered)],
                    transfo = "hel")
mod.pars
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/core_AG2008_mod_pars.pdf")
plot(mod.pars)
dev.off()

