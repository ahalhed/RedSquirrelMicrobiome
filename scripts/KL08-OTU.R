#---
#title: "Counting OTUs per grid/year"
#author: "Alicia Halhed"
#date: "05/07/2020"
#---

print("Set up (working directory, theme, and packages)")
# set working directory
setwd("/home/ahalhed/red-squirrel-w2020/R-env")

# attach required packages
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(vegan)

print("Initiate function for analysis")
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
ps

print("Read in the metadata")
rs_q2_metadata <- read.table("/home/ahalhed/red-squirrel-w2020/input/RS_meta.tsv", sep="\t")
colnames(rs_q2_metadata) <- c("SampleID", "Grid", "Location X", "Location Y", "Sex", "Age", "Month", "Season", "Year", "Squirrel.ID", "SireID", "DamID", "CollectionDate", "FoodSupplement", "BirthYear", "Location", "Date")

# start analysis
print("Starting initial data preparation")
print("Access and plot XY data")
XY_sub <- XY_year(rs_q2_metadata, "KL", 2008)
# plotting the locations
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
print("Subset the OTU table to find core OTUs")
OTU_core <- OTU_full[, cOTU$OTU]
print("Subset the OTU table to find rare OTUs")
OTU_rare <- OTU_full %>% 
  select(-one_of(cOTU$OTU))
# subset the samples from the core microbiome
print("Build the core community object (OTU table) for grid/year")
comm_core <- OTU_core %>% 
  subset(., rownames(.) %in% rownames(XY_sub)) %>%
  .[ colSums(.)>0, ]
print("Number of OTUs (columns) in core OTU table for KL 2008")
ncol(comm_core)
# subset the samples from the rare microbiome
print("Build the rare community object (OTU table) for grid/year")
comm_rare <- OTU_rare %>% 
  subset(., rownames(.) %in% rownames(XY_sub)) %>%
  .[ , colSums(.)>0 ]
print("Number of OTUs (columns) in rare OTU table for KL 2008")
ncol(comm_rare)

print("Building community matrix for full OTU table")
comm_full <- OTU_full %>% 
  subset(., rownames(.) %in% rownames(XY_sub)) %>%
  .[ colSums(.)>0, ]

print("Number of OTUs (columns) in full OTU table for KL 2008")
ncol(comm_full)
