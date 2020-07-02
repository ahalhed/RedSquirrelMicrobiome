#---
#title: "Counting OTUs per grid/year"
#author: "Alicia Halhed"
#date: "05/07/2020"
#---

print("Set up (working directory, theme, and packages)")
# set working directory
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env")

# attach required packages
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(vegan)

print("Initiate function for analysis")
# subset the metadata by grid and year
met_year <- function(metadata, grid, year) {
  met <- rownames_to_column(metadata, var = "SampleID")
  df1 <- subset(met, Grid == grid)
  df2 <- subset(df1, Year == year)
  df3 <- column_to_rownames(remove_rownames(df2), var = "SampleID")
  return(df3)
}
# subset the XY's by grid, year, month
XY_month <- function(metadata, grid, year, month) {
  met <- rownames_to_column(metadata, var = "SampleID")
  df1 <- subset(met, Grid == grid, 
                select = c("SampleID", "Location.X", "Location.Y", "Year", "Month"))
  df2 <- subset(df1, Year == year, 
                select = c("SampleID", "Location.X", "Location.Y", "Month"))
  df3 <- subset(df2, Month == month, 
                select = c("SampleID", "Location.X", "Location.Y"))
  df4 <- column_to_rownames(remove_rownames(df3), var = "SampleID")
  return(df4)
}
# build community object
comm_obj <- function(XY, c) {
  # subset the OTUs (c is OTU table being subset)
  comm <- c %>%
    subset(., rownames(.) %in% rownames(XY)) %>%
    .[ , colSums(.)>0 ]
  return(comm)
}
# subset the XY's by grid and year
XY_year <- function(metadata, grid, year) {
  df1 <- subset(metadata, Grid == grid, 
                select = c("Location.X", "Location.Y", "Year"))
  df2 <- subset(df1, Year == year, 
                select = c("Location.X", "Location.Y"))
  return(df2)
}

# get the data
print("Read in the Data")
print("Building phyloseq object")
ps <- qza_to_phyloseq(features = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/filtered-table-10.qza",
                      tree = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/trees/rooted_tree.qza",
                      taxonomy = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/taxonomy/SILVA-taxonomy-10.qza",
                      metadata = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/input/RS_meta.tsv")
ps

# based on the meta function from the microbiome package
# I don't want to load a whole package for one function
print("Read in the metadata")
rs_q2_metadata <- as(sample_data(ps), "data.frame")
rownames(rs_q2_metadata) <- sample_names(ps)

# start analysis
print("Starting initial data preparation")
print("Access XY data")
# for all months
XY_sub <- XY_year(rs_q2_metadata, "AG", 2008)
# for individual months
met <- met_year(rs_q2_metadata, "AG", 2008)
# loop to create individual month data frames
for (month in unique(met$Month)) {
  df <- XY_month(rs_q2_metadata, "AG", 2008, month)
  assign(paste('Month',month,sep = ' '),df)
  rm(df, month)
}
# make a list of the data frames generated from the loop
XY_list <- do.call("list",
                   # searching the global environment for the pattern
                   mget(grep("Month", names(.GlobalEnv), value=TRUE)))

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
comm_core <- comm_obj(XY_sub, OTU_core)
print("Number of OTUs (columns) in core OTU table for AG 2008")
ncol(comm_core)
# subset the samples from the rare microbiome
print("Build the rare community object (OTU table) for grid/year")
comm_rare <- comm_obj(XY_sub, OTU_rare)
print("Number of OTUs (columns) in rare OTU table for AG 2008")
ncol(comm_rare)
# subset the samples from the full microbiome
print("Building community matrix for full OTU table")
comm_full <- comm_obj(XY_sub, OTU_full)
print("Number of OTUs (columns) in full OTU table for AG 2008")
ncol(comm_full)

print("Build the community object (OTU table) for grid/year/month")
mFull <- lapply(XY_list, comm_obj, c=OTU_full)
mCore <- lapply(XY_list, comm_obj, c=OTU_core)
mRare <- lapply(XY_list, comm_obj, c=OTU_rare)

print("Number of OTUs (columns) in full OTU table for months")
lapply(mFull, ncol)
print("Number of OTUs (columns) in core OTU table for months")
lapply(mCore, ncol)
print("Number of OTUs (columns) in rare OTU table for months")
lapply(mRare, ncol)

