#---
#title: "Test Run for step.env with ID"
#author: "Alicia Halhed"
#date: "02/15/2020"
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

# start analysis
print("Starting initial data preparation")
print("Access and plot XY data")
XY_sub <- XY_year(rs_q2_metadata, "KL", 2008)
print("Computing Euclidean Distances")
eDIST <- dist(XY_sub)

# get OTUs (aka a community matrix)
print("Build community object (OTU table) for grid/year")
comm_obj <- otu_table(ps) %>% as.matrix %>% 
    as.data.frame %>% 
    t %>% as.data.frame %>% 
    subset(., rownames(.) %in% rownames(XY_sub)) %>%
    .[ rowSums(.)>0, ]
# sample ID's are rownames
# get the metadata subset
print("Extract metadata for grid/year")
meta_sub <- rs_q2_metadata %>% 
  subset(., SampleID %in% rownames(XY_sub)) %>%
  select_if(~ !any(is.na(.))) %>% 
  # age and birth year are collinear
  # should I add the squirrel id here? dam/sire id has misssingness
  select(Sex, Age, Month, Season, CollectionDate, BirthYear, Squirrel.ID)
# Remove objects we're done with
print("Removing phyloseq obejct and full metadata data frame")
rm(rs_q2_metadata, ps)

# unweighted PCNM
print("Unweighted PCNM")
UWpcnm <- pcnm(eDIST)
UWpcnm$vectors

# test with RDA
print("Testing with RDA (full model)")
abFrac <- rda(decostand(comm_obj, "hel") ~ ., meta_sub)
abFrac # Full model
anova(abFrac, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
RsquareAdj(abFrac)

# forward selection for parsimonious model
print("Forward selection for parsimonious model")
# env variables
print("Environmental variables")
abFrac0 <- rda(decostand(comm_obj, "hel") ~ 1, meta_sub) # Reduced model
step.env <- ordiR2step(abFrac0, scope = formula(abFrac))
step.env # an rda model, with the final model predictor variables
# this is a summary of the selection process
step.env$anova
anova(step.env)

