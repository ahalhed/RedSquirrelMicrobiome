#---
#title: "PCNM for Core Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "04/24/2020"
#---

print("Set up (working directory, theme, and packages)")
# set working directory
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial")

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
  met <- rownames_to_column(metadata, var = "SampleID")
  df1 <- subset(met, Grid == grid, 
         select = c("SampleID", "Location.X", "Location.Y", "Year"))
  df2 <- subset(df1, Year == year, 
         select = c("SampleID", "Location.X", "Location.Y"))
  df3 <- column_to_rownames(remove_rownames(df2), var = "SampleID")
  return(df3)
}

# maximum distance
max_dist <- function(dm) {
  df1 <- as.data.frame(as.matrix(dm))
  # functions is soft deprecated (replace with functions or lambdas)
  summ <- summarise_each(df1, ~ max(df1, na.rm=TRUE))
  m <- apply(summ, 1, max)
  return(m)
}

# get the data
print("Read in the Data")
print("Building phyloseq object")
ps <- qza_to_phyloseq(features = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/filtered-table.qza",
                      tree = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/trees/rooted_tree.qza",
                      taxonomy = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/taxonomy/GG-taxonomy.qza",
                      metadata = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/input/RS_meta.tsv")

# based on the meta function from the microbiome package
# I don't want to load a whole package for one function
print("Read in the metadata")
rs_q2_metadata <- as(sample_data(ps), "data.frame")
rownames(rs_q2_metadata) <- sample_names(ps)

# start analysis
print("Starting initial data preparation")
print("Access and plot XY data")
XY_sub <- XY_year(rs_q2_metadata, "KL", 2009)
print("Computing Euclidean Distances")
eDIST <- dist(XY_sub)
print("Maximum Euclidean Distance")
max_dist(eDIST)

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
# subset the samples from the core microbiome
print("Build the core community object (OTU table) for grid/year")
comm_obj <- OTU_core %>% 
  subset(., rownames(.) %in% rownames(XY_sub)) %>%
  .[ , colSums(.)>0 ]


# sample ID's are rownames
# get the metadata subset
print("Extract metadata for grid/year")
meta_sub <- rs_q2_metadata %>% 
  subset(., SampleID %in% rownames(XY_sub)) %>%
  select_if(~ !any(is.na(.))) %>% 
  # age and birth year are collinear
  # should I add the squirrel id here? dam/sire id has misssingness
  select(Sex, Age, Month, Season, CollectionDate, BirthYear)

# Remove objects we're done with
print("Removing pobjects that are no longer needed")
rm(rs_q2_metadata, cOTU, OTU95, RA_001, RA, RA_full, OTU_full, ps, OTU_core)

# unweighted PCNM
print("Unweighted PCNM")
UWpcnm <- pcnm(eDIST)
UWpcnm$vectors

# weighted PCNM
print("Weighted PCNM")
Wpcnm <- pcnm(eDIST, w = rowSums(comm_obj)/sum(comm_obj))
Wpcnm$vectors

# computing CCA with weighted PCNM
print("CCA with weighted PCNM")
# not including Birth year here due to collinearity with age (same information)
cca_sub <- cca(comm_obj ~ scores(Wpcnm) + Sex + Season + Age, meta_sub)
summary(cca_sub)

# Spatial partitioning of CCA (mso)
print("Multiscale ordination")
mso_sub <- mso(cca_sub, XY_sub)
# plot
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/plots/core_KL2009_mso.pdf")
msoplot(mso_sub, ylim = c(0, 45), main="2009 KL")
dev.off()

print("CCA retaining spatial patterns for MSO")
print("Constrained")
cca_con <- cca(log(comm_obj + 1) ~ Sex + Season + Age, meta_sub)
mso_sub2 <- mso(cca_con, XY_sub)
mso_sub2
print("Unconstrained")
cca_un <- cca(log(comm_obj + 1))
mso_sub3 <- mso(cca_un, XY_sub)
mso_sub3
# plot
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/plots/core_KL2009_mso2.pdf")
par(mfrow=c(1,2))
msoplot(mso_sub2, ylim = c(0,0.25), main="Constrained Ordination")
msoplot(mso_sub3, ylim = c(0,0.25), main="Unconstrained Ordination")
dev.off()

# Variance partitioning
print("Variance partitioning")
vp_mod1 <- varpart(comm_obj,  ~ ., scores(UWpcnm), data=meta_sub, transfo = "hel")
vp_mod1
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/plots/core_KL2009_vp_mod1.pdf")
plot(vp_mod1)
dev.off()

# test with RDA
print("Testing with RDA (full model)")
abFrac <- rda(decostand(comm_obj, "hel") ~ ., meta_sub)
abFrac # Full model
anova(abFrac, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
RsquareAdj(abFrac)

# Test fraction [a] using partial RDA:
print("Testing with partial RDA (fraction [a])")
aFrac <- rda(decostand(comm_obj, "hel") ~ . + Condition(scores(UWpcnm)), data = meta_sub)
# the anova was taking forever to run
anova(aFrac, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
RsquareAdj(aFrac)

# forward selection for parsimonious model
print("Forward selection for parsimonious model")
# env variables
print("Environmental variables")
abFrac0 <- rda(decostand(comm_obj, "hel") ~ 1, meta_sub) # Reduced model
step.env <- ordiR2step(abFrac0, scope = formula(abFrac))
step.env # an rda model, with the final model predictor variables
print("Summary of environmental selection process")
step.env$anova
print("ANOVA on full environmental selection")
anova(step.env)
# save plot
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/plots/core_KL2009_step_env.pdf")
plot(step.env)
dev.off()

# spatial variables
print("Spatial variables")
pcnm_df <- as.data.frame(scores(UWpcnm))
bcFrac <- rda(decostand(comm_obj, "hel") ~ ., pcnm_df) # Full model
bcFrac0 <- rda(decostand(comm_obj, "hel") ~ 1, pcnm_df) # Reduced model
step.space <- ordiR2step(bcFrac0, scope = formula(bcFrac))
step.space
print("Summary of spatial selection process")
step.space$anova
print("ANOVA on full spatial selection")
anova(step.space)
# save plot
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/plots/core_KL2009_step_space.pdf")
plot(step.space)
dev.off()

# Partition Bray-Curtis dissimilarities
print("Partition Bray-Curtis dissimilarities")
varpart(vegdist(comm_obj), ~ ., scores(UWpcnm), data = meta_sub)

# variation decomposition with parsimonious variables
print("Variation decomposition with parsimonious variables")
mod.pars <- varpart(comm_obj, ~ ., 
               pcnm_df[, names(step.space$terminfo$ordered)], 
               data = meta_sub[, names(step.env$terminfo$ordered)],
               transfo = "hel")
mod.pars
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/plots/core_KL2009_mod_pars.pdf")
plot(mod.pars)
dev.off()
