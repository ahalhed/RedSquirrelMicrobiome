#---
#title: "PCNM for Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "04/13/2020"
#---

print("Set up (working directory, theme, and packages)")
# set working directory
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial")

# attach required packages
library(qiime2R)
library(phyloseq)
library(vegan)
library(zCompositions)
# devtools::install_github('ggloor/CoDaSeq/CoDaSeq')
library(CoDaSeq)
library(tidyverse)

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
ps <- qza_to_phyloseq(features = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/filtered-table-10.qza",
                      tree = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/trees/rooted_tree.qza",
                      taxonomy = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/taxonomy/SILVA-taxonomy-10.qza",
                      metadata = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/input/RS_meta.tsv")

# based on the meta function from the microbiome package
# I don't want to load a whole package for one function
print("Read in the metadata")
rs_q2_metadata <- as(sample_data(ps), "data.frame")
rownames(rs_q2_metadata) <- sample_names(ps)

# example in https://github.com/ggloor/CoDaSeq/blob/master/Intro_tiger_ladybug.Rmd
print("Aitchison transformation")
# rows are OTUs
# impute the OTU table
OTUimp <- otu_table(ps) %>% as.data.frame %>% # all OTUs
  cmultRepl(., label=0, method="CZM")
# compute the aitchison values
OTUclr <- codaSeq.clr(OTUimp)
mean.clr <- apply(OTUclr, 2, mean)
var.clr <- apply(OTUclr, 2, var)
# unload packages we'r done with
detach("package:CoDaSeq", unload = TRUE)
detach("package:qiime2R", unload = TRUE)
detach("package:phyloseq", unload = TRUE)
detach("package:zCompositions", unload = TRUE)

# start analysis
print("Starting initial data preparation")
print("Access and plot XY data")
XY_sub <- XY_year(rs_q2_metadata, "AG", 2008)
# plotting the locations
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/plots/AG2008_XY.pdf")
XY_sub %>% ggplot(aes(x = `Location.X`, y = `Location.Y`)) + 
  geom_point() + 
  coord_fixed()
dev.off()
print("Computing Euclidean Distances")
eDIST <- dist(XY_sub)
print("Maximum Euclidean Distance")
max_dist(eDIST)

# get OTUs (aka a community matrix)
print("Build community object (OTU table) for grid/year")
comm_obj <- OTUclr %>% t %>% as.data.frame %>%
    subset(., rownames(.) %in% rownames(XY_sub)) %>%
    .[ , colSums(.)>0 ]
# sample ID's are rownames
# get the metadata subset
print("Extract metadata for grid/year")
meta_sub <- rs_q2_metadata %>% 
  subset(., rownames(.) %in% rownames(XY_sub)) %>%
  select_if(~ !any(is.na(.))) %>% 
  # age and birth year are collinear
  # remove location information
  subset(., select = -c(Location.X, Location.Y, Location, Date)) %>% 
  # select columns with more than one level
  .[sapply(., function(x) length(unique(x))>1)]
# Remove objects we're done with
print("Removing phyloseq obejct and full metadata data frame")
rm(rs_q2_metadata, ps)

# unweighted PCNM
print("Unweighted PCNM")
UWpcnm <- pcnm(eDIST)
UWpcnm$vectors
# plot with ordisurf
print("Plotting first three PCNM axes with ordisurf")
# replace grid-year with values used in this script
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/plots/AG2008_ordisurf123.pdf")
par(mfrow=c(1,3))
ordisurf(XY_sub, scores(UWpcnm, choi=1), bubble = 4, main = "PCNM 1")
ordisurf(XY_sub, scores(UWpcnm, choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_sub, scores(UWpcnm, choi=3), bubble = 4, main = "PCNM 3")
dev.off()

# weighted PCNM
print("Weighted PCNM")
Wpcnm <- pcnm(eDIST, w = rowSums(comm_obj)/sum(comm_obj))
Wpcnm$vectors

# computing CCA with weighted PCNM
print("CCA with weighted PCNM")
# not including Birth year here due to collinearity with age (same information)
cca_sub <- cca(comm_obj ~ scores(Wpcnm) + Sex + Season + Age, meta_sub)
summary(cca_sub)


# removing the MSO stuff, since it didn't like the aitchison stuff

# test with RDA
print("Testing with RDA (full model)")
abFrac <- rda(comm_obj ~ ., meta_sub)
abFrac # Full model
anova(abFrac, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
RsquareAdj(abFrac)

# Test fraction [a] using partial RDA:
print("Testing with partial RDA (fraction [a])")
aFrac <- rda(comm_obj ~ . + Condition(scores(UWpcnm)), data = meta_sub)
# the anova was taking forever to run
anova(aFrac, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
RsquareAdj(aFrac)

# forward selection for parsimonious model
print("Forward selection for parsimonious model")
# env variables
print("Environmental variables")
abFrac0 <- rda(comm_obj ~ 1, meta_sub) # Reduced model
step.env <- ordiR2step(abFrac0, scope = formula(abFrac))
step.env # an rda model, with the final model predictor variables
# this is a summary of the selection process
print("Summary of environmental selection process")
step.env$anova
print("ANOVA on full environmental selection")
anova(step.env)
# save plot
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/plots/AG2008_step_env.pdf")
plot(step.env)
dev.off()

# spatial variables
print("Spatial variables")
pcnm_df <- as.data.frame(scores(UWpcnm))
bcFrac <- rda(comm_obj ~ ., pcnm_df) # Full model
bcFrac0 <- rda(comm_obj ~ 1, pcnm_df) # Reduced model
step.space <- ordiR2step(bcFrac0, scope = formula(bcFrac))
step.space
# this is a summary of the selection process
print("Summary of spatial selection process")
step.space$anova
print("ANOVA on full spatial selection")
anova(step.space)
# save plot
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/plots/AG2008_step_space.pdf")
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
               data = meta_sub[, names(step.env$terminfo$ordered)])
mod.pars
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/plots/AG2008_mod_pars.pdf")
plot(mod.pars)
dev.off()

