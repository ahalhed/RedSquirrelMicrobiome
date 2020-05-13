#---
#title: "PCNM for Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "04/13/2020"
#---

print("Set up (working directory, theme, and packages)")
# set working directory
setwd("/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial")

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
XY_sub <- XY_year(rs_q2_metadata, "SU", 2008)
# plotting the locations
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/SU2008_XY.pdf")
XY_sub %>% ggplot(aes(x = `Location X`, y = `Location Y`)) + 
  geom_point() + 
  coord_fixed()
dev.off()
print("Computing Euclidean Distances")
eDIST <- dist(XY_sub)
print("Maximum Euclidean Distance")
max_dist(eDIST)

# get OTUs (aka a community matrix)
print("Build community object (OTU table) for grid/year")
comm_obj <- otu_table(ps) %>% as.matrix %>% 
    as.data.frame %>% 
    t %>% as.data.frame %>% 
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
  select(Age, CollectionDate, BirthYear)
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
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/SU2008_ordisurf123.pdf")
par(mfrow=c(1,3))
# these can be adjusted afterwards, once we know which are significant
# see section starting at 161
ordisurf(XY_sub, scores(UWpcnm, choi=1), bubble = 4, main = "PCNM 1")
ordisurf(XY_sub, scores(UWpcnm, choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_sub, scores(UWpcnm, choi=3), bubble = 4, main = "PCNM 3")
dev.off()

# weighted PCNM (labelled as WUWpcnm in scripts run - a find/replace typo)
print("Weighted PCNM")
Wpcnm <- pcnm(eDIST, w = rowSums(comm_obj)/sum(comm_obj))
Wpcnm$vectors

# computing CCA with weighted PCNM
print("CCA with weighted PCNM")
# not including Birth year here due to collinearity with age (same information)
cca_sub <- cca(comm_obj ~ scores(Wpcnm) + Age, meta_sub)
summary(cca_sub)
# Spatial partitioning of CCA (mso)
print("Multiscale ordination")
mso_sub <- mso(cca_sub, XY_sub)
# plot
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/SU2008_mso.pdf")
msoplot(mso_sub, ylim = c(0, 45), main="2008 SU")
dev.off()

print("CCA retaining spatial patterns for MSO")
print("Constrained")
cca_con <- cca(log(comm_obj + 1) ~ Age, meta_sub)
mso_sub2 <- mso(cca_con, XY_sub)
mso_sub2
print("Unconstrained")
cca_un <- cca(log(comm_obj + 1))
mso_sub3 <- mso(cca_un, XY_sub)
mso_sub3
# plot
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/SU2008_mso2.pdf")
par(mfrow=c(1,2))
msoplot(mso_sub2, ylim = c(0,0.25), main="Constrained Ordination")
msoplot(mso_sub3, ylim = c(0,0.25), main="Unconstrained Ordination")
dev.off()

# notes on data
# Grids JO and SU only has females from month 5, late spring
# Grids CH and LL only has females

# Variance partitioning
print("Variance partitioning")
vp_mod1 <- varpart(comm_obj,  ~ ., scores(UWpcnm), data=meta_sub, transfo = "hel")
vp_mod1
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/SU2008_vp_mod1.pdf")
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
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/SU2008_step_env.pdf")
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
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/SU2008_step_space.pdf")
plot(step.space)
dev.off()

# Partition Bray-Curtis dissimilarities
print("Partition Bray-Curtis dissimilarities")
varpart(vegdist(comm_obj), ~ ., scores(UWpcnm), data = meta_sub)

# variation decomposition with parsimonious variables
# probably going to fail because of insignificant environmental variables
print("Variation decomposition with parsimonious variables")
mod.pars <- varpart(comm_obj, ~ ., 
               pcnm_df[, names(step.space$terminfo$ordered)], 
               data = meta_sub[, names(step.env$terminfo$ordered)],
               transfo = "hel")
mod.pars
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/SU2008_mod_pars.pdf")
plot(mod.pars)
dev.off()

