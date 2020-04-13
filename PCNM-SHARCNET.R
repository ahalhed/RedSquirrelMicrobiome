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

print("Read in the Data")
print("Building phyloseq object")
ps <- qza_to_phyloseq(features = "/home/ahalhed/red-squirrel-w2020/filtered-table.qza",
                      tree = "/home/ahalhed/red-squirrel-w2020/trees/rooted_tree.qza",
                      taxonomy = "/home/ahalhed/red-squirrel-w2020/taxonomy/GG-taxonomy.qza",
                      metadata = "/home/ahalhed/red-squirrel-w2020/input/RS_meta.tsv")

print("Read in the metadata")
rs_q2_metadata <- read.table("~/OneDrive - University of Guelph/Alicia's Thesis/red-squirrel-w2020/RS_meta.tsv", sep="\t")
colnames(rs_q2_metadata) <- c("SampleID", "Grid", "Location X", "Location Y", "Sex", "Age", "Month", "Season", "Year", "Squirrel.ID", "SireID", "DamID", "CollectionDate", "FoodSupplement", "BirthYear", "Location", "Date")


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

# numeric metadata (depends on dplyr)
Grid_numeric <- function(metadata, grid, year) {
  df1 <- subset(metadata, Grid == grid) %>%
    remove_rownames() %>%
    column_to_rownames(., var = "SampleID") %>%
    subset(., Year == year)
  df2 <- df1 %>%
    # want sex and food supplement to be numeric
    # male 0, female 1; no 0, yes 1
    mutate(Sex = as.integer(recode(Sex, "Male" = 0, "Female" = 1)),
           FoodSupplement = as.integer(recode(FoodSupplement, 
                                   "No" = 0, "Yes" = 1)))
  df3 <- Filter(is.numeric, df2)
  drop <- c("Location X", "Location Y", "Date")
  df4 <- df3[,!(names(df3) %in% drop)]
  return(df4)
}


# XY data
XY_10_KL <- XY_year(rs_q2_metadata, "KL", 2010)
# numeric data for CCA
numeric_10_KL <- Grid_numeric(rs_q2_metadata, "KL", 2010)


# Compute Euclidean Distance
e_10_KL <- dist(XY_10_KL)
# checking what the highest values are in the dm
max_dist(e_10_KL)



# unweighted PCNM
pcnm_10_KL <- pcnm(e_10_KL)
# weighted PCNM
# get OTUs (aka a community matrix)
comm_10_KL <- otu_table(ps) %>% as.matrix %>% 
    as.data.frame %>% 
    t %>% as.data.frame %>% 
    subset(., rownames(.) %in% rownames(XY_10_KL)) %>%
    .[ rowSums(.)>0, ]
# sample ID's are rownames
comm_10_KL
# get the metadata subset
meta_10_KL <- rs_q2_metadata %>% 
  subset(., SampleID %in% rownames(XY_10_KL)) %>%
  select_if(~ !any(is.na(.))) %>% 
  # age and birth year are collinear
  select(Sex, Age, Month, Season, CollectionDate, BirthYear)
# calculate weighted PCNM
Wpcnm_10_KL <- pcnm(e_10_KL, 
                    w = rowSums(comm_10_KL)/sum(comm_10_KL))
# CCA - might need to change explanatories
# not including Birth year here due to collinearity with age (same information)
cca_10_KL <- cca(comm_10_KL ~ scores(Wpcnm_10_KL) + Sex + Season + Age, meta_10_KL)


# look at the PCNM values
pcnm_10_KL$vectors
Wpcnm_10_KL$vectors
# look at CCA summary
# what is significance? most important?
summary(cca_10_KL)


# Spatial partitioning of CCA (mso)


# multiscale ordination
mso_10_KL <- mso(cca_10_KL, XY_10_KL)
# plot
msoplot(mso_10_KL, ylim = c(0, 45), main="2010 KL")


# plotting the locations
XY_10_KL %>% ggplot(aes(x = `Location X`, y = `Location Y`)) + 
  geom_point() + 
  coord_fixed()



par(mfrow=c(1,3))
ordisurf(XY_10_KL, scores(pcnm_10_KL, choi=1), bubble = 4, main = "PCNM 1")
ordisurf(XY_10_KL, scores(pcnm_10_KL, choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_10_KL, scores(pcnm_10_KL, choi=3), bubble = 4, main = "PCNM 3")


# Variance partitioning
vp_mod1 <- varpart(comm_10_KL,  ~ ., scores(pcnm_10_KL), data=meta_10_KL, transfo = "hel")
vp_mod1
plot(vp_mod1)


# test with RDA
abFrac <- rda(decostand(comm_10_KL, "hel") ~ ., meta_10_KL)
abFrac # Full model
# RsquareAdj gives the same result as component [a] of varpart
RsquareAdj(abFrac)



# Test fraction [a] using partial RDA:
aFrac <- rda(decostand(comm_10_KL, "hel") ~ . + Condition(scores(pcnm_10_KL)), data = meta_10_KL)
# the anova was taking forever to run
anova(aFrac, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
RsquareAdj(aFrac)


# forward selection for parsimonious model
# env variables
abFrac0 <- rda(decostand(comm_10_KL, "hel") ~ 1, meta_10_KL) # Reduced model
# Here is where the magic happens, but almost automatically!
step.env <- ordiR2step(abFrac0, scope = formula(abFrac))
str(step.env) # this seems to be a very complicated
step.env # but it is actually just an rda model, with the final model predictor variables
anova(step.env) # so you can do the same stuff as before
plot(step.env)
step.env$anova # and if you are interested, this is a summary of the selection process

# spatial variables
rs_pcnm <- as.data.frame(scores(pcnm_10_KL))
bcFrac <- rda(decostand(comm_10_KL, "hel") ~ ., rs_pcnm) # Full model
bcFrac0 <- rda(decostand(comm_10_KL, "hel") ~ 1, rs_pcnm) # Reduced model
step.space <- ordiR2step(bcFrac0, scope = formula(bcFrac))
step.space$anova
plot(step.space)

# variation decomposition with parsimonious variables
mod.pars <- varpart(comm_10_KL, ~ ., 
               rs_pcnm[, names(step.space$terminfo$ordered)], 
               data = meta_10_KL[, names(step.env$terminfo$ordered)],
               transfo = "hel")
mod.pars
plot(mod.pars)

# Partition Bray-Curtis dissimilarities
varpart(vegdist(comm_10_KL), ~ ., scores(pcnm_10_KL), data = meta_10_KL)
