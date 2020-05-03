#---
#title: "Significant Ordisurf Plots - Rare Microbiome"
#author: "Alicia Halhed"
#date: "04/22/2020"
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
print("Read in the metadata")
rs_q2_metadata <- read.table("/home/ahalhed/red-squirrel-w2020/input/RS_meta.tsv", sep="\t")
colnames(rs_q2_metadata) <- c("SampleID", "Grid", "Location X", "Location Y", "Sex", "Age", "Month", "Season", "Year", "Squirrel.ID", "SireID", "DamID", "CollectionDate", "FoodSupplement", "BirthYear", "Location", "Date")

# start analysis
print("Starting initial data preparation")
print("Access and plot XY data")
XY_AG <- XY_year(rs_q2_metadata, "AG", 2008)
XY_CH <- XY_year(rs_q2_metadata, "CH", 2008)
XY_JO <- XY_year(rs_q2_metadata, "JO", 2008)
XY_KL8 <- XY_year(rs_q2_metadata, "KL", 2008)
XY_KL10 <- XY_year(rs_q2_metadata, "KL", 2010)
XY_LL <- XY_year(rs_q2_metadata, "LL", 2008)
XY_SU <- XY_year(rs_q2_metadata, "SU", 2008)

print("Computing Euclidean Distances")
d_AG <- dist(XY_AG)
d_CH <- dist(XY_CH)
d_JO <- dist(XY_JO)
d_KL8 <- dist(XY_KL8)
d_KL10 <- dist(XY_KL10)
d_LL <- dist(XY_LL)
d_SU <- dist(XY_SU)


# Remove objects we're done with
print("Removing full metadata data frame")
rm(rs_q2_metadata)

# unweighted PCNM
print("Unweighted PCNM")
AG <- pcnm(d_AG)
CH <- pcnm(d_CH)
JO <- pcnm(d_JO)
KL8 <- pcnm(d_KL8)
KL10 <- pcnm(d_KL10)
LL <- pcnm(d_LL)
SU <- pcnm(d_SU)

# plot with ordisurf
print("Plotting significant PCNM axes with ordisurf")

print("AG 2008")
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/rare_AG2008_ordiSIG.pdf")
par(mfrow=c(3,2))
ordisurf(XY_AG, scores(AG, choi=55), bubble = 4, main = "PCNM 55")
ordisurf(XY_AG, scores(AG, choi=56), bubble = 4, main = "PCNM 56")
ordisurf(XY_AG, scores(AG, choi=22), bubble = 4, main = "PCNM 22")
ordisurf(XY_AG, scores(AG, choi=27), bubble = 4, main = "PCNM 27")
ordisurf(XY_AG, scores(AG, choi=31), bubble = 4, main = "PCNM 31")
ordisurf(XY_AG, scores(AG, choi=73), bubble = 4, main = "PCNM 73")
dev.off()

print("CH 2008")
print("No significant core axes")

print("JO 2008")
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/rare_JO2008_ordiSIG.pdf")
par(mfrow=c(1,1))
ordisurf(XY_JO, scores(JO, choi=1), bubble = 4, main = "PCNM 1")
dev.off()

print("KL 2008")
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/rare_KL2008_ordiSIG.pdf")
par(mfrow=c(1,1))
ordisurf(XY_KL8, scores(KL8, choi=2), bubble = 4, main = "PCNM 2")
dev.off()

print("KL 2009")
print("No significant rare axes")

print("KL 2010")
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/rare_KL2010_ordiSIG.pdf")
par(mfrow=c(2,2))
ordisurf(XY_KL10, scores(KL10, choi=33), bubble = 4, main = "PCNM 33")
ordisurf(XY_KL10, scores(KL10, choi=4), bubble = 4, main = "PCNM 4")
ordisurf(XY_KL10, scores(KL10, choi=2), bubble = 4, main = "PCNM 2")
dev.off()

print("LL 2008")
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/rare_LL2008_ordiSIG.pdf")
par(mfrow=c(1,2))
ordisurf(XY_LL, scores(LL, choi=5), bubble = 4, main = "PCNM 5")
ordisurf(XY_LL, scores(LL, choi=20), bubble = 4, main = "PCNM 20")
dev.off()

print("SU 2008")
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/plots/rare_SU2008_ordiSIG.pdf")
par(mfrow=c(1,1))
ordisurf(XY_SU, scores(SU, choi=6), bubble = 4, main = "PCNM 6")
dev.off()


