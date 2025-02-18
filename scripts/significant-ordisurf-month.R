#---
#title: "Significant Ordisurf Plots - By Month"
#author: "Alicia Halhed"
#date: "04/22/2020"
#---

# salloc --time=0-00:10:00 --mem=4G --account=def-cottenie
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
print("Set up (working directory, theme, and packages)")
# set working directory
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env")

# attach required packages
library(qiime2R)
library(vegan)
library(tidyverse)

print("Initiate functions for analysis")
# subset the XY's by grid and year
XY_month <- function(metadata, grid, year, month) {
  df1 <- subset(metadata, Grid == grid, 
                select = c("SampleID", "Location X", "Location Y", "Year", "Month"))
  df2 <- subset(df1, Year == year, 
                select = c("SampleID", "Location X", "Location Y", "Month"))
  df3 <- subset(df2, Month == month, 
                select = c("SampleID", "Location X", "Location Y"))
  df4 <- column_to_rownames(remove_rownames(df3), var = "SampleID")
  return(df4)
}

# get the data
print("Read in the metadata")
rs_q2_metadata <- read_q2metadata("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/input/RS_meta.tsv")

# start analysis
print("Starting initial data preparation")
print("Access XY data")
XY_AG <- lapply(c(3:8), FUN = XY_month, metadata = rs_q2_metadata, grid = "AG", year = 2008)
XY_CH <- lapply(c(4,5), FUN = XY_month, metadata = rs_q2_metadata, grid = "CH", year = 2008)
XY_JO <- XY_month(rs_q2_metadata, "JO", 2008, 5) # only one month
XY_KL8 <- lapply(c(3:8), FUN = XY_month, metadata = rs_q2_metadata, grid = "KL", year = 2008)
XY_KL9 <- lapply(c(3:8), FUN = XY_month, metadata = rs_q2_metadata, grid = "KL", year = 2009)
XY_KL10 <- lapply(c(2:6), FUN = XY_month, metadata = rs_q2_metadata, grid = "KL", year = 2010)
XY_LL <- lapply(c(4,5), FUN = XY_month, metadata = rs_q2_metadata, grid = "LL", year = 2008)
XY_SU <- XY_month(rs_q2_metadata, "SU", 2008, 5) # only one month

print("Computing Euclidean Distances")
d_AG <- lapply(XY_AG, dist)
d_CH <- lapply(XY_CH, dist)
d_JO <- dist(XY_JO) # only one month
d_KL8 <- lapply(XY_KL8, dist)
d_KL9 <- lapply(XY_KL9, dist)
d_KL10 <- lapply(XY_KL10, dist)
d_LL <- lapply(XY_LL, dist)
d_SU <- dist(XY_SU) # only one month


# Remove objects we're done with
print("Removing full metadata data frame")
rm(rs_q2_metadata)

# unweighted PCNM
print("Unweighted PCNM")
AG <- lapply(d_AG, pcnm)
CH <- lapply(d_CH, pcnm)
JO <- pcnm(d_JO) # only one month
KL8 <- lapply(d_KL8, pcnm)
KL9 <- lapply(d_KL9, pcnm)
KL10 <- lapply(d_KL10, pcnm)
LL <- lapply(d_LL, pcnm)
SU <- pcnm(d_SU) # only one month

# plot with ordisurf
print("Plotting significant PCNM axes with ordisurf")
# replace grid-year with values used in this script
print("Core")
print("AG 2008")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/core_AG2008_MordiSIG.pdf")
par(mfrow=c(1,1))
ordisurf(XY_AG[[3]], scores(AG[[3]], choi=14), bubble = 4, main = "PCNM 14")
par(mfrow=c(2,2))
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=15), bubble = 4, main = "PCNM 15")
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=4), bubble = 4, main = "PCNM 4")
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=27), bubble = 4, main = "PCNM 27")
par(mfrow=c(1,1))
ordisurf(XY_AG[[6]], scores(AG[[6]], choi=9), bubble = 4, main = "PCNM 9")
dev.off()

print("CH 2008")
print("No significant core axes")


print("JO 2008")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/core_JO2008_MordiSIG.pdf")
par(mfrow=c(1,1))
ordisurf(XY_JO, scores(JO, choi=1), bubble = 4, main = "PCNM 1")
dev.off()

print("KL 2008")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/core_KL2008_MordiSIG.pdf")
par(mfrow=c(1,1))
ordisurf(XY_KL8[[1]], scores(KL8[[1]], choi=15), bubble = 4, main = "PCNM 15")
par(mfrow=c(2,2))
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=24), bubble = 4, main = "PCNM 24")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=18), bubble = 4, main = "PCNM 18")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=19), bubble = 4, main = "PCNM 19")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=17), bubble = 4, main = "PCNM 17")
par(mfrow=c(2,2))
ordisurf(XY_KL8[[4]], scores(KL8[[4]], choi=15), bubble = 4, main = "PCNM 15")
ordisurf(XY_KL8[[4]], scores(KL8[[4]], choi=5), bubble = 4, main = "PCNM 5")
ordisurf(XY_KL8[[4]], scores(KL8[[4]], choi=2), bubble = 4, main = "PCNM 2")
dev.off()

print("KL 2009")
print("No significant core axes")

print("KL 2010")
print("No significant core axes")

print("LL 2008")
print("No significant core axes")

print("SU 2008")
print("No significant core axes")


print("Full")
print("AG 2008")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/AG2008_MordiSIG.pdf")
par(mfrow=c(2,2))
ordisurf(XY_AG[[3]], scores(AG[[3]], choi=8), bubble = 4, main = "PCNM 8")
ordisurf(XY_AG[[3]], scores(AG[[3]], choi=7), bubble = 4, main = "PCNM 7")
ordisurf(XY_AG[[3]], scores(AG[[3]], choi=2), bubble = 4, main = "PCNM 2")
par(mfrow=c(1,2))
ordisurf(XY_AG[[4]], scores(AG[[4]], choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_AG[[4]], scores(AG[[4]], choi=12), bubble = 4, main = "PCNM 12")
par(mfrow=c(2,3))
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=1), bubble = 4, main = "PCNM 1")
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=25), bubble = 4, main = "PCNM 25")
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=36), bubble = 4, main = "PCNM 36")
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=27), bubble = 4, main = "PCNM 27")
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=6), bubble = 4, main = "PCNM 6")
dev.off()

print("CH 2008")
print("No significant core axes")

print("JO 2008")
print("No significant core axes")


print("KL 2008")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/KL2008_MordiSIG.pdf")
par(mfrow=c(1,2))
ordisurf(XY_KL8[[1]], scores(KL8[[1]], choi=17), bubble = 4, main = "PCNM 17")
ordisurf(XY_KL8[[1]], scores(KL8[[1]], choi=2), bubble = 4, main = "PCNM 2")
par(mfrow=c(1,2))
ordisurf(XY_KL8[[2]], scores(KL8[[2]], choi=1), bubble = 4, main = "PCNM 1")
ordisurf(XY_KL8[[2]], scores(KL8[[2]], choi=13), bubble = 4, main = "PCNM 13")
par(mfrow=c(2,3))
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=1), bubble = 4, main = "PCNM 1")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=5), bubble = 4, main = "PCNM 5")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=6), bubble = 4, main = "PCNM 6")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=16), bubble = 4, main = "PCNM 16")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=17), bubble = 4, main = "PCNM 17")
par(mfrow=c(2,2))
ordisurf(XY_KL8[[4]], scores(KL8[[4]], choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_KL8[[4]], scores(KL8[[4]], choi=12), bubble = 4, main = "PCNM 12")
ordisurf(XY_KL8[[4]], scores(KL8[[4]], choi=3), bubble = 4, main = "PCNM 3")
par(mfrow=c(1,2))
ordisurf(XY_KL8[[5]], scores(KL8[[5]], choi=17), bubble = 4, main = "PCNM 17")
ordisurf(XY_KL8[[5]], scores(KL8[[5]], choi=3), bubble = 4, main = "PCNM 3")
par(mfrow=c(2,2))
ordisurf(XY_KL8[[6]], scores(KL8[[6]], choi=3), bubble = 4, main = "PCNM 3")
ordisurf(XY_KL8[[6]], scores(KL8[[6]], choi=6), bubble = 4, main = "PCNM 6")
ordisurf(XY_KL8[[6]], scores(KL8[[6]], choi=12), bubble = 4, main = "PCNM 12")
ordisurf(XY_KL8[[6]], scores(KL8[[6]], choi=2), bubble = 4, main = "PCNM 2")
dev.off()

print("KL 2009")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/KL2009_MordiSIG.pdf")
par(mfrow=c(1,2))
ordisurf(XY_KL9[[2]], scores(KL9[[2]], choi=1), bubble = 4, main = "PCNM 1")
ordisurf(XY_KL9[[2]], scores(KL9[[2]], choi=4), bubble = 4, main = "PCNM 4")
par(mfrow=c(1,2))
ordisurf(XY_KL9[[5]], scores(KL9[[5]], choi=3), bubble = 4, main = "PCNM 3")
ordisurf(XY_KL9[[5]], scores(KL9[[5]], choi=1), bubble = 4, main = "PCNM 1")
par(mfrow=c(1,1))
ordisurf(XY_KL9[[6]], scores(KL9[[6]], choi=2), bubble = 4, main = "PCNM 2")
dev.off()

print("KL 2010")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/KL2010_MordiSIG.pdf")
par(mfrow=c(1,1))
ordisurf(XY_KL10[[1]], scores(KL10[[1]], choi=3), bubble = 4, main = "PCNM 3")
par(mfrow=c(3,4))
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=8), bubble = 4, main = "PCNM 8")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=4), bubble = 4, main = "PCNM 4")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=5), bubble = 4, main = "PCNM 5")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=17), bubble = 4, main = "PCNM 17")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=3), bubble = 4, main = "PCNM 3")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=7), bubble = 4, main = "PCNM 7")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=10), bubble = 4, main = "PCNM 10")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=1), bubble = 4, main = "PCNM 1")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=9), bubble = 4, main = "PCNM 9")
par(mfrow=c(1,2))
ordisurf(XY_KL10[[3]], scores(KL10[[3]], choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_KL10[[3]], scores(KL10[[3]], choi=9), bubble = 4, main = "PCNM 9")
par(mfrow=c(1,1))
ordisurf(XY_KL10[[5]], scores(KL10[[5]], choi=3), bubble = 4, main = "PCNM 3")
dev.off()

print("LL 2008")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/LL2008_MordiSIG.pdf")
par(mfrow=c(2,1))
ordisurf(XY_LL[[1]], scores(LL[[1]], choi=11), bubble = 4, main = "PCNM 11")
ordisurf(XY_LL[[1]], scores(LL[[1]], choi=7), bubble = 4, main = "PCNM 7")
par(mfrow=c(1,1))
ordisurf(XY_LL[[2]], scores(LL[[2]], choi=9), bubble = 4, main = "PCNM 9")
dev.off()

print("SU 2008")
print("No significant core axes")

print("Rare")
print("AG 2008")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/rare_AG2008_MordiSIG.pdf")
par(mfrow=c(2,2))
ordisurf(XY_AG[[3]], scores(AG[[3]], choi=8), bubble = 4, main = "PCNM 8")
ordisurf(XY_AG[[3]], scores(AG[[3]], choi=7), bubble = 4, main = "PCNM 7")
ordisurf(XY_AG[[3]], scores(AG[[3]], choi=2), bubble = 4, main = "PCNM 2")
par(mfrow=c(1,2))
ordisurf(XY_AG[[4]], scores(AG[[4]], choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_AG[[4]], scores(AG[[4]], choi=12), bubble = 4, main = "PCNM 12")
par(mfrow=c(2,2))
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=25), bubble = 4, main = "PCNM 25")
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=36), bubble = 4, main = "PCNM 36")
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=27), bubble = 4, main = "PCNM 27")
ordisurf(XY_AG[[5]], scores(AG[[5]], choi=6), bubble = 4, main = "PCNM 6")
dev.off()

print("CH 2008")
print("No significant core axes")

print("JO 2008")
print("No significant core axes")

print("KL 2008")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/rare_KL2008_MordiSIG.pdf")
par(mfrow=c(1,2))
ordisurf(XY_KL8[[1]], scores(KL8[[1]], choi=17), bubble = 4, main = "PCNM 17")
ordisurf(XY_KL8[[1]], scores(KL8[[1]], choi=2), bubble = 4, main = "PCNM 2")
par(mfrow=c(1,2))
ordisurf(XY_KL8[[2]], scores(KL8[[2]], choi=13), bubble = 4, main = "PCNM 13")
ordisurf(XY_KL8[[2]], scores(KL8[[2]], choi=1), bubble = 4, main = "PCNM 1")
par(mfrow=c(2,3))
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=1), bubble = 4, main = "PCNM 1")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=5), bubble = 4, main = "PCNM 5")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=6), bubble = 4, main = "PCNM 6")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=16), bubble = 4, main = "PCNM 16")
ordisurf(XY_KL8[[3]], scores(KL8[[3]], choi=17), bubble = 4, main = "PCNM 17")
par(mfrow=c(2,2))
ordisurf(XY_KL8[[4]], scores(KL8[[4]], choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_KL8[[4]], scores(KL8[[4]], choi=12), bubble = 4, main = "PCNM 12")
ordisurf(XY_KL8[[4]], scores(KL8[[4]], choi=3), bubble = 4, main = "PCNM 3")
par(mfrow=c(1,2))
ordisurf(XY_KL8[[5]], scores(KL8[[5]], choi=17), bubble = 4, main = "PCNM 17")
ordisurf(XY_KL8[[5]], scores(KL8[[5]], choi=3), bubble = 4, main = "PCNM 3")
par(mfrow=c(2,2))
ordisurf(XY_KL8[[6]], scores(KL8[[6]], choi=3), bubble = 4, main = "PCNM 3")
ordisurf(XY_KL8[[6]], scores(KL8[[6]], choi=6), bubble = 4, main = "PCNM 6")
ordisurf(XY_KL8[[6]], scores(KL8[[6]], choi=12), bubble = 4, main = "PCNM 12")
dev.off()

print("KL 2009")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/rare_KL2009_MordiSIG.pdf")
par(mfrow=c(1,2))
ordisurf(XY_KL9[[2]], scores(KL9[[2]], choi=1), bubble = 4, main = "PCNM 1")
ordisurf(XY_KL9[[2]], scores(KL9[[2]], choi=4), bubble = 4, main = "PCNM 4")
par(mfrow=c(1,1))
ordisurf(XY_KL9[[3]], scores(KL9[[3]], choi=3), bubble = 4, main = "PCNM 3")
par(mfrow=c(1,2))
ordisurf(XY_KL9[[5]], scores(KL9[[5]], choi=3), bubble = 4, main = "PCNM 3")
ordisurf(XY_KL9[[5]], scores(KL9[[5]], choi=1), bubble = 4, main = "PCNM 1")
par(mfrow=c(1,1))
ordisurf(XY_KL9[[6]], scores(KL9[[6]], choi=2), bubble = 4, main = "PCNM 2")
dev.off()

print("KL 2010")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/rare_KL2010_MordiSIG.pdf")
par(mfrow=c(1,1))
ordisurf(XY_KL10[[1]], scores(KL10[[1]], choi=3), bubble = 4, main = "PCNM 3")
par(mfrow=c(3,4))
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=8), bubble = 4, main = "PCNM 8")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=4), bubble = 4, main = "PCNM 4")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=5), bubble = 4, main = "PCNM 5")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=17), bubble = 4, main = "PCNM 17")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=3), bubble = 4, main = "PCNM 3")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=7), bubble = 4, main = "PCNM 7")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=10), bubble = 4, main = "PCNM 10")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=1), bubble = 4, main = "PCNM 1")
ordisurf(XY_KL10[[2]], scores(KL10[[2]], choi=9), bubble = 4, main = "PCNM 9")
par(mfrow=c(1,2))
ordisurf(XY_KL10[[3]], scores(KL10[[3]], choi=2), bubble = 4, main = "PCNM 2")
ordisurf(XY_KL10[[3]], scores(KL10[[3]], choi=9), bubble = 4, main = "PCNM 9")
par(mfrow=c(1,1))
ordisurf(XY_KL10[[5]], scores(KL10[[5]], choi=3), bubble = 4, main = "PCNM 3")
dev.off()

print("LL 2008")
pdf(file = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelMicrobiome/plots/rare_LL2008_MordiSIG.pdf")
par(mfrow=c(2,1))
ordisurf(XY_LL[[1]], scores(LL[[1]], choi=11), bubble = 4, main = "PCNM 11")
ordisurf(XY_LL[[1]], scores(LL[[1]], choi=7), bubble = 4, main = "PCNM 7")
par(mfrow=c(1,1))
ordisurf(XY_LL[[2]], scores(LL[[2]], choi=9), bubble = 4, main = "PCNM 9")
ordisurf(XY_LL[[2]], scores(LL[[2]], choi=6), bubble = 4, main = "PCNM 6")
dev.off()

print("SU 2008")
print("No significant core axes")
