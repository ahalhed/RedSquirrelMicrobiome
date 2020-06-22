#---
#title: "Plots4Ch2"
#author: "Alicia Halhed"
#date: "22/06/2020"
#output: html_document
#---
# set working directory to output the plots into
# run on graham cluster
# salloc --time=0-01:30:00 --mem=64G --account=def-cottenie
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
setwd("/home/ahalhed/red-squirrel/R-env/RedSquirrelSpatial/plots")
# attaching required packages for full analysis
# qiime2R to create phyloseq object
library(qiime2R)
library(phyloseq)
library(vegan)
library(lubridate)
library(ggpubr)
library(tidyverse)
# set theme for ggplots
theme_set(theme_bw())

# for accessing non-zero occurreneces of each OTU
nonzero <- function(x) sum(x!=0)


# subset the XY's by grid and year
# relies on tidyverse
XY_year <- function(metadata, grid, year) {
  # metadata is a data frame containing environmetal data
  # grid is the two letter grid label as found in the environmental data
  # year is the collection of interest
  rn <- rownames_to_column(metadata, var = "SampleID")
  df1 <- subset(rn, Grid == grid, 
                select = c("SampleID", "Location.X", "Location.Y", "Year"))
  df2 <- subset(df1, Year == year, 
                select = c("SampleID", "Location.X", "Location.Y"))
  df3 <- column_to_rownames(remove_rownames(df2), var = "SampleID")
  return(df3)
}


# subset an OTU table for grid/year based on XY data
# relies on tidyverse, phyloseq
Comm <- function(phy, XY) {
  # phy is a phyloseq object
  # XY is a two column data frame containing XY locations
  trans <- otu_table(phy) %>% 
    as.matrix %>% as.data.frame %>% 
    t %>% as.data.frame
  sub <- subset(trans, rownames(trans) %in% rownames(XY))
  nz <- sub[ , colSums(sub)>0 ]
  return(nz)
}


# compute PCNM scores and export as data frame
# relies on vegan
PCNM_df <- function(XY) {
  # XY is a two column data frame containing XY locations
  eDIST <- dist(XY)
  UWpcnm <- pcnm(eDIST)
  pcnm_df <- as.data.frame(scores(UWpcnm))
  return(pcnm_df)
}


# calculate jaccard distance from an community matrix formatted OTU table
# relies on tidyverse, vegan
jac <- function(OTU, met) {
  # OTU is the community matrix containing the OTUs
  # met is a data frame containing the environmental metadata
  m <- rownames_to_column(met, var = "SampleID")
  df <- OTU %>% data.matrix(rownames.force = T) %>%
    vegdist(method = "jaccard") %>% 
    as.matrix %>% as.data.frame %>%
    rownames_to_column("sampleid1") %>%
    pivot_longer(-sampleid1, names_to = "sampleid2", values_to = "Jaccard") %>%
    unique() %>% left_join(., m, 
                           by = c("sampleid1" = "SampleID")) %>%
    # .x for sample 1, .y for sample 2
    left_join(., m, by = c("sampleid2" = "SampleID"))
  return(df)
}


gr <- function(df){
  # df is a data frame containing metadata and jaccard distances
  # create a new grid column
  df$Grid <- ifelse(df$Grid.x == df$Grid.y, "Same Grid", "Different Grid")
  # add location information
  df$Location <- ifelse(df$sampleid1 == df$sampleid2, "Same Sample",
                        ifelse(df$`Grid.x` != df$`Grid.y`, "Different Grid",
                               ifelse(df$`Squirrel.ID.x` == df$`Squirrel.ID.y` &
                                        df$`Location.X.x` == df$`Location.X.y` &
                                        df$`Location.Y.x` == df$`Location.Y.y`,
                                      "Same Squirrel, Same Location",
                                      ifelse(df$Squirrel.ID.x != df$Squirrel.ID.y &
                                               df$`Location.X.x` == df$`Location.X.y` &
                                               df$`Location.Y.x` == df$`Location.Y.y`,
                                             "Different Squirrel, Same Location",
                                             ifelse(df$Squirrel.ID.x != df$Squirrel.ID.y &
                                                      df$`Location.X.x` != df$`Location.X.y` |
                                                      df$`Location.Y.x` != df$`Location.Y.y`,
                                                    "Different Squirrel, Different Location",
                                                    ifelse(df$Squirrel.ID.x == df$Squirrel.ID.y &
                                                             df$`Location.X.x` != df$`Location.X.y` |
                                                             df$`Location.Y.x` != df$`Location.Y.y`,
                                                           "Same Squirrel, Different Location", NA))))))
  # add date information
  df$CollectionDate <- ifelse(df$CollectionDate.x == df$CollectionDate.y, "Same Date", "Different Date")
  return(df)
}


# add interval column
INT <- function(df) {
  # df is a dataframe containing jaccard distances and specific comparison information
  # first find the date range with interval
  df1 <- interval(df$`CollectionDate.x`, df$`CollectionDate.y`)
  # find length of the interval
  df2 <- as.period(df1, unit = "days")
  # extract number of days only
  df3 <- gsub("([0-9]+).*$", "\\1", df2)
  df4 <- df3 %>% as.numeric %>% abs
  return(df4)
}


## Getting the data

# directory for these files
# setwd("/home/ahalhed/red-squirrel")
# needs 64G RAM to build
ps <- qza_to_phyloseq(features = "/home/ahalhed/red-squirrel/filtered-table.qza",
                      tree = "/home/ahalhed/red-squirrel/trees/rooted_tree.qza",
                      taxonomy = "/home/ahalhed/red-squirrel/taxonomy/GG-taxonomy.qza",
                      metadata = "/home/ahalhed/red-squirrel/input/RS_meta.tsv")
# load(file = "ps.RData")

# extract the metadata as a data frame
# based on code from meta function in the microbiome package
meta <- as(sample_data(ps), "data.frame")
rownames(meta) <- sample_names(ps)
# not using the meta function avoid loading the whole library
# just to use a single function one time

# transpose and convert to data frame
OTU_full <- otu_table(ps) %>% as.matrix %>% 
  as.data.frame %>% 
  t %>% as.data.frame 

## Core Microbial Community
# transform the counts to relative abundance
RA_full <- transform_sample_counts(ps, function(x) x/sum(x)) %>% 
  otu_table %>% as.matrix %>% as.data.frame %>% 
  t %>% as.data.frame

#Presence in 95% of samples; relative abundace > 0.1%
# 95% threshold
# minimum number of samples that an OTU must occur in 
# to be considered as part of the "core" microbiome
print("Calculating 95% threshold")
th <- 0.95*nrow(RA_full)
th

# Access the number of non-zero occurrences of each OTU
# see nonzero function in the Functions section
nz <- plyr::numcolwise(nonzero)(RA_full)

# total relative abundance
# average relative abundance across all 909 samples
RA <- RA_full %>% colSums() %>% as.data.frame
RA_001 <- RA %>%
  rownames_to_column(var = "OTU") %>% 
  rename("ColSum" = ".") %>% 
  mutate(RelativeAbundance = .$ColSum/nrow(OTU_full)) %>%
  filter(RelativeAbundance > 0.001)

# extract the OTUs that fall above the 95% threshold
# based on the number of non-zero occurrences of each OTU
OTU95 <- t(nz) %>% as.data.frame %>% 
  rownames_to_column(var = "OTU") %>% 
  rename("Occurrences" = "V1") %>%
  subset(., Occurrences > th)

# combine to ensure that those in 95% of the samples
# have a total relative abundance of 0.1%
cOTU <- inner_join(OTU95, RA_001)
# how many OTUs meet this
nrow(cOTU)

# Subset the OTU table to find core and rare OTUs
OTU_core <- OTU_full[, cOTU$OTU]
OTU_rare <- OTU_full %>% 
  select(-one_of(cOTU$OTU))

# Removing objects that are no longer needed
rm(cOTU, OTU95, RA_001, RA)


## Figure 1 - Core Community Cutoff
# create mini dataframes
# with the number of non zero occurrences
plo1 <- data.frame(t(nz))
# and the relative abundance
plo2 <- data.frame(colSums(RA_full[,]))
# then merge them together
plo <- merge(plo1,plo2,by="row.names",all.x=TRUE) %>%
  rename("OTU Label" = "Row.names" , 
         "Total Non-Zero Occurrences" = "t.nz.", 
         "Total Relative Abundance" = "colSums.RA_full.....") %>%
  # with a column for percentage
  mutate("Presence (% of Samples)" = 100*`Total Non-Zero Occurrences`/909, 
         "Relative Abundance (% of OTUs)" = 100*`Total Relative Abundance`/909)
# remove the mini data frames we're done with
rm(plo1, plo2)


# This plot shows the fraction of the OTUs included in the core microbiome
# create the framework for the plot
fig1 <- ggplot(plo, aes(x = `Presence (% of Samples)`, y = `Relative Abundance (% of OTUs)`)) + 
  geom_point() + 
  # adjust axes to remove buffer around ploting area
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) +
  # log transform the axis
  scale_y_log10() + 
  # non-text labelling for core area
  geom_vline(xintercept = 95, linetype="dashed", colour = "red") +
  geom_hline(yintercept=0.1, linetype="dashed", colour = "red") +
  annotate("rect", ymin = 0.1, ymax = 6, xmin = 95, xmax = 100, 
           alpha = 0.1, fill = "blue") + 
  # text annotations for different quadrants
  annotate("text", x = 101, y = 0.9, label = "Core OTUs", angle = -90) +
  annotate("text", x = 25, y = 5, label = "High abundance, low presence") + 
  annotate("text", x = 25, y = 0.00005, label = "Low abundance, low presence") +
  annotate("text", x = 101, y = 0.0005, label = "Low abundance, high presence", angle = -90)

# export plot 1 to a file
pdf(file = "figure1.pdf", width = 14)
fig1
dev.off()


## Figure 2 - spatial pattern example
# See code from _____ for creating the ordisurf plots (need to pick a month to include).

## Figure 3 - Adjusted R2

# read in the data
# values taken from output files from monthly PCNM analyses
adj <- read.csv("/home/ahalhed/red-squirrel/R-env/data/FullAdjR2.csv",
                na.strings = ("NS")) # this is the advantage of read.csv instead of read_csv
# pivot the data of interest into long format
adj <- adj %>% 
  # take only the spatial value rows (don't need environmental here)
  filter(VariableType == "Spatial") %>%
  # take the grid/date labels and adjusted R squared column for core/peripheral
  select(Grid, Year, Month, CoreR2Adj, RareR2Adj) %>%
  # rename columns (don't need to know the values are R2Adj)
  rename("Core" = "CoreR2Adj", "Peripheral" = "RareR2Adj") %>%
  # pivot to long format
  pivot_longer(-c(Grid, Year, Month), names_to = "Community", values_to = "Adjusted R2 Value")

# exporting figure 3
pdf("figure3.pdf", width = 10)
ggplot(adj, aes(Month, `Adjusted R2 Value`, colour = Community)) +
  geom_smooth(method = "lm") + geom_point(aes(shape = as.character(Year))) + 
  viridis::scale_color_viridis(discrete=TRUE) +
  labs(y = expression(paste("Adjusted R"^"2")), shape = "Collection Year") + 
  annotate("text", x = 3.75, y = 0.25, 
           label = "Total OTUs & Significant Spatial Patterns\n Peripheral Community: 541908, 23\n Core Community: 14, 11")
dev.off()


## Figure 4 - LOESS regression
# calculate jaccard distances
# these are relative slow running, so saving them in case the figures need to be modified
core_jaccard <- jac(OTU_core, meta)
write.table(core_jaccard, file='/home/ahalhed/red-squirrel/R-env/data/core-jaccard.tsv', quote=FALSE, sep='\t', row.names = F)
# core_jaccard <- read.table("/home/ahalhed/red-squirrel/R-env/data/core-jaccard.tsv", sep = "\t", header = T)
rare_jaccard <- jac(OTU_rare, meta)
write.table(rare_jaccard, file='/home/ahalhed/red-squirrel/R-env/data/rare-jaccard.tsv', quote=FALSE, sep='\t', row.names = F)
# rare_jaccard <- read.table("/home/ahalhed/red-squirrel/R-env/data/rare-jaccard.tsv", sep = "\t", header = T)

# mutate jaccard distance df to include date columns
core_jac <- core_jaccard %>%
  mutate(CollectionDate.x = as_date(.$CollectionDate.x),
         CollectionDate.y = as_date(.$CollectionDate.y))
rare_jac <- rare_jaccard %>%
  mutate(CollectionDate.x = as_date(.$CollectionDate.x),
         CollectionDate.y = as_date(.$CollectionDate.y))
# add the interval column between collection dates and locations
# core
core_jac <- gr(core_jac)
core_jac$int <- INT(core_jac)
# peripheral
rare_jac <- gr(rare_jac)
rare_jac$int <- INT(rare_jac)

# break up different groups
# core
core_dsSL <- filter(core_jac, Location %in% "Different Squirrel, Same Location")
core_dsDL <- filter(core_jac, Location %in% "Different Squirrel, Different Location")
core_ssSL <- filter(core_jac, Location %in% "Same Squirrel, Same Location")
core_ssDL <- filter(core_jac, Location %in% "Same Squirrel, Different Location")
# peripheral
rare_dsSL <- filter(rare_jac, Location %in% "Different Squirrel, Same Location")
rare_dsDL <- filter(rare_jac, Location %in%"Different Squirrel, Different Location")
rare_ssSL <- filter(rare_jac, Location %in% "Same Squirrel, Same Location")
rare_ssDL <- filter(rare_jac, Location %in% "Same Squirrel, Different Location")

# bind the four together (10% of different location, different location)
linesC <- rbind(core_ssDL, core_dsSL) %>%
  rbind(., core_ssSL) %>% rbind(., core_dsDL %>% sample_frac(.1))
linesR <- rbind(rare_ssDL, rare_dsSL) %>%
  rbind(., rare_ssSL) %>% rbind(., rare_dsDL %>% sample_frac(.1))
# within each of the years
linesYR <- linesR[which(linesR$Year.x == linesR$Year.y),]
linesYC <- linesC[which(linesC$Year.x == linesC$Year.y),]
# remove objects no longer needed
rm(linesC, linesR,
   core_dsSL, core_dsDL, core_ssSL, core_ssDL,
   rare_ssDL, rare_ssSL, rare_dsDL, rare_dsSL)

# peripheral by year
rareYP <- ggplot(linesYR, aes(x = int, y = Jaccard, color = Location)) +
  geom_smooth(method='loess', formula= y~x) + 
  labs(x = "Days between Sample Collection", y = "Jaccard Distance",
       color = "Samples Being Compared") + 
  viridis::scale_color_viridis(discrete=TRUE) +
  ggtitle("Peripheral Microbiome Within Collection Year")
rareYP

# core by year
coreYP <- ggplot(linesYC, aes(x = int, y = Jaccard, color = Location)) +
  geom_smooth(method='loess', formula= y~x) + 
  labs(x = "Days between Sample Collection", y = "Jaccard Distance",
       color = "Samples Being Compared") + 
  ggtitle("Core Microbiome Within Collection Year") +
  viridis::scale_color_viridis(discrete=TRUE)
coreYP

# export figure 2
pdf("figure2.pdf", width = 14)
ggarrange(coreYP, rareYP + ylim(0.4,1), common.legend = T, legend = "right")
dev.off()