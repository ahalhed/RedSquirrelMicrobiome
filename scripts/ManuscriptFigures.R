#---
#title: "Plots4Ch2"
#author: "Alicia Halhed"
#date: "22/06/2020"
#output: html_document
#---
# set working directory to output the plots into
# run on graham cluster interactively
# salloc --time=0-00:30:00 --mem=8G --account=def-cottenie
# module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
setwd("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/")
# attaching required packages for full analysis
# qiime2R to create phyloseq object
library(qiime2R)
library(phyloseq)
library(vegan)
library(zCompositions)
#devtools::install_github('ggloor/CoDaSeq/CoDaSeq')
library(CoDaSeq)
library(lubridate)
library(ggpubr)
library(tidyverse)
# set theme for ggplots
theme_set(theme_bw())

# calculate euclidean dissimilarity from an community matrix formatted OTU table
# relies on tidyverse, vegan
dis <- function(OTU, met) {
  # OTU is the community matrix containing the OTUs
  # met is a data frame containing the environmental metadata
  m <- rownames_to_column(met, var = "SampleID")
  df <- OTU %>% data.matrix(rownames.force = T) %>%
    vegdist(method = "euclidean") %>% 
    as.matrix %>% as.data.frame %>%
    rownames_to_column("sampleid1") %>%
    pivot_longer(-sampleid1, names_to = "sampleid2", values_to = "EucDis") %>%
    unique() %>% left_join(., m, 
                           by = c("sampleid1" = "SampleID")) %>%
    # .x for sample 1, .y for sample 2
    left_join(., m, by = c("sampleid2" = "SampleID"))
  return(df)
}


gr <- function(df){
  # df is a data frame containing metadata and Bray-Curtis distances
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


## get the data
print("Read in the Data")
print("Building phyloseq object")
ps <- qza_to_phyloseq(features = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/filtered-table-10.qza",
                      tree = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/trees/rooted_tree.qza",
                      metadata = "/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/input/RS_meta.tsv") %>%
  phyloseq(otu_table(t(otu_table(.)), taxa_are_rows = F), phy_tree(.), sample_data(.))
# taxonomy only
# phyloseq didn't like the initial labelling from SILVA
tax <- read_qza("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/taxonomy/SILVA-taxonomy-10.qza")$data %>%
  column_to_rownames(var = "Feature.ID")
# so I'm correcting that here to make it work
tax$Taxon <- tax$Taxon %>%
  str_replace_all("D_0__", "k__") %>%
  str_replace_all("D_1__", "p__") %>%
  str_replace_all("D_2__", "c__") %>%
  str_replace_all("D_3__", "o__") %>%
  str_replace_all("D_4__", "f__") %>%
  str_replace_all("D_5__", "g__") %>%
  str_replace_all("D_6__", "s__")
tax1 <- tax %>% 
  mutate("Kingdom" = word(.$Taxon, 1, sep = ";"), #k__
         "Phylum" = word(.$Taxon, 2, sep = ";"), #p__
         "Class" = word(.$Taxon, 3, sep = ";"), #c__
         "Order" = word(.$Taxon, 4, sep = ";"), #o__
         "Family" = word(.$Taxon, 5, sep = ";"), #f__
         "Genus" = word(.$Taxon, 6, sep = ";"), #g__
         "Species" = word(.$Taxon, 7, sep = ";")) %>% #s__
  select(Kingdom, Phylum, Class, Order, Family, Genus, Species)
# making the modified Q2 artifact into a phyloseq tax table
tax2 <- tax_table(tax1)
taxa_names(tax2) <- rownames(tax)
# combining the taxonomy with the data already in the phyloseq object
ps <- merge_phyloseq(ps, tax2)

# based on the meta function from the microbiome package
# I don't want to load a whole package for one function
print("Read in the metadata")
meta <- as(sample_data(ps), "data.frame")
rownames(meta) <- sample_names(ps)

print("Full OTU table")
print("Aitchison transformation")
# rows are OTUs, then transposed to OTUs as column
# impute the OTU table
OTUimp <- cmultRepl(otu_table(ps), label=0, method="CZM") # all OTUs
# compute the aitchison values
OTU_full <- codaSeq.clr(OTUimp) %>% as.data.frame

## Core and rare divide
print("Finding core microbiome")
print("Extract 95% Occupancy from BC Similarity Core")
# read in occupancy/abundance information
occ_abun <- read.csv("./data/core.csv")
cOTU <- occ_abun %>%
  # get the OTUs identified as core contributors to beta diversity
  # and greater than 95% occupancy (confirmed core)
  .[which(.$Community == "Confirmed Core"),]
# make the new data frames
print("Subset the OTU table to find core and rare OTUs")
OTU_core <- select(OTU_full, one_of(cOTU$otu))
OTU_rare <- select(OTU_full, -one_of(cOTU$otu))

# Removing objects that are no longer needed
rm(tax, tax1, tax2)


## Figure 1 - Core Community Cutoff
# This plot shows the fraction of the OTUs included in the core microbiome
# create the framework for the plot
fig1 <- ggplot(occ_abun, aes(y = otu_occ, x = otu_rel, color = Community, shape = Community)) + 
  geom_point() +
  # log transform the x axis
  scale_x_log10() +
  # set discrete viridis colour scheme
  scale_colour_viridis_d() + 
  # add axis labels
  labs(x = "Mean Relative Abundance of Each OTU (log10)", y = "Occupancy (Proportion of Samples)")

# export plot 1 to a file
pdf(file = "./plots/figure1.pdf", width = 14)
fig1
dev.off()


## Figure 2 - spatial pattern example
# See code from _____ for creating the ordisurf plots (need to pick a month to include).

## Figure 3 - Adjusted R2

# read in the data
# values taken from output files from monthly PCNM analyses
adj <- read.csv("./data/A-FullAdjR2.csv",
                na.strings = ("NS")) # this is the advantage of read.csv instead of read_csv
# pivot the data of interest into long format
adj <- adj %>% 
  # take only the spatial value rows (don't need environmental here)
  filter(VariableType == "Spatial") %>%
  # take the grid/date labels and adjusted R squared column for core/peripheral
  select(Grid, Year, Month, CoreR2Adj, RareR2Adj, FullR2Adj) %>%
  # rename columns (don't need to know the values are R2Adj)
  rename("Core" = "CoreR2Adj", "Rare" = "RareR2Adj", "Full" = "FullR2Adj") %>%
  # pivot to long format
  pivot_longer(-c(Grid, Year, Month), names_to = "Community", values_to = "Adjusted R2 Value")

fig3 <- ggplot(adj, aes(Month, `Adjusted R2 Value`, colour = Community)) +
  geom_smooth(method = "lm", aes(linetype = Community)) + 
  geom_jitter(aes(shape = as.character(Year))) + 
  scale_color_viridis_d() +
  labs(y = expression(paste("Adjusted R"^"2")), shape = "Collection Year") 
# exporting figure 3
pdf("./plots/figure3.pdf", width = 10)
fig3 #+ stat_regline_equation(label.y = c(0.21, 0.095, 0.075))
dev.off()

# is there a significant difference in the R2adj values based on the month and community of origin?
lm(`Adjusted R2 Value` ~ Community*Month, data = adj) %>% anova

## Figure 4 - LOESS regression
# calculate bray-curtis dissimilarity
# saving in case the figures need to be modified (gzip after)
core_dis <- dis(OTU_core, meta)
write.table(core_dis, file='./data/core-dis.tsv', quote=FALSE, sep='\t', row.names = F)
# core_dis <- read.table("/home/ahalhed/red-squirrel/R-env/data/core-dis.tsv.gz", sep = "\t", header = T)
rare_dis <- dis(OTU_rare, meta)
write.table(rare_dis, file='./data/rare-dis.tsv', quote=FALSE, sep='\t', row.names = F)
# rare_dis <- read.table("/home/ahalhed/red-squirrel/R-env/data/rare-dis.tsv.gz", sep = "\t", header = T)
full_dis <- dis(OTU_full, meta)
write.table(full_dis, file='./data/full-dis.tsv', quote=FALSE, sep='\t', row.names = F)
# full_dis <- read.table("/home/ahalhed/red-squirrel/R-env/data/full-dis.tsv", sep = "\t", header = T)

# mutate jaccard distance df to include date columns
core_dis <- core_dis %>%
  mutate(CollectionDate.x = as_date(.$CollectionDate.x),
         CollectionDate.y = as_date(.$CollectionDate.y))
rare_dis <- rare_dis %>%
  mutate(CollectionDate.x = as_date(.$CollectionDate.x),
         CollectionDate.y = as_date(.$CollectionDate.y))
full_dis <- full_dis %>%
  mutate(CollectionDate.x = as_date(.$CollectionDate.x),
         CollectionDate.y = as_date(.$CollectionDate.y))
# add the interval column between collection dates and locations
# core
core_dis <- gr(core_dis)
core_dis$int <- INT(core_dis)
# rare
rare_dis <- gr(rare_dis)
rare_dis$int <- INT(rare_dis)
# full
full_dis <- gr(full_dis)
full_dis$int <- INT(full_dis)
# Warning message:
# tz(): Don't know how to compute timezone for object of class factor; returning "UTC". This warning will become an error in the next major version of lubridate. 

# break up different groups
# core
core_dsSL <- filter(core_dis, Location %in% "Different Squirrel, Same Location")
core_dsDL <- filter(core_dis, Location %in% "Different Squirrel, Different Location")
core_ssSL <- filter(core_dis, Location %in% "Same Squirrel, Same Location")
core_ssDL <- filter(core_dis, Location %in% "Same Squirrel, Different Location")
# rare
rare_dsSL <- filter(rare_dis, Location %in% "Different Squirrel, Same Location")
rare_dsDL <- filter(rare_dis, Location %in% "Different Squirrel, Different Location")
rare_ssSL <- filter(rare_dis, Location %in% "Same Squirrel, Same Location")
rare_ssDL <- filter(rare_dis, Location %in% "Same Squirrel, Different Location")
# full
full_dsSL <- filter(full_dis, Location %in% "Different Squirrel, Same Location")
full_dsDL <- filter(full_dis, Location %in% "Different Squirrel, Different Location")
full_ssSL <- filter(full_dis, Location %in% "Same Squirrel, Same Location")
full_ssDL <- filter(full_dis, Location %in% "Same Squirrel, Different Location")

# bind the four together (10% of different location, different location)
linesC <- rbind(core_ssDL, core_dsSL) %>%
  rbind(., core_ssSL) %>% rbind(., core_dsDL %>% sample_frac(.1))
linesR <- rbind(rare_ssDL, rare_dsSL) %>%
  rbind(., rare_ssSL) %>% rbind(., rare_dsDL %>% sample_frac(.1))
linesF <- rbind(full_ssDL, full_dsSL) %>%
  rbind(., full_ssSL) %>% rbind(., full_dsDL %>% sample_frac(.1))
# within each of the years
linesYR <- linesR[which(linesR$Year.x == linesR$Year.y),]
linesYC <- linesC[which(linesC$Year.x == linesC$Year.y),]
linesYF <- linesF[which(linesF$Year.x == linesF$Year.y),]
# remove objects no longer needed
rm(linesC, linesR, linesF,
   core_dsSL, core_dsDL, core_ssSL, core_ssDL,
   rare_ssDL, rare_ssSL, rare_dsDL, rare_dsSL,
   full_ssDL, full_ssSL, full_dsDL, full_dsSL)

# peripheral by year
rareYP <- ggplot(linesYR, aes(x = int, y = EucDis, color = Location, linetype = Location)) +
  geom_smooth(method='loess', formula= y~x) + 
  labs(x = "Days between Sample Collection", y = "Euclidean Distance",
       color = "Samples Being Compared") + 
  scale_colour_viridis_d() +
  scale_linetype_manual("Samples Being Compared", values=c(1,2,4,3))

# core by year
coreYP <- ggplot(linesYC, aes(x = int, y = EucDis, color = Location, linetype = Location)) +
  geom_smooth(method='loess', formula= y~x) + 
  labs(x = "Days between Sample Collection", y = "Euclidean Distance",
       color = "Samples Being Compared") + 
  scale_colour_viridis_d() +
  scale_linetype_manual("Samples Being Compared", values=c(1,2,4,3))

# core by year
fullYP <- ggplot(linesYF, aes(x = int, y = EucDis, color = Location, linetype = Location)) +
  geom_smooth(method='loess', formula= y~x) + 
  labs(x = "Days between Sample Collection", y = "Euclidean Distance",
       color = "Samples Being Compared") + 
  ggtitle("Full Microbial Community") +
  scale_colour_viridis_d() +
  scale_linetype_manual("Samples Being Compared", values=c(1,2,4,3))

# export figure 4
pdf("./plots/figure4.pdf", height = 10, width = 12)
ggarrange(coreYP, rareYP, labels = c("A", "B"),
          nrow=2, common.legend = T)
dev.off()

# putting full in a supplemental figure
pdf("./plots/supp4full.pdf")
fullYP
dev.off()
