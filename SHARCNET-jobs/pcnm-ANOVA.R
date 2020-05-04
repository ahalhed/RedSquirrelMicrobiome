# load the file
load("/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrel.RData")

# attach required packages
library(tidyverse)
library(phyloseq)
library(vegan)

# phyloseq object (ps) previously read in (present in the environment)

# subset the XY's by grid and year
XY_year <- function(metadata, grid, year) {
  df1 <- subset(metadata, Grid == grid, 
                select = c("SampleID", "Location X", "Location Y", "Year"))
  df2 <- subset(df1, Year == year, 
                select = c("SampleID", "Location X", "Location Y"))
  df3 <- column_to_rownames(remove_rownames(df2), var = "SampleID")
  return(df3)
}
# GET the metadata
rs_q2_metadata <- read.table("/home/ahalhed/red-squirrel-w2020/input/RS_meta.tsv", sep="\t")
colnames(rs_q2_metadata) <- c("SampleID", "Grid", "Location X", "Location Y", "Sex", "Age", "Month", "Season", "Year", "Squirrel.ID", "SireID", "DamID", "CollectionDate", "FoodSupplement", "BirthYear", "Location", "Date")

# just to subset the metadata
XY_08_KL <- XY_year(rs_q2_metadata, "KL", 2008)

# get OTUs (aka a community matrix)
comm_08_KL <- otu_table(ps) %>% as.matrix %>% 
  as.data.frame %>% 
  t %>% as.data.frame %>% 
  subset(., row.names(.) %in% rownames(XY_08_KL)) %>%
  .[ rowSums(.)>0, ]
# sample ID's are rownames

# get the metadata subset
meta_08_KL <- rs_q2_metadata %>% 
  subset(., SampleID %in% rownames(XY_08_KL)) %>%
  remove_rownames() %>% column_to_rownames(var = "SampleID") %>%
  select_if(~ !any(is.na(.))) %>% 
  select(Sex, Age, Month, Season, CollectionDate, BirthYear)

# pcnm
pcnm_08_KL <- pcnm(dist(XY_08_KL))

# Test fraction [a+b], total environment, using RDA:
print("ANOVA for abFrac")
abFrac <- rda(decostand(comm_08_KL, "hel") ~ ., meta_08_KL)
anova(abFrac, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
print("abFrac adjusted R Squared")
RsquareAdj(abFrac)


# Test fraction [a] using partial RDA:
print("ANOVA for aFrac")
aFrac <- rda(decostand(comm_08_KL, "hel") ~ . + Condition(scores(pcnm_08_KL)), data = meta_08_KL)
anova(aFrac, step=200, perm.max=1000)
# RsquareAdj gives the same result as component [a] of varpart
print("aFrac adjusted R Squared")
RsquareAdj(aFrac)