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
abFrac <- rda(decostand(comm_08_KL, "hel") ~ ., meta_08_KL)

# spatial variables
#bcFrac <- rda(decostand(comm_08_KL, "hel") ~ ., scores(pcnm_08_KL)) # Full model, oops
print("Running RDA's")
rs_pcnm <- as.data.frame(scores(pcnm_08_KL))
bcFrac <- rda(decostand(comm_08_KL, "hel") ~ ., rs_pcnm) # Full model
bcFrac0 <- rda(decostand(comm_08_KL, "hel") ~ 1, rs_pcnm) # Reduced model
step.space <- ordiR2step(bcFrac0, scope = formula(bcFrac))
print("Summary of selection process")
step.space$anova
#                  R2.adj Df     AIC      F Pr(>F)   
#  + PCNM2         0.0026787  1 -105.77 1.6419  0.008 **
#  + PCNM33        0.0045491  1 -105.23 1.4472  0.012 * 
#  <All variables> 0.0049302  

# save the plot
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/step_space.pdf")
plot(step.space)
dev.off() # close the graphics device
