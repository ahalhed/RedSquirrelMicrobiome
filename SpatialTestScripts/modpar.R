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

# just to subset the metadata and get euclidean distance
XY_08_KL <- XY_year(rs_q2_metadata, "KL", 2008)
e_08_KL <- dist(XY_08_KL)

# get the metadata subset
meta_08_KL <- rs_q2_metadata %>% 
  subset(., SampleID %in% rownames(XY_08_KL)) %>%
  remove_rownames() %>% column_to_rownames(var = "SampleID") %>%
  select_if(~ !any(is.na(.))) %>% 
  select(Sex, Age, Month, Season, CollectionDate, BirthYear)

comm_08_KL <- otu_table(ps) %>% as.matrix %>% 
    as.data.frame %>% 
    t %>% as.data.frame %>% 
    subset(., rownames(.) %in% rownames(XY_08_KL)) %>%
    .[ rowSums(.)>0, ]

# unweighted PCNM
pcnm_08_KL <- pcnm(e_08_KL)

abFrac <- rda(decostand(comm_08_KL, "hel") ~ ., meta_08_KL) # full model
abFrac0 <- rda(decostand(comm_08_KL, "hel") ~ 1, meta_08_KL) # Reduced model
# Here is where the magic happens, but almost automatically!
step.env <- ordiR2step(abFrac0, scope = formula(abFrac))

rs_pcnm <- as.data.frame(scores(pcnm_08_KL))
bcFrac <- rda(decostand(comm_08_KL, "hel") ~ ., rs_pcnm) # Full model
bcFrac0 <- rda(decostand(comm_08_KL, "hel") ~ 1, rs_pcnm) # Reduced model
step.space <- ordiR2step(bcFrac0, scope = formula(bcFrac))

print("variation decomposition with parsimonious variables")
mod.pars <- varpart(comm_08_KL, ~ ., 
                    rs_pcnm[, names(step.space$terminfo$ordered)], 
                    data = meta_08_KL[, names(step.env$terminfo$ordered)],
                    transfo = "hel")
mod.pars
# save the plot
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/modpars_08KL.pdf")
plot(mod.pars)
dev.off() # close the graphics device
