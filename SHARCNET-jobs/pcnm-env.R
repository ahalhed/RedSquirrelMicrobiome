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

# env variables
print("Test of fraction [a+b], total environment, using RDA")
abFrac # Full model
#Call: rda(formula = decostand(comm_08_KL, "hel") ~ Sex + Age + Month + Season + CollectionDate + BirthYear, data = meta_08_KL)
#               Inertia Proportion Rank
#Total          0.6400     1.0000     
#Constrained    0.2531     0.3954   91
#Unconstrained  0.3869     0.6046  148
#Inertia is variance 
#Some constraints were aliased because they were collinear (redundant)

abFrac0 <- rda(decostand(comm_08_KL, "hel") ~ 1, meta_08_KL) # Reduced model
#Step: R2.adj= 0 
#Call: decostand(comm_08_KL, "hel") ~ 1 
#                R2.adjusted
#+ CollectionDate 0.0239596874
#<All variables>  0.0236843856
#+ Season         0.0118833967
#+ Month          0.0114943283
#+ Age            0.0017636551
#+ BirthYear      0.0017636551
#+ Sex            0.0006824704
#<none>           0.0000000000

# Here is where the magic happens, but almost automatically!
step.env <- ordiR2step(abFrac0, scope = formula(abFrac))
str(step.env) #mumbo jumbo in pcnm-env-30055778.out 

print("step.env")
step.env # but it is actually just an rda model, with the final model predictor variables
#              Inertia Rank
#Total            0.64     
#Unconstrained    0.64  239
#Inertia is variance 


print("ANOVA of step.env")
anova(step.env) # so you can do the same stuff as before


# save the plot
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/step_env.pdf")
plot(step.env)
dev.off() # close the graphics device

print("Summary of selection process")
step.env$anova # and if you are interested, this is a summary of the selection process
# Note, show how to access the RStudio object explorer
# not sure if this works on the server
step.env[["terminfo"]][["ordered"]]