# R script starts here
load("/home/ahalhed/red-squirrel-w2020/R-env/rs-dbRDA.RData")
library(phyloseq)
library(vegan)
library(tidyverse)

# ordinate for bray-curtis
rs_ordi_bc <- ordinate(ps, "CAP", "bray", ~Year+Grid+Age+Sex)
# ordinate for jaccard
rs_ordi_j <- ordinate(ps, "CAP", "jaccard", ~Year+Grid+Age+Sex)

# save workspace
save.image(file="/home/ahalhed/red-squirrel-w2020/R-env/rs-dbRDA.RData")