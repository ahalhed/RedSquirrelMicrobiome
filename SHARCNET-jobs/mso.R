# R script starts here
load("rs-PCNM.RData")
library(vegan)
library(phyloseq)
library(tidyverse)

# multiscale ordination
# mso_1 was run in an interactive session
#mso_1 <- mso(PCNM_CCA(XY_all[[1]], ps, numeric_all[[1]]), XY_all[[1]])
mso_2 <- mso(PCNM_CCA(XY_all[[2]], ps, numeric_all[[2]]), XY_all[[2]])
mso_3 <- mso(PCNM_CCA(XY_all[[3]], ps, numeric_all[[3]]), XY_all[[3]])
mso_4 <- mso(PCNM_CCA(XY_all[[4]], ps, numeric_all[[4]]), XY_all[[4]])
mso_5 <- mso(PCNM_CCA(XY_all[[5]], ps, numeric_all[[5]]), XY_all[[5]])
mso_6 <- mso(PCNM_CCA(XY_all[[6]], ps, numeric_all[[6]]), XY_all[[6]])
mso_7 <- mso(PCNM_CCA(XY_all[[7]], ps, numeric_all[[7]]), XY_all[[7]])
mso_8 <- mso(PCNM_CCA(XY_all[[8]], ps, numeric_all[[8]]), XY_all[[8]])

# save workspace
save.image(file="rs-PCNM.RData")