# R script starts here
load("/home/ahalhed/red-squirrel-w2020/R-env/rs-dbRDA.RData")
library(phyloseq)
library(vegan)
library(tidyverse)

# ordinate for jaccard
rs_ordi_j <- ordinate(ps, "CAP", "jaccard", ~Year+Grid+Age+Sex)
# plot
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/plots/dbRDA-j.pdf")
plot_ordination(ps, rs_ordi_j, "samples", color="Grid",
                shape = "Sex") +
  stat_ellipse(aes(group=Grid)) +
  geom_line(aes(group=Squirrel.ID)) +
  theme_bw() +
  facet_grid(Season~Year)
dev.off()

# save workspace
save.image(file="/home/ahalhed/red-squirrel-w2020/R-env/rs-dbRDA.RData")