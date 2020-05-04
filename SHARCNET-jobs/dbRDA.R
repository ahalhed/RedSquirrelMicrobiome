# R script starts here
load("/home/ahalhed/red-squirrel-w2020/R-env/rs-full.RData")
library(phyloseq)
library(vegan)
library(tidyverse)


rs_ordi_bc <- ordinate(ps, "CAP", "bray", ~Year+Grid+Age+Sex)
# plot
pdf(file = "/home/ahalhed/red-squirrel-w2020/R-env/plots/dbRDA-bc.pdf")
plot_ordination(ps, rs_ordi_bc, "samples", color="Grid",
                shape = "Sex") +
  stat_ellipse(aes(group=Grid)) +
  geom_line(aes(group=Squirrel.ID)) +
  theme_bw() +
  facet_grid(Season~Year)
dev.off()

# save workspace
save.image(file="/home/ahalhed/red-squirrel-w2020/R-env/rs-dbRDA.RData")