setwd("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/")
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(vegan)
library(zCompositions)
# devtools::install_github('ggloor/CoDaSeq/CoDaSeq')
library(CoDaSeq)
library(tidyverse)

theme_set(theme_bw())

otu <- read_qza("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/filtered-table-10.qza")$data
map <- read_q2metadata("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/input/RS_meta.tsv") # this is metadata
#~/OneDrive - University of Guelph/Alicia's Thesis/red-squirrel-w2020
#otu <- otu[which(rownames(otu) %in% occ_abun[which(occ_abun$fill == "core"),]$otu),]
otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- rownames_to_column(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance data frame
# using a compositional approach
# impute the OTU table (takes OTUs as columns)
OTUimp <- cmultRepl(t(otu), label=0, method="CZM") # all OTUs
# compute the aitchison values (subsequent analyses needs OTUs as columns)
OTUclr <- t(codaSeq.clr(OTUimp))

# Ranking OTUs based on their occupancy
# For caluclating raking index we included following conditions:
#   - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
#   - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), OTUclr) %>% 
  gather(SampleID, abun, -otu) %>%
  left_join(map, by = 'SampleID') %>%
  group_by(otu, CollectionDate) %>%
  summarise(time_freq=sum(abun>0)/length(abun),            # frequency of detection between time points
            coreTime=ifelse(time_freq == 1, 1, 0)) %>%     # 1 only if occupancy 1 with specific time, 0 if not
  group_by(otu) %>%
  summarise(sumF=sum(time_freq),
            sumG=sum(coreTime),
            nS=length(CollectionDate)*2,           
            Index=(sumF+sumG)/nS)                 # calculating weighting Index based on number of time points detected and 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,
            rank=Index) %>%
  arrange(desc(rank))

# Calculating the contribution of ranked OTUs to the euclidean similarity
EUCaddition <- NULL

# calculating EUC dissimilarity based on the 1st ranked OTU
otu_start=otu_ranked$otu[1]                   
start_matrix <- as.matrix(OTUclr[otu_start,])
start_matrix <- t(start_matrix)
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) dist(start_matrix[,x[1:2]]))
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
EUCaddition <- rbind(EUCaddition,df_s)
# calculating the EUC dissimilarity of the whole dataset 
x <-  apply(combn(ncol(OTUclr), 2), 2, function(x) dist(start_matrix[,x[1:2]]))   
x_names <- apply(combn(ncol(OTUclr), 2), 2, function(x) paste(colnames(OTUclr)[x], collapse=' - '))
df_full <- data.frame(x_names,x)
names(df_full)[2] <- length(rownames(OTUclr))
EUCfull <- left_join(EUCaddition,df_full, by='x_names')

rownames(EUCfull) <- EUCfull$x_names
temp_EUC <- EUCfull
temp_EUC$x_names <- NULL
temp_EUC_matrix <- as.matrix(temp_EUC)

EUC_ranked <- data.frame(rank = as.factor(row.names(t(temp_EUC_matrix))),t(temp_EUC_matrix)) %>% 
  gather(comparison, EUC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanEUC=mean(EUC)) %>%            # mean Bray-Curtis dissimilarity
  arrange(desc(-MeanEUC)) %>%
  mutate(proportionEUC=MeanEUC/max(MeanEUC))   # proportion of the dissimilarity explained by the n number of ranked OTUs
Increase=EUC_ranked$MeanEUC[-1]/EUC_ranked$MeanEUC[-length(EUC_ranked$MeanEUC)]
increaseDF <- data.frame(IncreaseEUC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
EUC_ranked <- left_join(EUC_ranked, increaseDF)
EUC_ranked <- EUC_ranked[-nrow(EUC_ranked),]

#Creating thresholds for core inclusion 

#Method: 
#A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
fo_difference <- function(pos){
  left <- (EUC_ranked[pos, 2] - EUC_ranked[1, 2]) / pos
  right <- (EUC_ranked[nrow(EUC_ranked), 2] - EUC_ranked[pos, 2]) / (nrow(EUC_ranked) - pos)
  return(left - right)
}
EUC_ranked$fo_diffs <- sapply(1:nrow(EUC_ranked), fo_difference)

elbow <- which.max(EUC_ranked$fo_diffs)

#B) Final increase in EUC similarity of equal or greater then 2% 
lastCall <- last(as.numeric(EUC_ranked$rank[(EUC_ranked$IncreaseEUC>=1.02)]))

#Creating occupancy abundance plot
occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(EUC_ranked$rank[(EUC_ranked$IncreaseEUC>=1.02)]))]] <- 'core'
# add 95% occupancy threshold for core
occ_abun$Community <- ifelse(occ_abun$otu_occ >= 0.95 & occ_abun$fill == "core", "Confirmed Core",
                             ifelse(occ_abun$otu_occ < 0.95 & occ_abun$fill == "core", "Core Candidate",
                                    "Confirmed Rare"))
# add a taxonomy column
tax <- read_qza("/home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/taxonomy/SILVA-taxonomy-10.qza")$data %>%
  rename("otu" = "Feature.ID")
# clean up/separatee taxonomy labels
tax$Taxon <- tax$Taxon %>%
  str_replace_all("D_0__", "") %>%
  str_replace_all("D_1__", "") %>%
  str_replace_all("D_2__", "") %>%
  str_replace_all("D_3__", "") %>%
  str_replace_all("D_4__", "") %>%
  str_replace_all("D_5__", "") %>%
  str_replace_all("D_6__", "")
occ_abunT <- tax %>% 
  mutate("Kingdom" = word(.$Taxon, 1, sep = ";"), #k__
         "Phylum" = word(.$Taxon, 2, sep = ";"), #p__
         "Class" = word(.$Taxon, 3, sep = ";"), #c__
         "Order" = word(.$Taxon, 4, sep = ";"), #o__
         "Family" = word(.$Taxon, 5, sep = ";"), #f__
         "Genus" = word(.$Taxon, 6, sep = ";"), #g__
         "Species" = word(.$Taxon, 7, sep = ";")) %>% #s__
  # join to core labels
  right_join(occ_abun)

# exporting the data frame with which are core
# to load into manuscript figure file
write.table(occ_abunT, file = "./data/coreA.csv", sep = ",", quote = F, row.names = F)

# so apparently with the aitchison, this doesn't work and everything is rare