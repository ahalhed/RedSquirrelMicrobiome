# Author: Alicia Halhed
# Species: North American Red Squirrel
# Sample: Fecal samples
# 16S rRNA
# working with qiime2-2019.10

# script starts here
# ________________________________________

# Working in the below one drive folder
cd /Users/aliciahalhed/OneDrive\ -\ University\ of\ Guelph/Alicia\'s\ Thesis/red-squirrel-data/rs-QIIME2-AH 
# Data in /Guelph/red-squirrel-data/original
# Add '&' to the end to put the process in the background
# To move a job currently running in the forground to the background, pause it with CTRL-Z and then run bg (moves last paused job to the backgroun)
# Use disown to allow terminal to close 

```{r}
library(Biostrings)
library(dplyr)
# read the fasta file into R (this is the unzipped version of the zip archive provided by the authors)
rs_fasta <- readDNAStringSet("/Users/aliciahalhed/OneDrive\ -\ University\ of\ Guelph/Alicia\'s\ Thesis/red-squirrel-data/rs-QIIME2-AH/RS_seqs.fasta")
seq_name <- names(rs_fasta)
sequence <- paste(rs_fasta)
# put the fasta file into a data frame
rs_fasta_df <- data.frame(seq_name, sequence)
# write the full data frame out to a file for QIIME2
# 'SampleData[Sequences]' format in Q2 doesn't seem to like the sequences being "chunked" but would take the continuous (what I did locally)
rs_fasta_df %>%
  mutate(seq_name = paste(">", rs_fasta_df$seq_name, sep="")) %>% 
  write.table("RS_seqsR.fasta", sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

# Activate QIIME2
conda activate qiime2-2019.10

# Convert FASTA file to a QIIME2 artifact (need to use an unzip version of the fasta file)
#import sequences for OTU picking
qiime tools import \
  --input-path ./input/RS_seqsR.fasta \
  --output-path seqs.qza \
  --type 'SampleData[Sequences]'

# dereplicate the sequences
qiime vsearch dereplicate-sequences \
  --i-sequences seqs.qza \
  --o-dereplicated-table table.qza \
  --o-dereplicated-sequences rep-seqs.qza 

qiime metadata tabulate \
  --m-input-file table.qza \
  --o-visualization table.qzv

# tabulate the metadata
qiime metadata tabulate \
  --m-input-file ./input/RS_meta.tsv \
  --o-visualization tabulated-metadata.qzv
# to look at a visualization
qiime tools view tabulated-metadata.qzv

# de novo OTU clustering
# the original publication selected de nove clustering at 99%, so that is what I am doing too
qiime vsearch cluster-features-de-novo \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table OTU-table-dn-99.qza \
  --o-clustered-sequences OTU-rep-seqs-dn-99.qza

# Generation of Tree for Phylogenetic diversity analysis

# de novo tree
# threads to 0 to use all available cores
qiime alignment mafft \
  --p-parttree \
  --p-n-threads 0 \
  --i-sequences OTU-rep-seqs-dn-99.qza \
  --o-alignment aligned_sequences.qza 
# https://docs.qiime2.org/2019.10/tutorials/phylogeny/#reducing-alignment-ambiguity-masking-and-reference-alignments
# below took ~29.5 hours to run (8 cores, 16G RAM)
qiime alignment mask \
  --i-alignment aligned_sequences.qza \
  --o-masked-alignment masked_sequences.qza

qiime phylogeny fasttree \
  --i-alignment aligned_sequences.qza \
  --o-tree ./trees/unrooted_tree.qza 

qiime phylogeny midpoint-root \
  --i-tree ./trees/unrooted_tree.qza \
  --o-rooted-tree ./trees/rooted_tree.qza 

# raxml currently running
# https://docs.qiime2.org/2019.10/plugins/available/phylogeny/raxml/
qiime phylogeny raxml \
  --i-alignment aligned_sequences.qza \
  --o-tree ./trees/unrooted_tree_um_raxml.qza 

qiime phylogeny midpoint-root \
  --i-tree ./trees/unrooted_tree_raxml.qza \
  --o-rooted-tree ./trees/rooted_tree_um_raxml.qza

# haven't run insertion tree yet (but will add it in if it gets run)

# NON-PHYLOGENETIC diversity analyses
# so apparently not all the samples have metadata... see ./output/server5-jan16-26377006.out
# I will therefore filter the table and then run the rest of the script
qiime feature-table filter-samples \
  --i-table OTU-table-dn-99.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-filtered-table filtered-table.qza

# NON-PHYLOGENETIC alpha diversity
qiime diversity core-metrics \
  --i-table filtered-table.qza \
  --p-sampling-depth 909 \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir core-metrics-non-phylogenetic

# Alpha diversity correlation
# Non phylogenetic
qiime diversity alpha-correlation \
  --i-alpha-diversity ./core-metrics-non-phylogenetic/observed_otus_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./alpha/nonphylogenetic_observed_otus_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./core-metrics-non-phylogenetic/evenness_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./alpha/nonphylogenetic_evenness_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./core-metrics-non-phylogenetic/shannon_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./alpha/nonphylogenetic_shannon_correlation.qzv


# Alpha diversity comparisons (for categorial variables)
# Non phylogenetic
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-non-phylogenetic/observed_otus_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./alpha/nonphylogenetic_observed_otus_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-non-phylogenetic/evenness_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./alpha/nonphylogenetic_evenness_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-non-phylogenetic/shannon_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./alpha/nonphylogenetic_shannon_sig.qzv

# rarefraction
# non phylogenetic
qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-non-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/core-metrics-non-phylogenetic-alpha-rarefaction.qzv

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-non-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics chao1 \
  --m-metadata-file ./input/RS_meta.tsv  \
  --p-steps 100 \
  --o-visualization ./alpha/nonphylogenetic_alpha_rarefaction_chao1.qzv

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-non-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics dominance \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/nonphylogenetic_alpha_rarefaction_dominance.qzv

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-non-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics simpson_e \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/nonphylogenetic_alpha_rarefaction_simpson_e.qzv 

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-non-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics goods_coverage \
  --m-metadata-file ./input/RS_meta.tsv  \
  --p-steps 100 \
  --o-visualization ./alpha/nonphylogenetic_alpha_rarefaction_goods_coverage.qzv 

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-non-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics simpson \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/nonphylogenetic_alpha_rarefaction_simpson.qzv 


# Beta diversity group significance 
# by location
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv  \
  --m-metadata-column Location \
  --o-visualization ./beta/beta-group-significance-non-phylogenetic-bray-Location.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Location \
  --o-visualization ./beta/beta-group-significance-non-phylogenetic-jaccard-Location.qzv

# by grid
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Grid \
  --o-visualization ./beta/beta-group-significance-non-phylogenetic-bray-Grid.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Grid \
  --o-visualization ./beta/beta-group-significance-non-phylogenetic-jaccard-Grid.qzv

# by Food supplement
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column FoodSupplement \
  --o-visualization ./beta/beta-group-significance-non-phylogenetic-bray-FoodSupplement.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column FoodSupplement \
  --o-visualization ./beta/beta-group-significance-non-phylogenetic-jaccard-FoodSupplement.qzv

# by Squirrel ID
# see notes on column typing
# https://forum.qiime2.org/t/summary-of-changes-to-metadata-in-qiime-2-2018-2-release/2884
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Squirrel.ID \
  --o-visualization ./beta/beta-group-significance-non-phylogenetic-bray-SquirrelID.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Squirrel.ID \
  --o-visualization ./beta/beta-group-significance-non-phylogenetic-jaccard-SquirrelID.qzv

# bioenv
qiime diversity bioenv \
  --i-distance-matrix ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./bioenv/non-phylogenetic-bioenv-bray.qzv

qiime diversity bioenv \
  --i-distance-matrix ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./bioenv/non-phylogenetic-bioenv-jaccard.qzv

# Mantel test to two distance matrices
# two-sided Mantel test to identify correlation between two distance matrices.
qiime diversity mantel \
  --i-dm1 ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --i-dm2 ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --p-label1 Jaccard \
  --p-label2 Bray-Curtis \
  --o-visualization nonphylogeneticmanteltest.qzv

# PHYLOGENETIC diversity analysis
# PHYLOGENETIC alpha diversity
# number is iterations and should be approximately the sample size
qiime diversity core-metrics-phylogenetic \
  --i-table filtered-table.qza \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-sampling-depth 909 \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir ./core-metrics-phylogenetic

# Alpha diversity correlation
# Phylogenetic
qiime diversity alpha-correlation \
  --i-alpha-diversity ./core-metrics-phylogenetic/observed_otus_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./alpha/phylogenetic_observed_otus_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./core-metrics-phylogenetic/evenness_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./alpha/phylogenetic_evenness_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./core-metrics-phylogenetic/shannon_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./alpha/phylogenetic_shannon_correlation.qzv


# Alpha diversity comparisons (for categorial variables)
# Phylogenetic
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-phylogenetic/observed_otus_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./alpha/phylogenetic_observed_otus_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-phylogenetic/evenness_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./alpha/phylogenetic_evenness_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-phylogenetic/shannon_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./alpha/phylogenetic_shannon_sig.qzv

# rarefraction
# phylogenetic
qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-phylogenetic/rarefied_table.qza \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-max-depth 900 \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/phylogenetic_alpha_rarefaction.qzv 

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/core-metrics-phylogenetic-alpha-rarefaction.qzv

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics chao1 \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/phylogenetic_alpha_rarefaction_chao1.qzv

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics dominance \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/phylogenetic_alpha_rarefaction_dominance.qzv

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics simpson_e \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/phylogenetic_alpha_rarefaction_simpson_e.qzv 

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics goods_coverage \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/phylogenetic_alpha_rarefaction_goods_coverage.qzv 

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics simpson \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/phylogenetic_alpha_rarefaction_simpson.qzv 


# Beta diversity group significance 
# by location
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Location \
  --output-dir ./beta/beta-group-significance-phylogenetic-bray-Location

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Location \
  --output-dir ./beta/beta-group-significance-phylogenetic-jaccard-Location

# by grid
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Grid \
  --output-dir ./beta/beta-group-significance-phylogenetic-bray-Grid

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Grid \
  --output-dir ./beta/beta-group-significance-phylogenetic-jaccard-Grid

# by Food supplement
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column FoodSupplement \
  --output-dir ./beta/beta-group-significance-phylogenetic-bray-FoodSupplement

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column FoodSupplement \
  --output-dir ./beta/beta-group-significance-phylogenetic-jaccard-FoodSupplement

# bioenv
qiime diversity bioenv \
  --i-distance-matrix ./core-metrics-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir ./bioenv/bioenv-bray

qiime diversity bioenv \
  --i-distance-matrix ./core-metrics-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir ./bioenv/bioenv-jaccard

# Mantel test to two distance matrices
# two-sided Mantel test to identify correlation between two distance matrices.
qiime diversity mantel \
  --i-dm1 ./core-metrics-phylogenetic/jaccard_distance_matrix.qza \
  --i-dm2 ./core-metrics-phylogenetic/bray_curtis_distance_matrix.qza \
  --p-label1 Jaccard \
  --p-label2 Bray-Curtis \
  --o-visualization ./phylogeneticmanteltest.qzv

# ASSIGN TAXONOMY
# Access the  Green genes reference database (smaller, runs faster)
wget -O "gg-13-8-99-nb-classifier.qza" "https://data.qiime2.org/2019.10/common/gg-13-8-99-nb-classifier.qza"
# assignment
qiime feature-classifier classify-sklearn \
 --i-classifier ./references/gg-13-8-99-nb-classifier.qza \
 --p-chunk-size 1000 \
 --p-n-jobs 8 \
 --i-reads rep-seqs.qza \
 --o-classification ./taxonomy/GG-taxonomy.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization dn-taxa-bar-plots.qzv
# Extracting Taxonomic Clasification
# Phylum
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table ./taxonomy/GG-table-l2.qza
# Class
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table ./taxonomy/GG-table-l3.qza
# Order
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table ./taxonomy/GG-table-l4.qza
# Family
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table ./taxonomy/GG-table-l5.qza
# Genus
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ./taxonomy/GG-table-l6.qza
# Species
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table ./taxonomy/GG-table-l7.qza

# Obtaining SILVA reference database (much larger database, will likely do a better job at classifying)
wget -O "silva-132-99-nb-classifier.qza" "https://data.qiime2.org/2019.10/common/silva-132-99-nb-classifier.qza"
# SILVA job currently in the queue
# Classifying taxonomies
qiime feature-classifier classify-sklearn \
  --i-classifier ./references/silva-132-99-nb-classifier.qza \
  --p-chunk-size 1000 \
  --p-n-jobs 8 \
  --i-reads rep-seqs.qza \
  --o-classification ./taxonomy/SILVA-taxonomy.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./taxonomy/SILVA-dn-taxa-bar-plots.qzv
# Extracting Taxonomic Clasification
# Phylum
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table ./taxonomy/SILVA-table-l2.qza
# Class
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table ./taxonomy/SILVA-table-l3.qza
# Order
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table ./taxonomy/SILVA-table-l4.qza
# Family
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table ./taxonomy/SILVA-table-l5.qza
# Genus
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ./taxonomy/SILVA-table-l6.qza
# Species
qiime taxa collapse \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/SILVA-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table ./taxonomy/SILVA-table-l7.qza

# Close QIIME2
conda deactivate
