#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=server5-jan16
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# non-phylogenetic diversity analysis

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

# moved some rarefaction to rs-server5-tree.sh, because they needed a tree that hadn't run yet
# ran out of time before beta diversity, so created a new script (rs-server5-2.sh) for the rest of this

# Beta diversity group significance 
# by location
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv  \
  --m-metadata-column Location \
  --output-dir ./beta/beta-group-significance-non-phylogenetic-bray-Location

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Location \
  --output-dir ./beta/beta-group-significance-non-phylogenetic-jaccard-Location

# by grid
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Grid \
  --output-dir ./beta/beta-group-significance-non-phylogenetic-bray-Grid

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Grid \
  --output-dir ./beta/beta-group-significance-non-phylogenetic-jaccard-Grid

# by Food supplement
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column FoodSupplement \
  --output-dir ./beta/beta-group-significance-non-phylogenetic-bray-FoodSupplement

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column FoodSupplement \
  --output-dir ./beta/beta-group-significance-non-phylogenetic-jaccard-FoodSupplement

# bioenv
qiime diversity bioenv \
  --i-distance-matrix ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir ./bioenv/non-phylogenetic-bioenv-bray

qiime diversity bioenv \
  --i-distance-matrix ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir ./bioenv/non-phylogenetic-bioenv-jaccard

# Mantel test to two distance matrices
# two-sided Mantel test to identify correlation between two distance matrices.
qiime diversity mantel \
  --i-dm1 ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --i-dm2 ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --p-label1 Jaccard \
  --p-label2 Bray-Curtis \
  --o-visualization nonphylogeneticmanteltest.qzv
# done