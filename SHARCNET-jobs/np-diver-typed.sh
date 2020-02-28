#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=np-diversity-jan18
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# replacement for the non-phylogenetic analyses in server5 and server5-2, run after the "numeric" ID's were set to categorical in the metadata

# NON-PHYLOGENETIC diversity analyses
# NON-PHYLOGENETIC alpha diversity
qiime diversity core-metrics \
  --i-table OTU-table-dn-99.qza \
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

# by Squirrel ID
# see notes on column typing
# https://forum.qiime2.org/t/summary-of-changes-to-metadata-in-qiime-2-2018-2-release/2884
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Squirrel.ID \
  --output-dir ./beta/beta-group-significance-non-phylogenetic-bray-SquirrelID

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Squirrel.ID \
  --output-dir ./beta/beta-group-significance-non-phylogenetic-jaccard-SquirrelID

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