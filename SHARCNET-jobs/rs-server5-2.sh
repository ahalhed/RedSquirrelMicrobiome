#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=server5-2
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

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