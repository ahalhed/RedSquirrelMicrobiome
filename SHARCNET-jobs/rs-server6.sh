#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=phylogenetic-server6
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# phylogenetic diversity analysis

# PHYLOGENETIC alpha diversity
# see note in rs-server5.sh script about the filtered table
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
# ran out of time after phylogenetic_shannon_sig.qzv, restarted with more time from here forward (5 hours for remaining)
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
  --o-visualization ./alpha/core-metrics-non-phylogenetic-alpha-rarefaction.qzv

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