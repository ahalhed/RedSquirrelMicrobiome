#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=sampled-server7
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# having issues with this script.... 
# diversity analysis for subset of data (based on squirrels resampled)
# Identify core features in table of ALL samples
qiime feature-table core-features \
  --i-table ./OTU-table-dn-99.qza \
  --p-steps 10 \
  --o-visualization core-features.qzv
qiime feature-table core-features \
  --i-table ./filtered-table.qza \
  --p-steps 10 \
  --o-visualization filtered-core-features.qzv
# picking the squirrels in the subset that were sampled multiple times
# see R code for where the ID's came from
qiime feature-table filter-samples \
  --i-table ./filtered-table.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-where "'SquirrelID' IN ('10050', '10061', '10082', '10086', '10112', '10119', '10121', '10125', '10126', '10130', '10131', '10133', '10137', '10138', '10140', '10145', '10153', '10161', '10165', '10168', '10176', '10177', '10181', '10182', '10183', '10184', '10200', '10212', '10213', '10228', '10239', '10260', '10262', '10264', '10265', '10273', '10274', '10275', '10278', '10285', '10304', '10305', '10307', '10318', '10320', '10322', '10324', '10334', '10338', '10340', '10341', '10342', '10355', '10369', '10370', '10374', '10375', '10376', '10377', '10389', '10390', '10391', '10393', '10395', '10408', '10410', '10411', '10416', '10461', '10517', '10569', '10571', '10588', '10627', '10636', '10640', '10642', '10645', '10646', '10658', '10678', '10679', '10686', '10699', '10701', '10713', '10724', '10725', '10736', '10746', '10747', '10748', '10752', '10754', '10759', '10760', '10761', '10782', '10784', '10818', '10828', '10939', '10960', '11150', '11164', '11171', '11176', '11177', '11181', '11201', '11202', '11203', '11210', '11211', '11254', '11256', '11266', '11272', '11280', '11282', '11319', '11324', '11344', '11488', '11489', '11882', '13659', '13660', '13663', '13680', '13682', '13688', '13692', '8036', '8082', '8095', '8297', '8349', '8385', '8386', '8404', '8407', '8497', '8505', '8510', '8513', '8529', '8557', '8565', '8572', '8576', '8579', '8582', '8588', '8592', '8619', '8636', '8643', '8656', '8657', '8659', '8662', '8680', '8690', '8691', '8700', '8703', '8720', '8726', '8727', '8728', '8729', '8730', '8732', '8735', '8745', '8755', '8771', '8815')" \
  --o-filtered-table ./resampled/filter-table.qza

# NON PHYLOGENETIC
# alpha diversity
# The rarefied table contains no samples or features. Verify your table is valid and that you provided a shallow enough sampling depth.
# 179 in the "IN"
qiime diversity core-metrics \
  --i-table ./resampled/filter-table.qza \
  --p-sampling-depth 150 \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir ./resampled/core-metrics-non-phylogenetic

# Alpha diversity correlation
qiime diversity alpha-correlation \
  --i-alpha-diversity ./resampled/core-metrics-non-phylogenetic/observed_otus_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./resampled/alpha/nonphylogenetic_observed_otus_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./resampled/core-metrics-non-phylogenetic/evenness_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./resampled/alpha/nonphylogenetic_evenness_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./resampled/core-metrics-non-phylogenetic/shannon_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./resampled/alpha/nonphylogenetic_shannon_correlation.qzv

# Alpha diversity comparisons (for categorial variables)
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./resampled/core-metrics-non-phylogenetic/observed_otus_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./resampled/alpha/nonphylogenetic_observed_otus_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./resampled/core-metrics-non-phylogenetic/evenness_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./resampled/alpha/nonphylogenetic_evenness_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./resampled/core-metrics-non-phylogenetic/shannon_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./resampled/alpha/nonphylogenetic_shannon_sig.qzv

# Beta diversity group significance 
# by location
qiime diversity beta-group-significance \
  --i-distance-matrix ./resampled/core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Location \
  --output-dir ./resampled/beta/beta-group-significance-non-phylogenetic-bray-Location

qiime diversity beta-group-significance \
  --i-distance-matrix ./resampled/core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv  \
  --m-metadata-column Location \
  --output-dir ./resampled/beta/beta-group-significance-non-phylogenetic-jaccard-Location

# bioenv
qiime diversity bioenv \
  --i-distance-matrix ./resampled/core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir ./resampled/bioenv/bioenv-bray

qiime diversity bioenv \
  --i-distance-matrix ./resampled/core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir ./resampled/bioenv/bioenv-jaccard

# Mantel test to two distance matrices
# two-sided Mantel test to identify correlation between two distance matrices.
qiime diversity mantel \
  --i-dm1 ./resampled/core-metrics-non-phylogenetic/jaccard_distance_matrix.qza \
  --i-dm2 ./resampled/core-metrics-non-phylogenetic/bray_curtis_distance_matrix.qza \
  --p-label1 Jaccard \
  --p-label2 Bray-Curtis \
  --o-visualization ./resampled/nonphylogeneticmanteltest.qzv

# PHYLOGENETIC diversity analysis
# alpha diversity
qiime diversity core-metrics \
  --i-table ./resampled/filter-table.qza \
  --p-sampling-depth 179 \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir ./resampled/core-metrics-phylogenetic

# Alpha diversity correlation
qiime diversity alpha-correlation \
  --i-alpha-diversity ./resampled/core-metrics-phylogenetic/observed_otus_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./resampled/alpha/phylogenetic_observed_otus_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./resampled/core-metrics-phylogenetic/evenness_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./resampled/alpha/phylogenetic_evenness_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./resampled/core-metrics-phylogenetic/shannon_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-method spearman \
  --o-visualization ./resampled/alpha/phylogenetic_shannon_correlation.qzv

# Alpha diversity comparisons (for categorial variables)
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./resampled/core-metrics-phylogenetic/observed_otus_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./resampled/alpha/phylogenetic_observed_otus_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./resampled/core-metrics-phylogenetic/evenness_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./resampled/alpha/phylogenetic_evenness_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./resampled/core-metrics-phylogenetic/shannon_vector.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./resampled/alpha/phylogenetic_shannon_sig.qzv

# Beta diversity group significance 
# by location
qiime diversity beta-group-significance \
  --i-distance-matrix ./resampled/core-metrics-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --m-metadata-column Location \
  --output-dir ./resampled/beta/beta-group-significance-phylogenetic-bray-Location

qiime diversity beta-group-significance \
  --i-distance-matrix ./resampled/core-metrics-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv  \
  --m-metadata-column Location \
  --output-dir ./resampled/beta/beta-group-significance-phylogenetic-jaccard-Location

# bioenv
qiime diversity bioenv \
  --i-distance-matrix ./resampled/core-metrics-phylogenetic/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir ./resampled/bioenv/bioenv-bray

qiime diversity bioenv \
  --i-distance-matrix ./resampled/core-metrics-phylogenetic/jaccard_distance_matrix.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --output-dir ./resampled/bioenv/bioenv-jaccard

# Mantel test to two distance matrices
# two-sided Mantel test to identify correlation between two distance matrices.
qiime diversity mantel \
  --i-dm1 ./resampled/core-metrics-phylogenetic/jaccard_distance_matrix.qza \
  --i-dm2 ./resampled/core-metrics-phylogenetic/bray_curtis_distance_matrix.qza \
  --p-label1 Jaccard \
  --p-label2 Bray-Curtis \
  --o-visualization ./resampled/phylogeneticmanteltest.qzv

# rs-server3.sh needs to finish first