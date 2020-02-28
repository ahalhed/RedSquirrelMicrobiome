# Author: Alicia Halhed
# Draft Date: November 18, 2019
# Species: North American Red Squirrel
# Sample: Fecal samples
# 16S rRNA

# script starts here
# ________________________________________

# Working in the below one drive folder
cd /Users/aliciahalhed/OneDrive\ -\ University\ of\ Guelph/Alicia\'s\ Thesis/red-squirrel-data/rs-QIIME2-AH 
# Data in /Guelph/red-squirrel-data/original
# Add '&' to the end to put the process in the background
# To move a job currently running in the forground to the background, pause it with CTRL-Z and then run bg (moves last paused job to the backgroun)
# Use disown to allow terminal to close 

# Move into directory with files
# cd ./Desktop/Guelph/red-squirrel-data
# Activate QIIME2
conda activate qiime2-2019.10

# for OTU picking
# code from: https://docs.qiime2.org/2019.10/tutorials/otu-clustering/
#import sequences for OTU picking
qiime tools import \
  --input-path ./sampled5/RS_sampling.fasta \
  --output-path ./sampled5/seqs.qza \
  --type 'SampleData[Sequences]'
# dereplicate the sequences
qiime vsearch dereplicate-sequences \
  --i-sequences ./sampled5/seqs.qza \
  --o-dereplicated-table ./sampled5/table.qza \
  --o-dereplicated-sequences ./sampled5/rep-seqs.qza 

qiime metadata tabulate \
  --m-input-file ./sampled5/table.qza \
  --o-visualization ./sampled5/table.qzv

# View sequences associated with each feature
qiime feature-table tabulate-seqs \
  --i-data ./sampled5/rep-seqs.qza \
  --o-visualization ./sampled5/rep-seqs.qzv
# de novo OTU clustering
# the original publication selected de nove clustering at 99%, so that is what I am doing too
qiime vsearch cluster-features-de-novo \
  --i-table ./sampled5/table.qza \
  --i-sequences ./sampled5/rep-seqs.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table ./sampled5/table-dn-99.qza \
  --o-clustered-sequences ./sampled5/rep-seqs-dn-99.qza
# table-dn-99.qza replaces 16S_abun5.qza (FeatureTable[Frequency] data format)
# rep-seqs-dn-99.qza replaces 16S_sequences.qza (FeatureData[Sequence] format)

# tabulate the metadata
qiime metadata tabulate \
  --m-input-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/tabulated-metadata.qzv

# to look at a visualization
qiime tools view ./sampled5/tabulated-metadata.qzv

# Generation of Tree for Phylogenetic diversity analysis
qiime alignment mafft \
  --p-parttree \
  --i-sequences ./sampled5/rep-seqs-dn-99.qza  \
  --o-alignment ./sampled5/aligned_sequencess.qza

qiime alignment mask \
  --i-alignment ./sampled5/aligned_sequencess.qza \
  --o-masked-alignment ./sampled5/masked_sequences.qza
# the above was put into the background and disowned (took ~ 40 minutes to run)
qiime phylogeny fasttree \
  --i-alignment ./sampled5/masked_sequences.qza \
  --o-tree ./sampled5/unrooted_tree.qza # ran in about an hour

qiime phylogeny midpoint-root \
  --i-tree ./sampled5/unrooted_tree.qza \
  --o-rooted-tree ./sampled5/rooted_tree.qza

# PHYLOGENETIC alpha diversity
# number is iterations and should be approximately the sample size
qiime diversity core-metrics-phylogenetic \
  --i-table ./sampled5/table-dn-99.qza \
  --i-phylogeny ./sampled5/rooted_tree.qza \
  --p-sampling-depth 909 \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --output-dir ./sampled5/core-metrics-phylogenetic5

# NON-PHYLOGENETIC alpha diversity
qiime diversity core-metrics \
  --i-table ./sampled5/table-dn-99.qza \
  --p-sampling-depth 909 \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --output-dir ./sampled5/core-metrics-non-phylogenetic5

# Alpha diversity correlation
# Phylogenetic
qiime diversity alpha-correlation \
  --i-alpha-diversity ./sampled5/core-metrics-phylogenetic5/observed_otus_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-method spearman \
  --o-visualization ./sampled5/alpha/phylogenetic_observed_otus_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./sampled5/core-metrics-phylogenetic5/evenness_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-method spearman \
  --o-visualization ./sampled5/alpha/phylogenetic_evenness_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./sampled5/core-metrics-phylogenetic5/shannon_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-method spearman \
  --o-visualization ./sampled5/alpha/phylogenetic_shannon_correlation.qzv
# Non phylogenetic
qiime diversity alpha-correlation \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/observed_otus_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-method spearman \
  --o-visualization ./sampled5/alpha/nonphylogenetic_observed_otus_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/evenness_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-method spearman \
  --o-visualization ./sampled5/alpha/nonphylogenetic_evenness_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/shannon_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-method spearman \
  --o-visualization ./sampled5/alpha/nonphylogenetic_shannon_correlation.qzv


# Alpha diversity comparisons (for categorial variables)
# Phylogenetic
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./sampled5/core-metrics-phylogenetic5/observed_otus_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/alpha/phylogenetic_observed_otus_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./sampled5/core-metrics-phylogenetic5/evenness_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/alpha/phylogenetic_evenness_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./sampled5/core-metrics-phylogenetic5/shannon_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/alpha/phylogenetic_shannon_sig.qzv
# Non phylogenetic
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/observed_otus_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/alpha/nonphylogenetic_observed_otus_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/evenness_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/alpha/nonphylogenetic_evenness_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/shannon_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/alpha/nonphylogenetic_shannon_sig.qzv
# rarefraction
# phylogenetic
qiime diversity alpha-rarefaction \
  --i-table ./sampled5/core-metrics-phylogenetic5/rarefied_table.qza \
  --i-phylogeny ./sampled5/rooted_tree.qza \
  --p-max-depth 900 \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-steps 100 \
  --o-visualization ./sampled5/alpha-rarefaction/phylogenetic_alpha_rarefaction.qzv 

qiime diversity alpha-rarefaction \
  --i-table ./sampled5/core-metrics-phylogenetic5/rarefied_table.qza \
  --i-phylogeny ./sampled5/rooted_tree.qza \
  --p-max-depth 909 \
  --p-metrics 'goods_coverage', 'doubles', 'heip_e', 'simpson', 'enspie', 'faith_pd', 'dominance', 'pielou_e', 'lladser_pe', 'observed_otus', 'shannon', 'robbins', 'gini_index', 'margalef', 'menhinick', 'berger_parker_d', 'brillouin_d', 'simpson_e', 'fisher_alpha', 'ace', 'chao1', 'mcintosh_e', 'mcintosh_d', 'singles', 'michaelis_menten_fit' \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/alpha-rarefaction/phylogenetic_alpha_rarefaction.qzv 

# non phylogenetic
qiime diversity alpha-rarefaction \
  --i-table ./sampled5/core-metrics-non-phylogenetic5/rarefied_table.qza \
  --p-max-depth 900 \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-steps 100 \
  --o-visualization ./sampled5/alpha-rarefaction/core-metrics-non-phylogenetic-alpha-rarefaction.qzv

qiime diversity alpha-rarefaction \
  --i-table ./sampled5/core-metrics-non-phylogenetic5/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./sampled5/rooted_tree.qza \
  --p-metrics chao1 \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-steps 100 \
  --o-visualization ./sampled5/alpha-rarefaction/nonphylogenetic_alpha_rarefaction_chao1.qzv

qiime diversity alpha-rarefaction \
  --i-table ./sampled5/core-metrics-non-phylogenetic5/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./sampled5/rooted_tree.qza \
  --p-metrics dominance \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-steps 100 \
  --o-visualization ./sampled5/alpha-rarefaction/nonphylogenetic_alpha_rarefaction_dominance.qzv

qiime diversity alpha-rarefaction \
  --i-table ./sampled5/core-metrics-non-phylogenetic5/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./sampled5/rooted_tree.qza \
  --p-metrics simpson_e \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-steps 100 \
  --o-visualization ./sampled5/alpha-rarefaction/nonphylogenetic_alpha_rarefaction_simpson_e.qzv 

qiime diversity alpha-rarefaction \
  --i-table ./sampled5/core-metrics-non-phylogenetic5/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./sampled5/rooted_tree.qza \
  --p-metrics goods_coverage \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-steps 100 \
  --o-visualization ./sampled5/alpha-rarefaction/nonphylogenetic_alpha_rarefaction_goods_coverage.qzv 

qiime diversity alpha-rarefaction \
  --i-table ./sampled5/core-metrics-non-phylogenetic5/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./sampled5/rooted_tree.qza \
  --p-metrics simpson \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-steps 100 \
  --o-visualization ./sampled5/alpha-rarefaction/nonphylogenetic_alpha_rarefaction_simpson.qzv 


# need alpha rarefraction for phylogenetic

# Beta diversity group significance 
# by location
qiime diversity beta-group-significance \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-column Location \
  --output-dir ./sampled5/beta/beta-group-significance-non-phylogenetic5-bray-Location

qiime diversity beta-group-significance \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/jaccard_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-column Location \
  --output-dir ./sampled5/beta/beta-group-significance-non-phylogenetic5-jaccard-Location

# by grid
qiime diversity beta-group-significance \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-column Grid \
  --output-dir ./sampled5/beta/beta-group-significance-non-phylogenetic5-bray-Grid

qiime diversity beta-group-significance \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/jaccard_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-column Grid \
  --output-dir ./sampled5/beta/beta-group-significance-non-phylogenetic5-jaccard-Grid

# by Food supplement
qiime diversity beta-group-significance \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-column FoodSupplement \
  --output-dir ./sampled5/beta/beta-group-significance-non-phylogenetic5-bray-FoodSupplement

qiime diversity beta-group-significance \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/jaccard_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-column FoodSupplement \
  --output-dir ./sampled5/beta/beta-group-significance-non-phylogenetic5-jaccard-FoodSupplement

# bioenv
qiime diversity bioenv \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --output-dir ./sampled5/bioenv/bioenv5-bray

qiime diversity bioenv \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/jaccard_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --output-dir ./sampled5/bioenv/bioenv5-jaccard

# Mantel test to two distance matrices
# two-sided Mantel test to identify correlation between two distance matrices.
qiime diversity mantel \
  --i-dm1 ./sampled5/core-metrics-non-phylogenetic5/jaccard_distance_matrix.qza \
  --i-dm2 ./sampled5/core-metrics-non-phylogenetic5/bray_curtis_distance_matrix.qza \
  --p-label1 Jaccard \
  --p-label2 Bray-Curtis \
  --o-visualization ./sampled5/manteltest.qzv

# Identify core features in table of ALL samples
qiime feature-table core-features \
  --i-table ./sampled5/table-dn-99.qza \
  --p-steps 10 \
  --o-visualization ./sampled5/core-features.qzv


# Subset the data to repeat the diversity analyses for groups of interest
# repeat the above alpha/beta/mantel stuff for the subset
qiime feature-table filter-samples \
  --i-table ./sampled5/table-dn-99.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-where "Season IN ('Early Spring', 'Late Spring', 'Summer')" \
  --o-filtered-table ./sampled5/Season_filter-table.qza

# picking the squirrels in the subset that were sampled multiple times
qiime feature-table filter-samples \
  --i-table ./sampled5/table-dn-99.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-where "SquirrelID IN ('8745', '10275', '10342', '10374', '11176')" \
  --o-filtered-table ./sampled5/Squirrel_filter-table.qza
# alpha diversity
qiime diversity core-metrics \
  --i-table ./sampled5/Squirrel_filter-table.qza \
  --p-sampling-depth 909 \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --output-dir ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest

# Alpha diversity correlation
qiime diversity alpha-correlation \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/observed_otus_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-method spearman \
  --o-visualization ./sampled5/alpha/squirreltest/nonphylogenetic_observed_otus_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/evenness_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-method spearman \
  --o-visualization ./sampled5/alpha/squirreltest/nonphylogenetic_evenness_correlation.qzv

qiime diversity alpha-correlation \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/shannon_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-method spearman \
  --o-visualization ./sampled5/alpha/squirreltest/nonphylogenetic_shannon_correlation.qzv

# Alpha diversity comparisons (for categorial variables)
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/observed_otus_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/alpha/squirreltest/nonphylogenetic_observed_otus_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/evenness_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/alpha/squirreltest/nonphylogenetic_evenness_sig.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/shannon_vector.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/alpha/squirreltest/nonphylogenetic_shannon_sig.qzv


# Beta diversity group significance 
# by location
qiime diversity beta-group-significance \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-column Location \
  --output-dir ./sampled5/beta/squirreltest/beta-group-significance-non-phylogenetic5-bray-Location

qiime diversity beta-group-significance \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/jaccard_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-column Location \
  --output-dir ./sampled5/beta/squirreltest/beta-group-significance-non-phylogenetic5-jaccard-Location

# bioenv
qiime diversity bioenv \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --output-dir ./sampled5/bioenv/squirreltest/bioenv5-bray

qiime diversity bioenv \
  --i-distance-matrix ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/jaccard_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --output-dir ./sampled5/bioenv/squirreltest/bioenv5-jaccard


# Mantel test to two distance matrices
# two-sided Mantel test to identify correlation between two distance matrices.
qiime diversity mantel \
  --i-dm1 ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/jaccard_distance_matrix.qza \
  --i-dm2 ./sampled5/core-metrics-non-phylogenetic5/core-metrics-non-phylogenetic-squirreltest/bray_curtis_distance_matrix.qza \
  --p-label1 Jaccard \
  --p-label2 Bray-Curtis \
  --o-visualization ./sampled5/squirrel-manteltest.qzv


# Assign taxonomy
# Access the  Green genes reference database (smaller, runs faster)
wget -O "gg-13-8-99-nb-classifier.qza" "https://data.qiime2.org/2019.10/common/gg-13-8-99-nb-classifier.qza"
# assignment
qiime feature-classifier classify-sklearn \
 --i-classifier ./gg-13-8-99-nb-classifier.qza \
 --i-reads ./sampled5/rep-seqs.qza \
 --o-classification ./sampled5/GG-taxonomy.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/GG-taxonomy.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/dn-taxa-bar-plots.qzv
# Extracting Taxonomic Clasification
# Phylum
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/GG-taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table ./sampled5/GG-table-l2.qza
# Class
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/GG-taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table ./sampled5/GG-table-l3.qza
# Order
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/GG-taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table ./sampled5/GG-table-l4.qza
# Family
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/GG-taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table ./sampled5/GG-table-l5.qza
# Genus
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/GG-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ./sampled5/GG-table-l6.qza
# Species
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/GG-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table ./sampled5/GG-table-l7.qza

# Obtaining SILVA reference daatbase (much larger database, will likely do a better job at classifying)
wget -O "silva-132-99-nb-classifier.qza" "https://data.qiime2.org/2019.10/common/silva-132-99-nb-classifier.qza"

# Classifying taxonomies
# paused -job 17183
# got killed 9 error on this one
qiime feature-classifier classify-sklearn \
  --i-classifier ./silva-132-99-nb-classifier.qza \
  --i-reads ./sampled5/rep-seqs.qza \
  --o-classification ./sampled5/SILVA-taxonomy.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/GG-taxonomy.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --o-visualization ./sampled5/SILVA-dn-taxa-bar-plots.qzv
# Extracting Taxonomic Clasification
# Phylum
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/SILVA-taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table ./sampled5/SILVA-table-l2.qza
# Class
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/SILVA-taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table ./sampled5/SILVA-table-l3.qza
# Order
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/SILVA-taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table ./sampled5/SILVA-table-l4.qza
# Family
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/SILVA-taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table ./sampled5/SILVA-table-l5.qza
# Genus
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/SILVA-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ./sampled5/SILVA-table-l6.qza
# Species
qiime taxa collapse \
  --i-table ./sampled5/table.qza \
  --i-taxonomy ./sampled5/SILVA-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table ./sampled5/SILVA-table-l7.qza

# Close QIIME2
conda deactivate
