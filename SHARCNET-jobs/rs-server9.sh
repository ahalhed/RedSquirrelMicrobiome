#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --ntasks-per-node=16
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=silva-feb21
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# 24 hours was too short, so extended to 72 hours and increased CPUs

# Assign taxonomy with SILVA
# Obtaining SILVA reference daatbase (much larger database, will likely do a better job at classifying)
# wget -O "silva-132-99-nb-classifier.qza" "https://data.qiime2.org/2019.10/common/silva-132-99-nb-classifier.qza"

# Classifying taxonomies
qiime feature-classifier classify-sklearn \
  --i-classifier ./references/silva-132-99-nb-classifier.qza \
  --p-reads-per-batch 1000 \
  --p-n-jobs -1 \
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

# still running