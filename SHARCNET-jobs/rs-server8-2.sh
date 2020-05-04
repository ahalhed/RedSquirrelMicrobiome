#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=server8-2-jan19
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# modification of rs-server8.sh, needed to use the filtered table (all samples need matching metadata)
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table filtered-table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./taxonomy/dn-taxa-bar-plots.qzv
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

# going to need to rerun the SILVA with this correction too