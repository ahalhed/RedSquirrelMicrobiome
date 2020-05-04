#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=server8-jan17
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# Assign taxonomy with GREEN GENES
# Access the  Green genes reference database (smaller, runs faster)
# wget -O "gg-13-8-99-nb-classifier.qza" "https://data.qiime2.org/2019.10/common/gg-13-8-99-nb-classifier.qza"
# assignment (took just under 7 hours to run)
# added --p-n-jobs so it will use all available cores and --p-chunk-size to help with running in parallel
qiime feature-classifier classify-sklearn \
 --i-classifier ./references/gg-13-8-99-nb-classifier.qza \
 --p-reads-per-batch 1000 \
 --p-n-jobs -1 \
 --i-reads rep-seqs.qza \
 --o-classification ./taxonomy/GG-taxonomy.qza
# Generating taxonomy visualization
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --m-metadata-file ./input/RS_meta.tsv \
  --o-visualization ./taxonomy/dn-taxa-bar-plots.qzv
# Extracting Taxonomic Clasification
# Phylum
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table ./taxonomy/GG-table-l2.qza
# Class
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 3 \
  --o-collapsed-table ./taxonomy/GG-table-l3.qza
# Order
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 4 \
  --o-collapsed-table ./taxonomy/GG-table-l4.qza
# Family
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table ./taxonomy/GG-table-l5.qza
# Genus
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ./taxonomy/GG-table-l6.qza
# Species
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy ./taxonomy/GG-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table ./taxonomy/GG-table-l7.qza
