#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=insertiontree-jan23
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# create tree
qiime fragment-insertion sepp \
  --i-representative-sequences rep-seqs.qza \
  --i-reference-database ./references/sepp-refs-gg-13-8.qza \
  --o-tree ./trees/insertion-tree.qza \
  --o-placements insertion-placements.qza
# access the filtered sequences (some sequences might get dropped in the process of building the tree)
qiime fragment-insertion filter-features \
  --i-table table.qza \
  --i-tree ./trees/insertion-tree.qza \
  --o-filtered-table filtered_table.qza \
  --o-removed-table removed_table.qza
