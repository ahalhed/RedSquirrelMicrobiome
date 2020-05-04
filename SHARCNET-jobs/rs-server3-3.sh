#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu 32G
#SBATCH --job-name=raxmltree-jan25
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------


# Generation of Tree for Phylogenetic diversity analysis
# going to try the tree with RaxML
# https://docs.qiime2.org/2019.10/plugins/available/phylogeny/raxml/
qiime phylogeny raxml \
  --i-alignment aligned_sequences.qza \
  --o-tree ./trees/unrooted_tree_um_raxml.qza 

qiime phylogeny midpoint-root \
  --i-tree ./trees/unrooted_tree_raxml.qza \
  --o-rooted-tree ./trees/rooted_tree_um_raxml.qza 
