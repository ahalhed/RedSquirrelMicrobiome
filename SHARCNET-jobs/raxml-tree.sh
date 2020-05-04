#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=7-00:00:00
#SBATCH --ntasks-per-node=16 
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name=raxml-feb14
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# https://docs.qiime2.org/2019.10/plugins/available/phylogeny/raxml/
qiime phylogeny raxml \
  --i-alignment aligned_sequences.qza \
  --p-n-threads 16 \
  --o-tree ./trees/unrooted_tree_um_raxml.qza 

qiime phylogeny midpoint-root \
  --i-tree ./trees/unrooted_tree_um_raxml.qza \
  --o-rooted-tree ./trees/rooted_tree_um_raxml.qza