#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=server3-jan16
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------


# Generation of Tree for Phylogenetic diversity analysis
# threads to 0 to use all available cores
qiime alignment mafft \
  --p-parttree \
  --p-n-threads 0 \
  --i-sequences OTU-rep-seqs-dn-99.qza \
  --o-alignment aligned_sequences.qza
# the alignment ran in the available time, but non of the below
# started a new job (rs-server3-2.sh) with the same run parameters for the rest

qiime alignment mask \
  --i-alignment aligned_sequences.qza \
  --o-masked-alignment masked_sequences.qza

qiime phylogeny fasttree \
  --i-alignment masked_sequences.qza \
  --o-tree ./trees/unrooted_tree.qza 

qiime phylogeny midpoint-root \
  --i-tree ./trees/unrooted_tree.qza \
  --o-rooted-tree ./trees/rooted_tree.qza 

# time ran out