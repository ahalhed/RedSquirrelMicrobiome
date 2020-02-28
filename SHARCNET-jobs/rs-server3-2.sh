#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=server3-2-jan23
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------


# Generation of Tree for Phylogenetic diversity analysis

# the alignment ran in the available time, but non of the below
# started a new job (rs-server3-2.sh) with the same run parameters but with 48 hours for the rest
# ran out ovf time again (masking didn't even run)
# based on the time issue and the discussion in https://docs.qiime2.org/2019.10/tutorials/phylogeny/#reducing-alignment-ambiguity-masking-and-reference-alignments
# chose to not mask the alignment moving forward and simply run the trees

# it is possible to multithread this but it may result in a different tree than single threading, so I am just not multithreading
# http://www.microbesonline.org/fasttree/#OpenMP
qiime phylogeny fasttree \
  --i-alignment aligned_sequences.qza \
  --o-tree ./trees/unrooted_tree_um.qza 

qiime phylogeny midpoint-root \
  --i-tree ./trees/unrooted_tree.qza \
  --o-rooted-tree ./trees/rooted_tree_um.qza 
