#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=server4-jan16
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# the below site was helpful
# https://library.qiime2.org/plugins/q2-fragment-insertion/16/

# load reference tree
#wget \
#  -O "sepp-refs-gg-13-8.qza" \
#  "https://data.qiime2.org/2019.10/common/sepp-refs-gg-13-8.qza"

# create tree
qiime fragment-insertion sepp \
  --i-representative-sequences rep-seqs.qza \
  --i-reference-database ./references/sepp-refs-gg-13-8.qza \
  --o-tree ./trees/insertion-tree.qza \
  --o-placements ./trees/insertion-placements.qza
# access the filtered sequences
qiime fragment-insertion filter-features \
  --i-table table.qza \
  --i-tree ./trees/insertion-tree.qza \
  --o-filtered-table filtered_table.qza \
  --o-removed-table removed_table.qza

# non zero exit status (need to debug this)
# Command '['run-sepp.sh', '/tmp/qiime2-archive-mxvd8nvo/1ead1605-6598-477d-8c2b-a571f76ffebd/data/dna-sequences.fasta', 'q2-fragment-insertion', '-x', '1', '-A', '1000', '-P', '5000', '-a', '/tmp/qiime2-archive-4mw8sejl/a14c6180-506b-4ecb-bacb-9cb30bc3044b/data/aligned-dna-sequences.fasta', '-t', '/tmp/qiime2-archive-4mw8sejl/a14c6180-506b-4ecb-bacb-9cb30bc3044b/data/tree.nwk', '-r', '/tmp/qiime2-archive-4mw8sejl/a14c6180-506b-4ecb-bacb-9cb30bc3044b/data/raxml-info.txt']' returned non-zero exit status 1.
