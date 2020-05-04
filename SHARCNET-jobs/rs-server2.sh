#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=server2-jan15
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------


# de novo OTU clustering
# the original publication selected de nove clustering at 99%, so that is what I am doing too
qiime vsearch cluster-features-de-novo \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --p-perc-identity 0.99 \
  --p-threads 250 \
  --o-clustered-table OTU-table-dn-99.qza \
  --o-clustered-sequences OTU-rep-seqs-dn-99.qza

 #done