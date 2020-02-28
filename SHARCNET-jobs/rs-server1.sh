#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=server1
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# this script loads in the data

# load QIIME2 BEFORE submitting the job
# module load miniconda3
# conda activate qiime2-2019.10

# Note about the file "RS_seqsR.fasta"
# this is NOT an unzipped version of the "RS_seqs.zip" archive
# This file was built in R from that archive to make the sequences continuous (see below information for explanation)
# https://docs.qiime2.org/2019.10/tutorials/importing/#sequences-without-quality-information-i-e-fasta

# load the sequences to a QIIME2 artifact
qiime tools import \
  --input-path ./input/RS_seqsR.fasta \
  --output-path seqs.qza \
  --type 'SampleData[Sequences]'

# dereplicate the sequences
qiime vsearch dereplicate-sequences \
  --i-sequences seqs.qza \
  --o-dereplicated-table table.qza \
  --o-dereplicated-sequences rep-seqs.qza 

qiime metadata tabulate \
  --m-input-file table.qza \
  --o-visualization table.qzv

# removed the tabulate-seqs because of increased run time (and the visualization wasn't necessary)
# it put the job over the 1 hour (got killed), so I ended up tabulating the metadata outside the job
# 30 minutes likely would have been plenty of time

# tabulate the metadata
qiime metadata tabulate \
  --m-input-file ./input/RS_meta.tsv \
  --o-visualization tabulated-metadata.qzv

# close qiime2
# conda deactivate
# done