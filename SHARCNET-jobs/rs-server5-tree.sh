#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu 16G
#SBATCH --dependency=afterok:27793601
#SBATCH --job-name=server-5tree
#SBATCH --output=./output/%x-%j.out

#script starts here
#----------------------------------

# I didn't realize these non-phylogenetic rarefactions needed the trees, so this part has been separated from the original rs-server5.sh script

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-non-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics chao1 \
  --m-metadata-file ./input/RS_meta.tsv  \
  --p-steps 100 \
  --o-visualization ./alpha/nonphylogenetic_alpha_rarefaction_chao1.qzv

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-non-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics dominance \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/nonphylogenetic_alpha_rarefaction_dominance.qzv

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-non-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics simpson_e \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/nonphylogenetic_alpha_rarefaction_simpson_e.qzv 

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-non-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics goods_coverage \
  --m-metadata-file ./input/RS_meta.tsv  \
  --p-steps 100 \
  --o-visualization ./alpha/nonphylogenetic_alpha_rarefaction_goods_coverage.qzv 

qiime diversity alpha-rarefaction \
  --i-table ./core-metrics-non-phylogenetic/rarefied_table.qza \
  --p-max-depth 900 \
  --i-phylogeny ./trees/rooted_tree.qza \
  --p-metrics simpson \
  --m-metadata-file ./input/RS_meta.tsv \
  --p-steps 100 \
  --o-visualization ./alpha/nonphylogenetic_alpha_rarefaction_simpson.qzv 

# rs-server3.sh needs to finish first