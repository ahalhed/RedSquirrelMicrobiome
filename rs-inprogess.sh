# the below can apparently be used to do the alignment and phylogeny in one step
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --output-dir mafft-fasttree-output
#need to look into it a bit more

# not working stuff (need to work on this more) - associates with rs-analysis-AH5.sh

# Generate heatmap
# denovo IDs missing from metadata (fresh OTU picking might help)
qiime feature-table heatmap \
  --i-table ./sampled5/table-dn-99.qza \
  --p-color-scheme RdYlBu \
  --o-visualization ./sampled5/5heatmap.qzv
# the above ran for like three hours then gave the following error:
#Plugin error from feature-table:
#Image size of 899726x1334 pixels is too large. It must be less than 2^16 in each direction.
#Debug info has been saved to /var/folders/wd/s68yr8892s137xwylsfl43nw0000gn/T/qiime2-q2cli-err-slq0s2xp.log
qiime feature-table heatmap \
  --i-table ./sampled5/table-dn-99.qza \
  --m-feature-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-feature-metadata-column Season \
  --p-title 'Season Bray Curtis' \
  --p-metric braycurtis \
  --p-color-scheme RdYlBu \
  --o-visualization ./sampled5/season_5heatmap.qzv
# got 'Killed: 9' error - something to do with too much memory

# making an insertion tree
# get the reference database (green genes, preformtted)
wget \
  -O "sepp-refs-gg-13-8.qza" \
  "https://data.qiime2.org/2019.10/common/sepp-refs-gg-13-8.qza"
# make insertion tree
qiime fragment-insertion sepp \
  --i-representative-sequences ./sampled5/rep-seqs-dn-99.qza \
  --i-reference-database ./sepp-refs-gg-13-8.qza \
  --o-tree ./tree_insertion.qza \
  --o-placements ./tree_placements.qza \
  --p-threads 4

# longitudinal analysis
# need to play around with the volatility plots some more
# by year
qiime longitudinal volatility \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-file ./sampled5/core-metrics-phylogenetic5/unweighted_unifrac_pcoa_results.qza \
  --p-state-column Year \
  --p-individual-id-column SquirrelID \
  --p-default-group-column Grid \
  --o-visualization ./sampled5/pc_vol_year.qzv
# by month
qiime longitudinal volatility \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-file ./sampled5/core-metrics-phylogenetic5/unweighted_unifrac_pcoa_results.qza \
  --p-state-column Month \
  --p-individual-id-column SquirrelID \
  --p-default-group-column Grid \
  --o-visualization ./sampled5/pc_vol_month.qzv
# by Date
qiime longitudinal volatility \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-file ./sampled5/core-metrics-phylogenetic5/unweighted_unifrac_pcoa_results.qza \
  --p-state-column Date \
  --p-individual-id-column SquirrelID \
  --p-default-group-column Grid \
  --o-visualization ./sampled5/pc_vol.qzv
# distance based analysis
qiime longitudinal first-distances \
  --i-distance-matrix ./sampled5/core-metrics-phylogenetic5/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --p-state-column Date \
  --p-individual-id-column SquirrelID \
  --o-first-distances ./sampled5/from_first_unifrac.qza
# visualization
qiime longitudinal volatility \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv \
  --m-metadata-file ./sampled5/from_first_unifrac.qza \
  --p-state-column Date \
  --p-individual-id-column SquirrelID \
  --p-default-group-column Grid \
  --o-visualization ./sampled5/from_first_unifrac_vol.qzv

# lme
qiime longitudinal linear-mixed-effects \
  --m-metadata-file ./sampled5/RS_metadata_sub2.tsv  \
  --m-metadata-file ./sampled5/from_first_unifrac.qza \
  --p-metric Age \
  --p-state-column Date \
  --p-individual-id-column SquirrelID \
  --p-group-columns Location \
  --o-visualization ./from_first_unifrac_lme.qzv


# set up for SILVA
# is a better, more updated reference databased than green genes

# exporting to biom/tsv
qiime tools export --input-path OTU-table-dn-99.qza --output-path ./exported
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
