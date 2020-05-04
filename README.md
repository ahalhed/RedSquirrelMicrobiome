# RedSquirrelSpatial
Files for red squirrel microbiome spatial analysis.

## What can you find in this directory?
The four character alphanumeric identifier in each filename (AG08, for example) corresponds to the sampling grid and collection year for the samples of interest in that script. There are 8 possible combinations (3 years, 6 grids). Files labelled with "core" were run with a community matrix containing only the core OTUs (present in >95% of collected samples and comprising a minimum of 0.1% of the observed microbial abundance). Those labelled "rare" were run with a community matrix containing those OTUs NOT present in the core microbiome.

### Script files
The files in the scripts folder are a combination of shell and R scripts, labelled the same with different extensions for corresponding scripts. These particular script files are post-QIIME2 analysis. 

There is a single QIIME2 script in this directory, labelled *rs-q2.sh*. This particular script was run in chunks as individual batch jobs, as many steps were quite long running. These individual jobs are stored in a separate private GitHub repository. The resulting QIIME2 visualizations and artifacts are stored in a directory on the GRAHAM cluster of SHARCNET (ComputeCanada).

### Plots
In addition to STDOUT, discussed below, most scripts output a variety of figures corresponding to different steps of the analysis. The labelling of these files is consistent with the type of plot, microbial community (full (no added label), core, or rare), and grid/year combination.

### Output files
Files labelled with the ".out" extension are text files containing the STDOUT from all SHARCNET batch jobs. These files are named according to the grid/year combination whose spatial analysis STDOUT is contained within the file and for the SHARCNET job number. The files labelled with "bc" are output for partitioning Bray-Curtis dissimilarities. This code was included in the full scripts, but didn't run in many cases due to failure of variation decomposition with parsimonious variables due to insignificant environmental variables. The full scripts have since been reorganized to put the partitioning before the variation decomposition to mitigate the need to run additional jobs in future replicates of the tasks.

