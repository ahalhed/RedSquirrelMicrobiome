#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-03:00:00
#SBATCH --mem-per-cpu 64G
#SBATCH --job-name=coreTest
#SBATCH --output=./outputs/%x-%j.out


# cd /home/ahalhed/AliciaMSc/squirrel
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
# replace AG08 with specific grid/year combo being run
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/scripts/core.R