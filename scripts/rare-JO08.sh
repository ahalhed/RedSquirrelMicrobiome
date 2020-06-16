#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-01:30:00
#SBATCH --mem-per-cpu 64G
#SBATCH --job-name=rare-JO08
#SBATCH --output=./outputs/%x-%j.out

#---
#title: "PCNM for Rare Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "02/24/2020"

# script starts here
#---

#set up
# initiated frpm: /home/ahalhed/red-squirrel-w2020/R-env/
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
# replace JO08 with specific grid/year combo being run
Rscript /home/ahalhed/projects/def-cottenie/ahalhed/red-squirrel/R-env/RedSquirrelSpatial/scripts/rare-JO08.R
