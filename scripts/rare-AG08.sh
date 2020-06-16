#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu 64G
#SBATCH --job-name=rare-AG08
#SBATCH --output=./outputs/%x-%j.out

#---
#title: "PCNM for Rare Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "04/24/2020"

# script starts here
#---

#set up
# initiated frpm: /home/ahalhed/red-squirrel-w2020/R-env/
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
# replace AG08 with specific grid/year combo being run
Rscript /home/ahalhed/projects/def-cottenie/ahalhed/red-squirrel/R-env/RedSquirrelSpatial/scripts/rare-AG08.R
