#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-12:00:00
#SBATCH --mem-per-cpu 64G
#SBATCH --job-name=core-KL10
#SBATCH --output=./outputs/%x-%j.out

#---
#title: "PCNM for Core Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "04/24/2020"

# script starts here
#---

#set up
# initiated frpm: /home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
# replace KL10 with specific grid/year combo being run
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/scripts/core-KL10.R