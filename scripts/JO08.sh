#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu 128G
#SBATCH --job-name=JO08-PCNM
#SBATCH --output=./%x-%j.out

#---
#title: "PCNM for Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "04/13/2020"

#script starts here
#---

#set up
# cd /home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
# replace JO08 with specific grid/year combo being run
Rscript /home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/scripts/JO08-2.R
