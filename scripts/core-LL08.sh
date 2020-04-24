#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-00:30:00
#SBATCH --mem-per-cpu 64G
#SBATCH --job-name=core-LL08
#SBATCH --output=/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/%x-%j.out

#---
#title: "PCNM for Core Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "02/24/2020"

# script starts here
#---

#set up
# initiated frpm: /home/ahalhed/red-squirrel-w2020/R-env/
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
# replace LL08 with specific grid/year combo being run
Rscript /home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/scripts/core-LL08.R