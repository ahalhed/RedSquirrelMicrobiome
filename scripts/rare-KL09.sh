#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu 64G
#SBATCH --job-name=rare-KL09
#SBATCH --output=/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/%x-%j.out

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
Rscript /home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/scripts/rare-KL09.R
