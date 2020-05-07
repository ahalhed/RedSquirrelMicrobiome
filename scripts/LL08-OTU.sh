#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-00:45:00
#SBATCH --mem-per-cpu 64G
#SBATCH --job-name=LL08-OTU
#SBATCH --output=./output/%x-%j.out

#---
#title: "Counting OTUs per grid/year"
#author: "Alicia Halhed"
#date: "05/07/2020"

#script starts here
#---

#set up
# cd /home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
# replace LL08 with specific grid/year combo being run
Rscript /home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/scripts/LL08-OTU.R