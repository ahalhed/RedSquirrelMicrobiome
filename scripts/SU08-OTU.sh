#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-02:00:00
#SBATCH --mem-per-cpu 300G
#SBATCH --job-name=SU08-OTU
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
# replace SU08 with specific grid/year combo being run
Rscript /home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/scripts/SU08-OTU.R