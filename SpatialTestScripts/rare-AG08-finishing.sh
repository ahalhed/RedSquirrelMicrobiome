#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-12:00:00
#SBATCH --mem-per-cpu 128G
#SBATCH --job-name=rare-AG08-finishing
#SBATCH --output=/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/%x-%j.out

#---
#title: "Finishing the rare AG script"
#author: "Alicia Halhed"
#date: "02/24/2020"

#script starts here
#---

#set up
# cd /home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
# replace AG08 with specific grid/year combo being run
Rscript /home/ahalhed/AliciaMSc/squirrel/SpatialTestScripts/rare-AG08-interacted.R
