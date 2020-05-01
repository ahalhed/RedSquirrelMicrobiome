#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu 4G
#SBATCH --job-name=Rare-Ordi-Sig
#SBATCH --output=/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/%x-%j.out

#---

#set up
# initiated frpm: /home/ahalhed/red-squirrel-w2020/R-env/
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
R

# run R script
# replace KL08 with specific grid/year combo being run
Rscript /home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/scripts/significant-ordisurf-rare.R
