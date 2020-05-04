#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-00:15:00
#SBATCH --mem-per-cpu 4G
#SBATCH --job-name=Core-Ordi-Sig
#SBATCH --output=/home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/%x-%j.out

#---

#set up
# initiated frpm: /home/ahalhed/red-squirrel-w2020/R-env/
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

# run R script
Rscript /home/ahalhed/red-squirrel-w2020/R-env/RedSquirrelSpatial/scripts/significant-ordisurf-core.R
