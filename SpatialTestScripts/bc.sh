#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu 128G
#SBATCH --job-name=bc-10KL
#SBATCH --output=./PCNM/%x-%j.out

#script starts here
#----------------------------------

# initiated from the R-env directory
#set up
# cd /home/ahalhed/red-squirrel-w2020/R-env/
module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.0
R

# run R script
Rscript /home/ahalhed/red-squirrel-w2020/R-env/scripts/bc.R
