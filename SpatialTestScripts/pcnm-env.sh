#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu 64G
#SBATCH --job-name=pcnm-env
#SBATCH --output=../output/%x-%j.out

#script starts here
#----------------------------------

#set up
cd /home/ahalhed/red-squirrel-w2020/R-env/

module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.0
R

# run R script
Rscript /home/ahalhed/red-squirrel-w2020/R-env/pcnm-env.R
