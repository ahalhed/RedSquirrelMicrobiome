#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=0-00:30:00
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name=Figures
#SBATCH --output=./outputs/%x-%j.out

#---
#title: "PCNM for Squirrel Microbiome (SHARCNET)"
#author: "Alicia Halhed"
#date: "06/22/2020"

#script starts here
#---
# attach required modules 
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
# run R script
Rscript /home/ahalhed/projects/def-cottenie/Microbiome/RedSquirrelMicrobiome/R-env/RedSquirrelSpatial/scripts/ManuscriptFigures.R
