#!/bin/bash

#SBATCH --job-name=Occupancy
#SBATCH --array=1-100
#SBATCH --time=00:30:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=500 # megabytes 
#SBATCH --partition=batch,long
#
#SBATCH --mail-user=frederik.delaender@unamur.be
#SBATCH --mail-type=ALL
#
#SBATCH --comment=Occupancy

ml purge

ml R

Rscript model-simulations.R $SLURM_ARRAY_TASK_ID
