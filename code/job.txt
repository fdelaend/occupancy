#!/bin/bash

#SBATCH --job-name=Occupancy
#SBATCH --array=1-100
#SBATCH --time=00:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=500 # megabytes 
#SBATCH --partition=batch,long
#
#SBATCH --mail-user=frederik.delaender@unamur.be
#SBATCH --mail-type=ALL
#
#SBATCH --comment=Feasibility

ml purge

ml releases/2020b

ml R/4.0.4-foss-2020b

Rscript Feasibility_computations_cluster.R $SLURM_ARRAY_TASK_ID
