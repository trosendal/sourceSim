#!/bin/bash

#SBATCH --account=sourcesim
#SBATCH --time=01:00:00
#SBATCH --partition=core
#SBATCH --ntasks=1

export OMP_NUM_THREADS=1
Rscript experiment3.R
