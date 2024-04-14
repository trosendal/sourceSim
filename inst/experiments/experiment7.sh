#!/bin/bash

#SBATCH --account=sourcesim
#SBATCH --partition=core
#SBATCH --ntasks=1

export OMP_NUM_THREADS=1
Rscript experiment7.R
