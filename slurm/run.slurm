#!/bin/bash
#SBATCH --time=150:00:00
#SBATCH --job-name="sim101"
#SBATCH --mem-per-cpu=3GB
#SBATCH --output=outfiles/sim101.out
#SBATCH --cpus-per-task=128

module load stack/2024-06  gcc/12.2.0
module load cmake/3.27.7
module load r/4.4.0

export OMP_NUM_THREADS=1

Rscript --vanilla anc_rand_simulation_new.R