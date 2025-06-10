#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --job-name="fig"
#SBATCH --mem-per-cpu=150GB
#SBATCH --output=outfiles/fig.out
#SBATCH --cpus-per-task=1

module load stack/2024-06  gcc/12.2.0
module load cmake/3.27.7
module load r/4.4.0

export OMP_NUM_THREADS=1

Rscript --vanilla figures.R