#!/bin/bash
#SBATCH --account=def-jwhitton
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16000M
#SBATCH --time=0-48:00
#SBATCH --job-name=nquack_fix
#SBATCH --output=%x.out

module load r samtools

Rscript ./scripts/nQuack_cluster.R