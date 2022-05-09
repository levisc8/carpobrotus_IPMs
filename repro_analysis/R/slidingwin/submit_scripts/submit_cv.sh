#!/bin/bash
#SBATCH --chdir /work/$USER
#SBATCH -o /work/%u/%x-%j.out

# Specify job name
#SBATCH -J allModelComparison

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=levisc8@gmail.com

# hard runtime
#SBATCH -t 250:00:00

# N_Cores
#SBATCH -c 8

# memory per core (hard limit)
#SBATCH --mem-per-cpu=4G


#create a single output directory per job
output_dir="/work/$USER/$SLURM_JOB_NAME-$SLURM_JOB_ID"
mkdir -p "$output_dir"

module load foss/2019b R/4.0.0

vr=$1
climate=$2

output="$output_dir"/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${climate}_${vr}
mkdir -p "$output"/cv_results

Rscript "$HOME"/iceplant_shrinkage/shrink_compare.R --vital="$vr" --climate="$climate" "$output"
