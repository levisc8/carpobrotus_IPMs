#!/bin/bash
#$ -S /bin/bash
#$ -wd /work/$USER
#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -j y

#Specify job name
#$ -N shrinkageEstimation

# hard runtime
#$ -l h_rt=120

# memory per core (hard limit)
#$ -l h_vmem=8G

# N_Cores
#$ -pe smp 4

#needed when submitting a non-parallel job
#$ -binding linear:4

#create a single output directory per job
output_dir="/work/$USER/$JOB_NAME-$JOB_ID"
mkdir -p "$output_dir"

module load foss/2019b R/4.0.0

climate=$1 
func=$2 
demo=$3 
occ=$4 
output="$output_dir"/${JOB_NAME}_${JOB_ID}_${climate}_${func}

mkdir -p "$output"/stan_code
mkdir -p "$output"/models
mkdir -p "$output"/figures

Rscript "$HOME"/iceplant_shrinkage/shrinkage_test.R --climate="$climate" --func="$func" --demo="$demo" --occ="$occ" "$output"
