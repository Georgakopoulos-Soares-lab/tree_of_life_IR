#!/bin/bash

#SBATCH -J HolySnake
#SBATCH --time=48:00:00
#SBATCH -p gg
#SBATCH -N 2
#SBATCH -n 288

j=288
if [[ ! -n ${SSH_CONNECTION} ]];
then
  snakemake --snakefile extract_nonbdna.smk \
                --latency-wait 5 \
                --keep-going \
                --cores $j \
                --keep-incomplete
else
  snakemake --snakefile extract_nonbdna.smk \
        --executor local \
        --jobs $j \
        --keep-going \
        --keep-incomplete \
        --rerun-incomplete \
        --scheduler greedy \
        --latency-wait 45
fi
