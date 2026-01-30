#!/bin/bash

#SBATCH -A BCS25073
#SBATCH -J NonBDNA
#SBATCH -p gg
#SBATCH -N 2
#SBATCH -n 288
#SBATCH --time=48:00:00
#SBATCH --output=nonbdna_pipeline/tss_tes_log/slurm-%j_%a.out
#SBATCH --error=nonbdna_pipeline/tss_tes_log/slurm-%j_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=nc29578@my.utexas.edu

export POLARS_MAX_THREADS=1
snakemake -s tss_tes_pipeline.smk --configfile nonbdna_pipeline/config_IR.yaml --cores 288
