#!/bin/bash

#SBATCH -J SnakeBootstrapDomainLevel
#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -p gg
#SBATCH -A BCS25073
#SBATCH --output=bootstrap_domain_log/bootstrap_domain_%j_%x.out
#SBATCH --error=bootstrap_domain_log/bootstrap_domain_%j_%x.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nc29578@my.utexas.edu

j=${1:-1}
config=$2
unlock=${3:-0}
latency=${4:-30}
snakemake=${5:-"bootstrap_domain.smk"}

function UnlockDir() {
	if [[ $1 -eq 1 ]];
	then
		snakemake --unlock \
			  --snakefile $snakemake \
			  --configfile $config
	fi
}

UnlockDir $unlock
snakemake --snakefile $snakemake \
        --configfile $config \
        --keep-incomplete \
        --rerun-incomplete \
	--rerun-triggers mtime \
        --scheduler greedy \
        --keep-going \
        --latency-wait $latency \
	--cluster 'sbatch -N 1 -n 144 -p gg --time=48:00:00 -A BCS25073' \
        --jobs 36
	# --executor slurm \
echo "Process has been completed succesfully."
