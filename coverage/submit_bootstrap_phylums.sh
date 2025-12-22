#!/bin/bash

#SBATCH -J BootstrapPhylumLevelSnake 
#SBATCH --time=48:00:00
#SBATCH -N 3
#SBATCH -n 432
#SBATCH -p gg
#SBATCH -A BCS25073
#SBATCH --output=bootstrap_domain_log/bootstrap_phylum_%j_%x.out
#SBATCH --error=bootstrap_domain_log/bootstrap_phylum_%j_%x.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nc29578@my.utexas.edu

j=${1:-1}
config=$2
unlock=${3:-0}
latency=${4:-30}

function UnlockDir() {
	if [[ $1 -eq 1 ]];
	then
		snakemake --unlock \
			  --snakefile bootstrap_phylums.smk
			  --configfile $config
	fi 
}

UnlockDir $unlock
snakemake --snakefile bootstrap_phylums.smk
            --configfile $config \
	    --executor local \
            --keep-incomplete \
            --rerun-incomplete \
	    --rerun-triggers mtime \
            --scheduler greedy \
            --keep-going \
            --latency-wait $latency \
            --jobs $j
echo "Process has been completed succesfully."
