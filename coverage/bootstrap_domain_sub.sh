#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -J SnakemakeNikol

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

if [[ ! -n "$SSH_CONNECTION" ]];
then
  echo "Local environment detected."
  echo "Initializing bioinformatics genomic compartment bootstrap analysis. (mode ${mode}; cores ${j}) [LOCAL]. Authored by Nikol Chantzi <3."
	snakemake --snakefile $snakemake \
            --configfile $config \
            --rerun-incomplete \
	    --rerun-triggers mtime \
            --reason \
            --scheduler greedy \
            --keep-going \
            --latency-wait $latency \
            --cores $j
else
  echo "SSH Connection detected."
  echo "Initializing bioinformatics genomic compartment bootstrap analysis. (mode ${mode}; cores ${j}) [SERVER]. Authored by Nikol Chantzi <3."
	snakemake --snakefile $snakemake \
            --configfile $config \
            --rerun-incomplete \
	    --rerun-triggers mtime \
            --reason \
            --keep-going \
            --jobs $j \
            --latency-wait $latency \
            --cluster-config config/cluster_bootstrap.yaml \
            --cluster "sbatch -p {cluster.partition} \
                -t {cluster.time} \
                -N {cluster.nodes} \
		-p {cluster.partition} \
                -J {cluster.jobName} \
                -o MindiBootstrapJobDetails/{cluster.jobName}-%x-%j.out \
                -e MindiBootstrapJobDetails/{cluster.jobName}-%x-%j.err"
fi
