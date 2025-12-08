#!/bin/bash

j=${1:-1}
config=$2
unlock=${3:-0}
latency=${4:-30}
snakemake=${5:-"bootstrap_phylums.smk"}

function UnlockDir() {
	if [[ $1 -eq 1 ]];
	then
		snakemake --unlock --snakefile $snakemake --configfile $config --cores $j
	fi
}
UnlockDir $3
if [[ ! -n "$SSH_CONNECTION" ]];
then
  echo "Local environment detected."
  echo "Initializing bioinformatics genomic compartment bootstrap analysis. (mode ${mode}; cores ${j}) [LOCAL]. Authored by Nikol Chantzi <3."
	snakemake --snakefile $snakemake \
            --configfile $config
            --rerun-incomplete \
	    --rerun-triggers mtime \
            --reason \
            --use-conda \
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
            --reason \
            --keep-going \
            --jobs $j \
            --latency-wait $latency \
            --cluster-config config/cluster_bootstrap.yaml \
            --cluster "sbatch -p {cluster.partition} \
                -t {cluster.time} \
                -p {cluster.partition} \
                -N {cluster.nodes} \
                -J {cluster.jobName} \
                -o MindiBootstrapJobDetails/{cluster.jobName}-%x-%j.out \
                -e MindiBootstrapJobDetails/{cluster.jobName}-%x-%j.err"
fi
