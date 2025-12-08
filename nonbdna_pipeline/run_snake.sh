#!/bin/bash

#SBATCH -J HolySnake
#SBATCH --time=48:00:00
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 1

j=$1
if [[ -d ".snakemake" ]];
then
	echo "SNAKEMAKE <CHANNI>"
	# rm -rf .snakemake
	# snakemake --cleanup-metadata --snakefile extract_nonbdna.smk
fi

if [[ ! -n ${SSH_CONNECTION} ]];
then
  snakemake --snakefile extract_nonbdna.smk \
                --latency-wait 5 \
                --keep-going \
                --cores $j \
                --keep-incomplete
else
  snakemake --snakefile extract_nonbdna.smk \
	    --keep-incomplete \
	    --rerun-triggers mtime \
	    --keep-going \
	    --latency-wait 45 \
	    --cluster-config cluster_settings.yaml \
	    --cluster "sbatch -p {cluster.partition} \
	    		      -t {cluster.time} \
			      -N {cluster.nodes} \
		              -J {cluster.jobName} \
		              -o jobOut/{cluster.jobName}-%j.out \
			      -e jobOut/{cluster.jobName}-%j.err" -j $j
fi
