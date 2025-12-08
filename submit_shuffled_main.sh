#!/bin/bash

#SBATCH -J NonBDNA
#SBATCH --account=izg5139_cr_default
#SBATCH --partition=standard
#SBATCH --mem=35GB
#SBATCH --time=48:00:00
#SBATCH --output=debug_nonbdna_shuffled/slurm-%A_%a.out
#SBATCH --error=debug_nonbdna_shuffled/slurm-%A_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=nmc6088@psu.edu

mkdir -p debug_nonbdna_shuffled
PATTERN="STR,MR"
# BID=${SLURM_ARRAY_TASK_ID}
BID=$2
echo "Processing Bucket ${BID}."
# Create directories
mkdir -p extractions_STR_MR_shuffled
mkdir -p log_debug_nonbdna_shuffled
# Start job
python main.py --schedule $1 --bucket_id ${BID} --pattern ${PATTERN} \
		--logdir log_debug_nonbdna_shuffled \
		--outdir extractions_STR_MR_shuffled
echo "Bucket ${BID} is complete."
