#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p gh
#SBATCH -J Shuffler
#SBATCH --out=shuffle_log/shuffle_%x_%j.out
#SBATCH --err=shuffle_log/shuffle_%x_%j.err

schedule=${1:-schedule.json}
level=${2:-2}
OUTDIR=$SCRATCH/shuffled_genomes/level_${level}
mkdir -p $OUTDIR

bucket=${SLURM_ARRAY_TASK_ID}
echo "Processing bucket ${bucket}..."
python shuffle.py --schedule $schedule --outdir $OUTDIR --bucket $bucket --level $level
echo "Bucket ${bucket} has been processed succesfully."
