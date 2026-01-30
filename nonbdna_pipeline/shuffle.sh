#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH -N 2
#SBATCH -n 288
#SBATCH -p gg
#SBATCH -J Shuffler
#SBATCH --out=shuffle_log/shuffle_%x_%j.out
#SBATCH --err=shuffle_log/shuffle_%x_%j.err
#SBATCH -A BCS25073

schedule=${1:-schedule.json}
TOTAL_BUCKETS=$2
level=${3:-2}
## PARAMS ##
DATE=$(date +"%m-%d-%Y")
# OUTDIR=$SCRATCH/shuffled_genomes_level_${level}_${DATE}
OUTDIR="$SCRATCH/shuffled_genomes_level_2_01-23-2026"
mkdir -p "$OUTDIR"
##
LAST_BID=$((TOTAL_BUCKETS-1))

echo "Total buckets: $TOTAL_BUCKETS"
echo "Last bucket ID: $LAST_BID"
echo "Launching all buckets..."

for bucket in $(seq 0 $LAST_BID);
do
    echo "Starting bucket $bucket"
    srun --exclusive \
	    -N1 \
	    -n1 \
	    python shuffle.py --schedule $schedule --outdir $OUTDIR --bucket_id $bucket --level $level &
done

wait
echo "All buckets completed."

