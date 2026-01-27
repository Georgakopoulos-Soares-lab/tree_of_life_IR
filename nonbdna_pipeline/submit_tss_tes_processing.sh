#!/bin/bash

#SBATCH -J MINDI_tss_tes_processing
#SBATCH --output=logs/tss_tes_processing_%A_%a.out
#SBATCH --error=logs/tss_tes_processing_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH -p gg
#SBATCH -N 3
#SBATCH -n 412

SCHEDULE=$1 
TOTAL_BUCKETS=$2 
LAST_BID=$((TOTAL_BUCKETS-1))
PATTERN=${3:-"IR"}
POLARITY=${4:-0}
WINDOW_SIZE=${5:-500}
PARTITION_COL=${6:-"spacer_length"}

for BUCKET in $(seq 0 $LAST_BID);
do
    srun --exclusive -N 1 -n 1 python -m nonbdna_pipeline.tss_tes_processing $SCHEDULE \
        --bucket_id $BUCKET \
        --window_size $WINDOW_SIZE \
        --partition_col $PARTITION_COL
        --pattern $PATTERN \
        --strand_polarity $POLARITY &
done
wait