#!/bin/bash

#SBATCH -J NonBDNA
#SBATCH -p gg
#SBATCH -N 2
#SBATCH -n 288
#SBATCH --time=48:00:00
#SBATCH --output=nonbdna_main_shuffled/slurm-%j_%a.out
#SBATCH --error=nonbdna_main_shuffled/slurm-%j_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=nc29578@my.utexas.edu

SCHEDULE=$1
PATTERN=${2:-"IR"}
BID=${3:-140}
LAST_BID=$((BID - 1))
OUTDIR=$SCRATCH/"nonbdna_data_extractions"
level=${4:-2}

mkdir -p $OUTDIR
DATE=$(date +"%m-%d-%Y")
echo "Total buckets: $BID"
echo "Last bucket id: $LAST_BID"
echo "Pattern: $PATTERN"
echo "Shuffling level: $level."
##
echo "Initializing process..."
## 
for BUCKET in $(seq 0 $LAST_BID);
do
	srun --exclusive \
		-N1 \
		-n1 \
		-t 48:00:00 \
		python main.py \
	 	--schedule $SCHEDULE \
	       	--bucket_id $BUCKET \
	       	--pattern $PATTERN \
	       	--outdir "$OUTDIR/extractions_shuffled_level_${DATE}_${level}_${PATTERN}" &
done

wait
echo "Process finished!"
