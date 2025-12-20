#!/bin/bash

#SBATCH -J DensityNonB
#SBATCH --time=48:00:00
#SBATCH -N 3
#SBATCH -n 432
#SBATCH -p gg
#SBATCH --output=density_IR_log/density_%j_%x_%a.out
#SBATCH --error=density_IR_log/density_%j_%x_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nc29578@my.utexas.edu

MODE="density"
POLARITY_MODE="GC"
PARTITION_COL="spacer_length"

export POLARS_MAX_THREADS=2
export PYTHONUNBUFFERED=1
SCHEDULE=$1
DESIGN=${2:-"design.csv"}
BID=${3:-432}
COMPARTMENT=${4:-"gene"}
OUTDIR=${5:-"enrichment_out_IR"}
WINDOW_SIZE=${6:-500}

mkdir -p $OUTDIR
DATE=$(date +"%m-%d-%Y")
LAST_BID=$((BID-1))
echo "Date $DATE"
echo "Total buckets: $BID"
echo "Last bucket id: $LAST_BID"
echo "Pattern: $PATTERN"
echo "Initializing process..."

##
for BUCKET in $(seq 0 $LAST_BID);
do
	srun --exclusive \
		-N1 \
		-n1 \
		-t 48:00:00 \
		python density_utils.py ${SCHEDULE} \
			--out ${OUTDIR} \
        		--bucket_id ${BUCKET} \
        		--mode ${MODE} \
        		--design ${DESIGN} \
        		--compartment ${COMPARTMENT} \
			--biotype 1 \
        		--polarity_mode ${POLARITY_MODE}  \
        		--window_size ${WINDOW_SIZE} \
			--partition_col ${PARTITION_COL} & 
done

wait
echo "Process has been completed succesfully."
