#!/bin/bash

#SBATCH -J DensityNonB
#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH -p gg
#SBATCH --output=density_IR_log/density_%j_%x_%a.out
#SBATCH --error=density_IR_log/density_%j_%x_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nc29578@my.utexas.edu

BUCKET=${SLURM_ARRAY_TASK_ID}
if [[ -z ${BUCKET} ]];
then
  BUCKET=$3
  if [[ -z $BUCKET ]];
then
	BUCKET=0
fi
fi

MODE="density"
POLARITY_MODE="GC"
PARTITION_COL="spacer_length"

export PYTHONUNBUFFERED=1
SCHEDULE=$1
DESIGN=${2:-"design.csv"}
COMPARTMENT=${4:-"gene"}
OUTDIR=${5:-"enrichment_out_IR"}
WINDOW_SIZE=${6:-500}

mkdir -p $OUTDIR
echo "Processing bucket ${BUCKET}..."

python density_utils.py ${SCHEDULE} \
	--out ${OUTDIR} \
        --bucket_id ${BUCKET} \
        --mode ${MODE} \
        --design ${DESIGN} \
        --compartment ${COMPARTMENT} \
        --polarity_mode ${POLARITY_MODE}  \
        --window_size ${WINDOW_SIZE} \
	--partition_col ${PARTITION_COL}
echo "Bucket ${BUCKET} has been processed succesfully."
