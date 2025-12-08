#!/bin/bash

#SBATCH -J Coverage
#SBATCH --time=48:00:00
#SBATCH -p gh
#SBATCH -N 1
#SBATCH --output=coverage_IR_log/density_%j_%x_%a.out
#SBATCH --error=coverage_IR_log/density_%j_%x_%a.err
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

echo "Processing bucket ${BUCKET}..."
POLARITY=0
POLARITY_MODE="GC"
PARTITION_COL="spacer_length"

export PYTHONUNBUFFERED=1
SCHEDULE=$1
DESIGN=$2
OUT=${4:-"coverage_out_IR"}
mkdir -p ${OUT}
python gff_utils.py ${SCHEDULE} \
          --bucket_id ${BUCKET} \
	  --design ${DESIGN} \
          --polarity ${POLARITY}  \
          --strand_mode ${POLARITY_MODE} \
	  --partition_col ${PARTITION_COL} \
	  --out ${OUT} \
	  --overload_biotype 1
echo "Bucket ${BUCKET} has seized to operate."
