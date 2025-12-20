#!/bin/bash

#SBATCH -J CoverageNonB
#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH -n 144
#SBATCH -p gg
#SBATCH --output=coverage_IR_log/density_%j_%x_%a.out
#SBATCH --error=coverage_IR_log/density_%j_%x_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nc29578@my.utexas.edu

# PARAMS
export PYTHONUNBUFFERED=1
export POLARS_MAX_THREADS=2

# 
POLARITY=0
POLARITY_MODE="GC"
PARTITION_COL="spacer_length"

# INPUT
SCHEDULE=$1
DESIGN=$2
BUCKET=${3:-288}
OUT=${4:-"coverage_out_IR"}

mkdir -p $OUT
DATE=$(date +"%m-%d-%Y")
# RUN PARAMS
echo "Date $DATE"
echo "Total buckets: $BID"
echo "Initializing process..."

# RUN PIPELINE
python gff_utils.py ${SCHEDULE} \
          		--bucket_id ${BUCKET} \
	  		--design ${DESIGN} \
          		--polarity ${POLARITY}  \
          		--strand_mode ${POLARITY_MODE} \
	  		--partition_col ${PARTITION_COL} \
	  		--out ${OUT} \
	  		--overload_biotype 1
echo "Process has been completed succesfully."
