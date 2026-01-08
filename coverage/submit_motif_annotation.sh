#!/bin/bash

#SBATCH -J MotifAnnotation
#SBATCH --time=48:00:00
#SBATCH -N 2
#SBATCH -n 288
#SBATCH -p gg
#SBATCH --output=annotate_motifs_IR/annotation_%j_%x_%a.out
#SBATCH --error=annotate_motifs_IR/annotation_%j_%x_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nc29578@my.utexas.edu

# PARAMS
export PYTHONUNBUFFERED=1
export POLARS_MAX_THREADS=2

#
SCHEDULE=$1
DESIGN=$2
BID=${3:-288}
OUT=${4:-$SCRATCH/"annotate_motifs_IR"}
mkdir -p $OUT

DATE=$(date +"%m-%d-%Y")
LAST_BID=$((BID-1))

# RUN PARAMS
echo "Date $DATE"
echo "Total buckets: $BID"
echo "Last bucket id: $LAST_BID"
echo "Outsourcing motif annotation to <-- $OUT..."
echo "Initializing process..."

# RUN PIPELINE
for BUCKET in $(seq 0 $LAST_BID);
do
	srun --exclusive \
		-N1 \
		-n1 \
		-t 48:00:00 \
			python annotate_motifs_terminator.py ${SCHEDULE} \
          		--bucket_id ${BUCKET} \
	  		    --design ${DESIGN} \
	  		    --outdir ${OUT} &
done
wait

cat ${OUT}/*.csv > ${OUT}/motif_annotation.all.csv
echo "Process has been completed succesfully."
