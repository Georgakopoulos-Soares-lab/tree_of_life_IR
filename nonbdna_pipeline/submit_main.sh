#!/bin/bash

#SBATCH -J NonBDNA
#SBATCH --partition=gg
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=48:00:00
#SBATCH --output=nonbdna_main/slurm-%j-%a.out
#SBATCH --error=nonbdna_main/slurm-%j-%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=nc29578@my.utexas.edu

SCHEDULE=$1
PATTERN=${2:-"STR"}
BID=${SLURM_ARRAY_TASK_ID}
OUTDIR=$SCRATCH/"nonbdna_data_extractions"
mkdir -p $OUTDIR

if [[ -z $BID ]];
then
	echo "TASK_ID not provided. Will fetch from user."
	BID=$3

	if [[ -z $BID ]];
	then
		echo "No user provided TASK_ID. Exiting..."
		exit 1
	fi
fi

mkdir -p "extractions_${PATTERN}"
echo "Processing Bucket ${BID} for PATTERN ${PATTERN}."
python main.py --schedule $SCHEDULE \
	       --bucket_id $BID \
	       --pattern $PATTERN \
	       --outdir "$OUTDIR/extractions_${PATTERN}"
echo "Bucket ${BID} is complete."
