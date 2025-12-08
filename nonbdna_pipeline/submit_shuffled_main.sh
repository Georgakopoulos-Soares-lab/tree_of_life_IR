#!/bin/bash

#SBATCH -J NonBDNA
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=48:00:00
#SBATCH --output=nonbdna_main_shuffled/slurm-%j_%a.out
#SBATCH --error=nonbdna_main_shuffled/slurm-%j_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=nc29578@my.utexas.edu

SCHEDULE=$1
PATTERN=${2:-"IR"}
BID=${SLURM_ARRAY_TASK_ID}
OUTDIR=$SCRATCH/"nonbdna_data_extractions"
level=3

mkdir -p $OUTDIR

echo "Processing Bucket ${BID}."
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

python main.py --schedule $SCHEDULE \
	       --bucket_id $BID \
	       --pattern $PATTERN \
	       --outdir "$OUTDIR/extractions_shuffled_level_${level}_${PATTERN}"
echo "Bucket ${BID} is complete."
