#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH -p gh
#SBATCH -J AGAT
#SBATCH --output=agat_out/%j_%x.out
#SBATCH --error=agat_err/%j_%x.err


# # # # # # #
OUTDIR=$SCRATCH/gff_AGAT_files
mkdir -p $OUTDIR
# # # # # # # 

BID=${SLURM_ARRAY_TASK_ID}
if [[ -z $BID ]];
then
	BID=$2
fi
echo "Working on bucket $BID."
# micromamba activate agat
python agatify_gff.py $1 --bucket_id $BID --destination $OUTDIR
