#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-582
#SBATCH --output=SLURM-outputs/dbcan-cazyme-annotation-%A-%a.out
#SBATCH --err=SLURM-outputs/dbcan-cazyme-annotation-%A-%a.err
#SBATCH --mail-type=ALL

CONFIG="$1"
SAMPLE=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $2}' $CONFIG)
FILEPATH=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $3}' $CONFIG)

conda run -n dbcan run_dbcan CAZyme_annotation --input_raw_data "${FILEPATH}" --output_dir "${SAMPLE}" --db_dir /ref/gdlab/data/dbcan_db/2025/ --mode prok 2> CAZyme-annotation.log