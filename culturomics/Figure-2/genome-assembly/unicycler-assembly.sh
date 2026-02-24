#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --array=1-631%50
#SBATCH --output=SLURM-outputs/unicycler-%A-%a.out
#SBATCH --err=SLURM-outputs/unicycler-%A-%a.err
#SBATCH --mail-type=ALL

eval "$( spack load --sh py-unicycler@0.5.1 )"
CONFIG="config.tsv"
SAMPLE=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $2}' $CONFIG)
READ_ONE=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $3}' $CONFIG)
READ_TWO=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $4}' $CONFIG)
# Perform genome assembly (short read) for each of the isolates
unicycler -1 "${READ_ONE}" -2 "${READ_TWO}" -t "${SLURM_CPUS_PER_TASK}"  -o genome-assembly/"${SAMPLE}"