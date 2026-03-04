#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-631
#SBATCH --output=SLURM-outputs/quast-%A-%a.out
#SBATCH --err=SLURM-outputs/quast-%A-%a.err
#SBATCH --mail-type=ALL

eval "$( spack load --sh py-quast@5.2.0 )"

CONFIG="genome-assembly-qc/quast/config.tsv"
SAMPLE=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $2}' $CONFIG)
FILEPATH=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $3}' $CONFIG)
# Perform quality control assessment of genome assemblies
quast.py "${FILEPATH}" -l "${SAMPLE}" -o genome-assembly-qc/quast/"${SAMPLE}"