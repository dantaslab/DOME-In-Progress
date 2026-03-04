#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --array=1-658
#SBATCH --output=SLURM-outputs/fastp-%A-%a.out
#SBATCH --err=SLURM-outputs/fastp-%A-%a.err
#SBATCH --mail-type=ALL

eval "$( spack load --sh fastp )"
CONFIG="config.tsv"
SAMPLE=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $2}' $CONFIG)
READ_ONE=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $3}' $CONFIG)
READ_TWO=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $4}' $CONFIG)
# Perform fastp quality control for each of the isolates
fastp --in1 "${READ_ONE}" --in2 "${READ_TWO}" --out1 trimmed-reads/"${SAMPLE}_R1.trimmed.fastq.gz" --out2 trimmed-reads/"${SAMPLE}_R2.trimmed.fastq.gz" --thread "${SLURM_CPUS_PER_TASK}" -j trimmed-reads/"${SAMPLE}.json" -h trimmed-reads/"${SAMPLE}.html"
