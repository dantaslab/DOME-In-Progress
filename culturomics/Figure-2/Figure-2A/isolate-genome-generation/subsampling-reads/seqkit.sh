#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=1-631
#SBATCH --output=SLURM-outputs/seqkit-%A-%a.out
#SBATCH --err=SLURM-outputs/seqkit-%A-%a.err
#SBATCH --mail-type=ALL

eval "$( spack load --sh seqkit )"
CONFIG="$1"
INPUT="$2"
OUTPUT="$3"
SAMPLE=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $2}' "$CONFIG")
NUM_READS_50X=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $3}' "$CONFIG")

# Perform read subsampling for each of the isolates
zcat "${INPUT}/${SAMPLE}_R1.trimmed.fastq.gz" | seqkit sample -s 150 -n "${NUM_READS_50X}" -o "${OUTPUT}/${SAMPLE}_R1.trimmed.subsampled.fastq.gz"
zcat "${INPUT}/${SAMPLE}_R2.trimmed.fastq.gz" | seqkit sample -s 150 -n "${NUM_READS_50X}" -o "${OUTPUT}/${SAMPLE}_R2.trimmed.subsampled.fastq.gz"%         