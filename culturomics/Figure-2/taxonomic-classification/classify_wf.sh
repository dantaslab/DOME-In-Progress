#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=SLURM-outputs/gtdbtk-%A-%a.out
#SBATCH --err=SLURM-outputs/gtdbtk-%A-%a.err
#SBATCH --mail-type=ALL

export GTDBTK_DATA_PATH=/ref/gdlab/data/gtdbtk_db/release220/

ASSEMBLY_DIR="$1"
OUTPUT_DIR="$2"

conda run -n gtdbtk@2.4.1 gtdbtk classify_wf --genome_dir "${ASSEMBLY_DIR}" -x fasta --out_dir "${OUTPUT_DIR}" --cpus "${SLURM_CPUS_PER_TASK}" --pplacer_cpus "${SLURM_CPUS_PER_TASK}" --mash_db /ref/gdlab/data/gtdbtk_db/release220/
