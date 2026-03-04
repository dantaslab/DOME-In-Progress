#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=SLURM-outputs/gtdbtk-denovo.out
#SBATCH --err=SLURM-outputs/gtdbtk-denovo.err
#SBATCH --mail-type=ALL

export GTDBTK_DATA_PATH="/ref/gdlab/data/gtdbtk_db/release220/"


ASSEMBLY_DIR="$1"
OUTPUT_DIR="$2"
BAC_FILE="$3"

conda run -n gtdbtk@2.4.1 gtdbtk de_novo_wf --genome_dir "${ASSEMBLY_DIR}" -x fasta --out_dir "${OUTPUT_DIR}" --cpus "${SLURM_CPUS_PER_TASK}" --bacteria --outgroup_taxon o__Pseudomonadales --custom_taxonomy_file "${BAC_FILE}" --force --skip_gtdb_refs
