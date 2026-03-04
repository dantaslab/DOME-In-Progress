#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-631
#SBATCH --output=SLURM-outputs/checkm-%A-%a.out
#SBATCH --err=SLURM-outputs/checkm-%A-%a.err
#SBATCH --mail-type=ALL

export CHECKM_DATA_PATH=/ref/gdlab/data/checkm_db/20150116

eval "$( spack load --sh py-checkm-genome@1.2.3 )"

CONFIG="config.tsv"
SAMPLE=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $2}' $CONFIG)
FILEPATH=$(awk -v ARRAYTASKID=$SLURM_ARRAY_TASK_ID '$1==ARRAYTASKID {print $3}' $CONFIG)
# Perform checkm genome assembly quality control for each of the isolates
checkm lineage_wf -f "genome-assembly-qc/checkm/checkm-results/${SAMPLE}/${SAMPLE}.tsv" -t "${SLURM_CPUS_PER_TASK}" -x fasta --tab_table "${FILEPATH}" "/scratch/gdlab/caelanjmiller/DOME/isolate-genomics/third-iteration/genome-assembly-qc/checkm/checkm-results/${SAMPLE}" --reduced_tree