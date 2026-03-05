#!/bin/bash

#SBATCH --job-name=s04_featurecounts
#SBATCH --array=1-36
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm/featurecounts/featurecounts_%a.out
#SBATCH --error=slurm/featurecounts/featurecounts_%a.err

eval $( spack load --sh subread )

sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p 240601_DOME2/RNA/250626_data/samplenames.txt)

basedir="240601_DOME2/RNA/250626_data"
indir="240601_DOME2/RNA/250626_data/bowtieJE2"
outdir="${basedir}/d03_featurecounts/JE2_${sample}"

#index is where the output from bowtie2-build is stored
index="240601_DOME2/RNA/250626_data/bowtiebuild/JE2"
bam="240601_DOME2/RNA/250626_data/bowtieJE2/${sample}"

mkdir -p ${outdir}

featureCounts -p -t gene -g gene_id -s 1 -a 240601_DOME2/RNA/250626_staphassembly_JE2/GCF_002085525.1/genomic.gtf -o ${outdir}/sample.txt ${bam}_sorted.bam
