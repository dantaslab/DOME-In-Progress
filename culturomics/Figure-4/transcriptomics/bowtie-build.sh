#!/bin/bash

#SBATCH --job-name=bowtieBuild
#SBATCH --mem=4G
#SBATCH --output=slurm/bowtie/x_bowtieBuild_%A_%a.out
#SBATCH --error=slurm/bowtie/y_bowtieBuild_%A_%a.out

eval $( spack load --sh bowtie2 )

time bowtie2-build 240601_DOME2/RNA/250626_staphassembly_JE2/GCF_002085525.1/GCF_002085525.1_ASM208552v1_genomic.fna 240601_DOME2/RNA/250626_data/bowtiebuild/JE2
