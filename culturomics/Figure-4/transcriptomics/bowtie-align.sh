#!/bin/bash

#SBATCH --job-name=bowtie2
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=2-36
#SBATCH --output=slurm/bowtie2/x_bowtie2_%A_%a.out
#SBATCH --error=slurm/bowtie2/y_bowtie2_%A_%a.out

eval $( spack load --sh bowtie2 )
#load samtools/1.14
eval $( spack load --sh /6p5wlkk )

basedir="240601_DOME2/RNA"
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p 240601_DOME2/RNA/250626_data/samplenames.txt)

#index is where the output from bowtie2-build is stored
index="240601_DOME2/RNA/250626_data/bowtiebuild/JE2"
outdir="240601_DOME2/RNA/250626_data/bowtieJE2"

#make output directory
mkdir -p ${outdir}

time bowtie2 -p ${SLURM_CPUS_PER_TASK} \
    -x ${index} \
    -S ${outdir}/${sample}_mappedreads.sam \
    -1 240601_DOME2/RNA/250626_data/d01a_trimmomatic/${sample}_FW_trim.fastq.gz \
    -2 240601_DOME2/RNA/250626_data/d01a_trimmomatic/${sample}_RV_trim.fastq.gz

# convert sam to bam, then generate sorted bam file and index it.
samtools view -bS ${outdir}/${sample}_mappedreads.sam > ${outdir}/${sample}_raw.bam
samtools sort ${outdir}/${sample}_raw.bam > ${outdir}/${sample}_sorted.bam
samtools index ${outdir}/${sample}_sorted.bam
# generate stats
samtools depth -aa ${outdir}/${sample}_sorted.bam > ${outdir}/${sample}_depth.txt
samtools stats ${outdir}/${sample}_sorted.bam > ${outdir}/${sample}_stats.txt

rm ${outdir}/${sample}_mappedreads.sam
rm ${outdir}/${sample}_raw.bam
