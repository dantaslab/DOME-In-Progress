#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --array=2-36
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm/fastqc/z_fastqc_%A_%a.out
#SBATCH --error=slurm/fastqc/z_fastqc_%A_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

eval $( spack load --sh /ldosddd )

basedir="240601_DOME2/RNA/250626_data"
rawin="240601_DOME2/RNA/250626_data/concatreads"
cleanin="240601_DOME2/RNA/250626_data/d01a_trimmomatic"
rawout="${basedir}/d02_read_qc/raw_reads/fastqc"
cleanout="${basedir}/d02_read_qc/clean_reads/fastqc"

mkdir -p ${rawout}
mkdir -p ${cleanout}

export JAVA_ARGS="-Xmx8000M"

sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p 240601_DOME2/RNA/250626_data/samplenames.txt)
samplereal=$(sed -n ${SLURM_ARRAY_TASK_ID}p 240601_DOME2/RNA/250626_data/samplenames.txt)

time sh -c \
"fastqc ${rawin}/${sample}_R1_catted.fastq.gz ${rawin}/${sample}_R2_catted.fastq.gz -o ${rawout} -t ${SLURM_CPUS_PER_TASK}; \
fastqc ${cleanin}/${samplereal}_FW_trim.fastq.gz ${cleanin}/${samplereal}_RV_trim.fastq.gz ${cleanin}/${samplereal}_UP_trim.fastq.gz -o ${cleanout} -t ${SLURM_CPUS_PER_TASK}"
