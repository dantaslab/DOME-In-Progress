#!/bin/bash

#SBATCH --job-name=multiqc
#SBATCH --mem=20G
#SBATCH --output=slurm/multiQC/z_multiqc_%A_%a.out
#SBATCH --error=slurm/multiQC/z_multiqc_%A_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

eval $( spack load --sh py-multiqc)

basedir="240601_DOME2/RNA/250626_data"
rawin="240601_DOME2/RNA/250626_data/concatreads"
indir="${basedir}/d02_read_qc"
cleanin="240601_DOME2/RNA/250626_data/d01a_trimmomatic"
rawout="${indir}/raw_reads/multiqc"
cleanout="${indir}/clean_reads/multiqc"

time sh -c \
"multiqc ${rawin} -o ${rawout} -n raw_multiqc; \
multiqc ${cleanin} -o ${cleanout} -n clean_multiqc"
