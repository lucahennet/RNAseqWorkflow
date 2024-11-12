#!/bin/bash

#SBATCH --job-name=get_sample_list
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1M
#SBATCH --time=00:02:00
#SBATCH --partition=pibu_el8

# setup of the environment and the directories
mkdir -p metadata

WORKDIR="/data/users/${USER}/RNAseqWorkflow"
FASTQ_FOLDER="/data/courses/rnaseq_course/toxoplasma_de/reads"
OUTPUT_FILE="${WORKDIR}/metadata/samplelist.tsv"

# creates a table with 3 columns: sample name, read 1 and read 2 to a metadata directory
for FILE in $FASTQ_FOLDER/*_*1.fastq.gz
do 
    PREFIX="${FILE%_*.fastq.gz}"
    SAMPLE=`basename $PREFIX`
    echo -e "${SAMPLE}\t$FILE\t${FILE%?.fastq.gz}2.fastq.gz" >> "$OUTPUT_FILE"
done
