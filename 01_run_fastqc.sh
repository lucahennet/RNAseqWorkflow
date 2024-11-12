#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc
#SBATCH --array=1-16

# setup of the environment and the directories
WORKDIR="/data/users/lhennet/RNAseqWorkflow"
OUTDIR="$WORKDIR/fastqc_output"
SAMPLELIST="$WORKDIR/metadata/samplelist.tsv"
CONTAINER_PATH="/containers/apptainer/fastqc-0.12.1.sif"

# extracts the sample name and read1 and read2 file paths from samplelist.tsv for the array jobs
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# creates a txt file for each sample in the fastqc_output directory
OUTFILE="$OUTDIR/${SAMPLE}.txt"

mkdir -p $OUTDIR

# fastqc analysis from the container using the --bind flag to give the container the access to the sample files
apptainer exec \
 --bind /data/courses/rnaseq_course/toxoplasma_de/reads:/data/courses/rnaseq_course/toxoplasma_de/reads \
 $CONTAINER_PATH fastqc -t 2 -o $OUTDIR $READ1 $READ2 >> $OUTFILE 2>&1
