#!/usr/bin/env bash

#SBATCH --job-name=samtools_indexing
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/lhennet/RNAseqWorkflow/log/samtools_indexing_%J.out
#SBATCH --error=/data/users/lhennet/RNAseqWorkflow/log/samtools_indexing_%J.err
#SBATCH --array=1-16

# setup of the environment and the directories
WORKDIR="/data/users/lhennet/RNAseqWorkflow"
LOGDIR="$WORKDIR/log"
SAMPLELIST="$WORKDIR/metadata/samplelist.tsv"
OUTDIR="$WORKDIR/mapping"

mkdir -p $LOGDIR
mkdir -p $OUTDIR

# extracts the sample name and read1 and read2 file paths from samplelist.tsv for the array jobs
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# indexing the sorted bam files using samtools (creation of .bai files)
apptainer exec \
 -B /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
  samtools index $OUTDIR/${SAMPLE}_sorted.bam