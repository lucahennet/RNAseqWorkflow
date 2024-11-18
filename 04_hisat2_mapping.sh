#!/usr/bin/env bash

#SBATCH --job-name=hisat2_mapping
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/lhennet/RNAseqWorkflow/log/mapping_%J.out
#SBATCH --error=/data/users/lhennet/RNAseqWorkflow/log/mapping_%J.err
#SBATCH --array=1-16

# setup of the environment and the directories
WORKDIR="/data/users/lhennet/RNAseqWorkflow"
LOGDIR="$WORKDIR/log"
SAMPLELIST="$WORKDIR/metadata/samplelist.tsv"
OUTDIR="$WORKDIR/mapping"
INDEXGENOME="$WORKDIR/hisat_indexing"

mkdir -p $LOGDIR
mkdir -p $OUTDIR
mkdir -p $OUTDIR/summary

# extracts the sample name and read1 and read2 file paths from samplelist.tsv for the array jobs
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# mapping of the genome with hisat2
apptainer exec \
 -B /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
 hisat2 -x $INDEXGENOME/GRCm39 -1 $READ1 -2 $READ2 -S $OUTDIR/$SAMPLE.sam -p 4 --rna-strandness RF \
 --summary-file $OUTDIR/summary/${SAMPLE}_summary.txt