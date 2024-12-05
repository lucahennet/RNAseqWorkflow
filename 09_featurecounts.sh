#!/usr/bin/env bash

#SBATCH --job-name=featurecounts
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/lhennet/RNAseqWorkflow/log/featurecounts_%J.out
#SBATCH --error=/data/users/lhennet/RNAseqWorkflow/log/featurecounts_%J.err

# setup of the environment and the directories
WORKDIR="/data/users/lhennet/RNAseqWorkflow"
LOGDIR="$WORKDIR/log"
OUTDIR="$WORKDIR/feature_counts"
REFGENOME="$WORKDIR/reference_genome"

mkdir -p $LOGDIR
mkdir -p $OUTDIR

module load Subread/2.0.3-GCC-10.3.0

# production of a table of counts based on the annotation file with the number of reads per genes in each sample.
# The parameter Q stands for quality. It specifies the minimum mapping quality score for a read to be included in 
# the counting process. Here a score of at least 10 will be considered.
featureCounts \
  -a $REFGENOME/Mus_musculus.GRCm39.113.gtf \
  -o $OUTDIR/gene_counts.txt \
  -F "GTF" \
  -t "exon" \
  -g "gene_id" \
  -p \
  -s 2 \
  -T 4 \
  -Q 10 \
  $WORKDIR/mapping/*_sorted.bam