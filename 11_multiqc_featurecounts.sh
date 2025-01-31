#!/bin/bash

#SBATCH --job-name=multiqc_featurecounts
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/lhennet/RNAseqWorkflow/log/multiqc_featurecounts_%J.out
#SBATCH --error=/data/users/lhennet/RNAseqWorkflow/log/multiqc_featurecounts_%J.err

# setup of the environment and the directories
WORKDIR="/data/users/${USER}/RNAseqWorkflow"
OUTDIR="$WORKDIR/feature_counts"
MULTIQC_OUTDIR="$WORKDIR/multiqc_feature_counts"
CONTAINER_PATH="/containers/apptainer/multiqc-1.19.sif"
LOGDIR="$WORKDIR/log"

mkdir -p $LOGDIR
mkdir -p $OUTDIR
mkdir -p $MULTIQC_OUTDIR

# multiqc analysis
apptainer exec \
 -B /data/ $CONTAINER_PATH multiqc $OUTDIR -o $MULTIQC_OUTDIR --interactive --title "MultiQC Report - featureCounts"