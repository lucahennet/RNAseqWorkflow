#!/usr/bin/env bash

#SBATCH --job-name=multiqc
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10M
#SBATCH --time=00:20:00
#SBATCH --partition=pibu_el8

# setup of the environment and the directories
WORKDIR="/data/users/${USER}/RNAseqWorkflow"
OUTDIR="$WORKDIR/fastqc_output"
MULTIQC_OUTDIR="$WORKDIR/multiqc_data"
CONTAINER_PATH="/containers/apptainer/multiqc-1.19.sif"

mkdir -p $OUTDIR
mkdir -p $MULTIQC_OUTDIR

# multiqc analysis
apptainer exec \
 -B $OUTDIR $CONTAINER_PATH multiqc $OUTDIR -o $MULTIQC_OUTDIR