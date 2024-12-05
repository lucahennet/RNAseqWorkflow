#!/usr/bin/env bash

#SBATCH --job-name=qualimap
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:20:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/lhennet/RNAseqWorkflow/log/qualimap_%J.out
#SBATCH --error=/data/users/lhennet/RNAseqWorkflow/log/qualimap_%J.err
#SBATCH --array=1

# setup of the environment and the directories
WORKDIR="/data/users/${USER}/RNAseqWorkflow"
LOGDIR="$WORKDIR/log"
OUTDIR="$WORKDIR/qualimap"
SAMPLELIST="$WORKDIR/metadata/samplelist.tsv"
REFGENOME="$WORKDIR/reference_genome"

mkdir -p $OUTDIR
mkdir -p $LOGDIR

module load Qualimap/2.3-foss-2021a-R-4.3.2

# unzipping the annotation file if it was not already done before
cd $REFGENOME

gunzip Mus_musculus.GRCm39.113.gtf.gz

if [ ! -f "Mus_musculus.GRCm39.113.gtf" ]; then
    gunzip -c Mus_musculus.GRCm39.113.gtf.gz > Mus_musculus.GRCm39.113.gtf
fi

# extracts the sample name from samplelist.tsv for the array jobs
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`

# post-alignment quality control using qualimap
unset DISPLAY

qualimap rnaseq \
    -bam $WORKDIR/mapping/${SAMPLE}.bam \
    -gtf $REFGENOME/Mus_musculus.GRCm39.113.gtf \
    -outdir $OUTDIR/$SAMPLE \
    -outfile "${SAMPLE}" \
    -outformat "HTML" \
    -p strand-specific-reverse \
    --java-mem-size=12G >& "$OUTDIR/${SAMPLE}_qualimap.log"