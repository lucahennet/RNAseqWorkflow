#!/usr/bin/env bash

#SBATCH --job-name=hisat2_indexing
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/lhennet/RNAseqWorkflow/log/indexing_%J.out
#SBATCH --error=/data/users/lhennet/RNAseqWorkflow/log/indexing_%J.err

# setup of the environment and the directories
WORKDIR="/data/users/lhennet/RNAseqWorkflow"
LOGDIR="$WORKDIR/log"
OUTDIR="$WORKDIR/hisat_indexing"
REFGENOME="$WORKDIR/reference_genome"

mkdir -p $LOGDIR
mkdir -p $OUTDIR
mkdir -p $REFGENOME

# download the reference genome and the annotations if was not already done before
cd $REFGENOME

FASTA_URL="https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
FASTA_FILE="Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

if [ ! -f "$FASTA_FILE" ] && [ ! -f "${FASTA_FILE%.gz}" ]; then wget $FASTA_URL; fi

GTF_URL="https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz"
GTF_FILE="Mus_musculus.GRCm39.113.gtf.gz"

if [ ! -f "$GTF_FILE" ]; then wget $GTF_URL; fi

# unzipping the genome if it was not already done before
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

if [ ! -f "Mus_musculus.GRCm39.dna.primary_assembly.fa" ]; then
    gunzip -c Mus_musculus.GRCm39.dna.primary_assembly.fa.gz > Mus_musculus.GRCm39.dna.primary_assembly.fa
fi

# indexing the genome with hisat2
apptainer exec \
 -B /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
 hisat2-build Mus_musculus.GRCm39.dna.primary_assembly.fa $OUTDIR/GRCm39