#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=32G
#SBATCH --time=06:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=hisat2_indexing
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/hisat2_GRCm39/indexing/indexing_%j.out
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/hisat2_GRCm39/indexing/indexing_%j.err

# Define paths
GENOME_FA="/data/users/mkummer/RNA-Sequencing-Course/reference_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF_FILE="/data/users/mkummer/RNA-Sequencing-Course/reference_genome/Mus_musculus.GRCm39.115.gtf"
RESULTS_DIR="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/indexing"
SIF_PATH="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# Ensure results directory exists
mkdir -p "$RESULTS_DIR"

# Run HISAT2 indexing inside container
apptainer exec \
  --bind "$GENOME_FA":"$GENOME_FA" \
  --bind "$GTF_FILE":"$GTF_FILE" \
  --bind "$RESULTS_DIR":"$RESULTS_DIR" \
  "$SIF_PATH" \
  hisat2-build \
  -p 2 \
  "$GENOME_FA" \
  "$RESULTS_DIR/GRCm39"
