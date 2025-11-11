#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=06:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=hisat2_indexing
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/hisat2/indexing/indexing_%j.out
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/hisat2/indexing/indexing_%j.err

# Define paths
GENOME_FA="/data/users/mkummer/RNA-Sequencing-Course/reference_genome/Mus_musculus_c57bl6nj.C57BL_6NJ_v3.dna.primary_assembly.1.fa"
GTF_FILE="/data/users/mkummer/RNA-Sequencing-Course/reference_genome/Mus_musculus_c57bl6nj.C57BL_6NJ_v3.115.gtf"
RESULTS_DIR="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2/indexing"
SIF_PATH="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# Ensure results directory exists
mkdir -p "$RESULTS_DIR"

# Run HISAT2 indexing inside container
apptainer exec \
  --bind "$GENOME_FA":"$GENOME_FA" \
  --bind "$GTF_FILE":"$GTF_FILE" \
  --bind "$RESULTS_DIR":"$RESULTS_DIR" \
  "$SIF_PATH" \
  bash -c "
    cd $RESULTS_DIR &&
    hisat2_extract_splice_sites.py $GTF_FILE > splice_sites.txt &&
    hisat2_extract_exons.py $GTF_FILE > exons.txt &&
    hisat2-build \
      --ss splice_sites.txt \
      --exon exons.txt \
      $GENOME_FA \
      C57BL_6NJ_v3
  "
