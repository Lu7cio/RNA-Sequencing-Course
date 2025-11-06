#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/fastqc_analysis_quality_check_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/fastqc_analysis_error_%j.e

module load FastQC/0.11.9-Java-11

# --- User-configurable defaults ---
READS_DIR="/data/courses/rnaseq_course/toxoplasma_de/reads_Lung"
RESULTS_DIR="/data/users/mkummer/RNA-Sequencing-Course/results"
SIF_PATH="/containers/apptainer/fastqc-0.12.1.sif"
FASTQ_PATTERN="*.fastq.*"

mkdir -p "$RESULTS_DIR"

for file in "$READS_DIR"/$FASTQ_PATTERN; do
  echo "Processing $file..."
  apptainer exec \
  --bind "$READS_DIR":"$READS_DIR" \
  --bind "$RESULTS_DIR":"$RESULTS_DIR" \
  "$SIF_PATH" fastqc -t 2 -o "$RESULTS_DIR" "$file"
done
