#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/multiqc/multiqc_analysis_quality_check_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/multiqc/multiqc_analysis_error_%j.e

#Define raw data paths, continers path and output paths and raw data pattern
input_fastq_dir="/data/users/mkummer/RNA-Sequencing-Course/results/fastqc"
output_multiqc_dir="/data/users/mkummer/RNA-Sequencing-Course/results/multiqc"
SIF_PATH="/containers/apptainer/multiqc-1.19.sif"

#Ensure results directory exists
mkdir -p "$output_multiqc_dir"

# Run MultiQC on all FastQC .zip files
apptainer exec \
  --bind "$input_fastq_dir":"$input_fastq_dir" \
  --bind "$output_multiqc_dir":"$output_multiqc_dir" \
  "$SIF_PATH" multiqc -o "$output_multiqc_dir" "$input_fastq_dir"

