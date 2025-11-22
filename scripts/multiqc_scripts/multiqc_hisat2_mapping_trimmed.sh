#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/multiqc_hisat2_mapping_trimmed_bam/multiqc_analysis__hisat2_mapping_trimmed_bam_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/multiqc_hisat2_mapping_trimmed_bam/multiqc_analysis_hisat2_mapping_trimmed_bam_error_%j.e

#Define raw data paths, continers path and output paths and raw data pattern
input_BAM_dir="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/samtools/converted_bam"
output_multiqc_dir="/data/users/mkummer/RNA-Sequencing-Course/results/multiqc_hisat2_mapping_trimmed_bam"
SIF_PATH="/containers/apptainer/multiqc-1.19.sif"
BAM_SUMMARY_PATTERN="*_summary.txt"

#Ensure results directory exists
mkdir -p "$output_multiqc_dir"

# Get list of files and select the one for this array task for parallel processing
files=("$input_BAM_dir"/$BAM_SUMMARY_PATTERN)
file="${files[$SLURM_ARRAY_TASK_ID - 1]}"

# Run MultiQC on all FastQC .zip files
apptainer exec \
  --bind "$input_BAM_dir":"$input_BAM_dir" \
  --bind "$output_multiqc_dir":"$output_multiqc_dir" \
  "$SIF_PATH" multiqc -o "$output_multiqc_dir" "$input_BAM_dir" 