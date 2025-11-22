#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/fastqc_fastp_trimmed/fastqc_fastp_trimmed_analysis_quality_check_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/fastqc_fastp_trimmed/fastqc_fastp_trimmed_analysis_error_%j.e
#SBATCH --array=1-30


#Define fastp trimmed data file paths, containers path and output paths and raw data pattern
fastp_trimming_dir="/data/users/mkummer/RNA-Sequencing-Course/results/fastp_trimming"
RESULTS_DIR="/data/users/mkummer/RNA-Sequencing-Course/results/fastqc_fastp_trimmed"
SIF_PATH="/containers/apptainer/fastqc-0.12.1.sif"
FASTQ_PATTERN="*_trimmed.fastq.*"

#Ensure results directory exists
mkdir -p "$RESULTS_DIR"

# Get list of files and select the one for this array task for parallel processing
files=("$fastp_trimming_dir"/$FASTQ_PATTERN)
file="${files[$SLURM_ARRAY_TASK_ID - 1]}"

#Loop over all fastp trimmed fastq files in the fastp trimming directory and run FastQC inside the Apptainer container
if [ -f "$file" ]; then
  echo "Processing $file..."
  apptainer exec \
  --bind "$fastp_trimming_dir":"$fastp_trimming_dir" \
  --bind "$RESULTS_DIR":"$RESULTS_DIR" \
  "$SIF_PATH" fastqc -t 1 -o "$RESULTS_DIR" "$file"
fi

#-t 1 : use 1 thread
#--bind : bind mount directories from host to container
#"$file" : input fastq file to process
#"$SIF_PATH" : path to the Apptainer container with FastQC installed
#"$FASTQ_PATTERN" : pattern to match fastq files
#"$fastp_trimming_dir" : directory containing the input fastq files
#-o "$RESULTS_DIR" : specify output directory for FastQC results




