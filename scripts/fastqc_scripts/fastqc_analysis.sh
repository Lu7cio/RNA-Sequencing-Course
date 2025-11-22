#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/fastqc/fastqc_analysis_quality_check_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/fastqc/fastqc_analysis_error_%j.e
#SBATCH --array=1-30


#Define raw data paths, continers path and output paths and raw data pattern
READS_DIR="/data/courses/rnaseq_course/toxoplasma_de/reads_Lung"
RESULTS_DIR="/data/users/mkummer/RNA-Sequencing-Course/results/fastqc"
SIF_PATH="/containers/apptainer/fastqc-0.12.1.sif"
FASTQ_PATTERN="*.fastq.*"

#Ensure results directory exists
mkdir -p "$RESULTS_DIR"

# Get list of files and select the one for this array task for parallel processing
files=("$input_fastq_dir"/$FASTQ_PATTERN)
file="${files[$SLURM_ARRAY_TASK_ID - 1]}"

#Loop over all raw data fastq files in the reads directory and run FastQC inside the Apptainer container
if [ -f "$file" ]; then
  echo "Processing $file..."
  apptainer exec \
  --bind "$READS_DIR":"$READS_DIR" \
  --bind "$RESULTS_DIR":"$RESULTS_DIR" \
  "$SIF_PATH" fastqc -t 2 -o "$RESULTS_DIR" "$file"
fi


#-t 2 : use 2 threads
#--bind : bind mount directories from host to container
#"$file" : input fastq file to process
#"$SIF_PATH" : path to the Apptainer container with FastQC installed
#"$READS_DIR" : directory containing the input fastq files
#"$FASTQ_PATTERN" : pattern to match fastq files
#"$input_fastq_dir" : directory containing the input fastq files
#-o "$RESULTS_DIR" : specify output directory for FastQC results




