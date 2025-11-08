#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/fastp_trimming/fastp_trimming_quality_check_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/fastp_trimming/fastp_trimming_error_%j.e
#SBATCH --array=1-30

#Define raw data paths, continers path and output paths and raw data pattern
input_fastq_dir="/data/courses/rnaseq_course/toxoplasma_de/reads_Lung"
output_fastp_dir="/data/users/mkummer/RNA-Sequencing-Course/results/fastp_trimming"
SIF_PATH="/containers/apptainer/fastp_0.24.1.sif"
FASTQ_PATTERN="*.fastq.*"


#Ensure results directory exists
mkdir -p "$output_fastp_dir"

# Get list of files and select the one for this array task for parallel processing
files=("$input_fastq_dir"/$FASTQ_PATTERN)
file="${files[$SLURM_ARRAY_TASK_ID - 1]}"


if [ -f "$file" ]; then
  echo "Processing $file..."
  filename=$(basename "$file" .fastq.gz)
  apptainer exec \
    --bind "$input_fastq_dir":"$input_fastq_dir" \
    --bind "$output_fastp_dir":"$output_fastp_dir" \
    "$SIF_PATH" fastp \
      -i "$file" \
      -o "$output_fastp_dir/${filename}_trimmed.fastq.gz" \
      --html "$output_fastp_dir/${filename}_fastp.html" \
      --json "$output_fastp_dir/${filename}_fastp.json" \
      --detect_adapter_for_pe \
      --thread 1
else
  echo "No file found for task ID $SLURM_ARRAY_TASK_ID"
fi