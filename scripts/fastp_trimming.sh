#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/fastp_trimming/fastp_trimming_quality_check_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/fastp_trimming/fastp_trimming_error_%j.e
#SBATCH --array=1-15

#Define raw data paths, continers path and output paths and raw data pattern
input_fastq_dir="/data/courses/rnaseq_course/toxoplasma_de/reads_Lung"
output_fastp_dir="/data/users/mkummer/RNA-Sequencing-Course/results/fastp_trimming"
SIF_PATH="/containers/apptainer/fastp_0.24.1.sif"
FASTQ_PATTERN_1="*_1.fastq.gz"
FASTQ_PATTERN_2="*_2.fastq.gz"

#Ensure results directory exists
mkdir -p "$output_fastp_dir"


# Get list of R1,R2 files and select one for this array task
files_R1=("$input_fastq_dir"/$FASTQ_PATTERN_1)
R1="${files_R1[$SLURM_ARRAY_TASK_ID - 1]}"
files_R2=("$input_fastq_dir"/$FASTQ_PATTERN_2)
R2="${files_R2[$SLURM_ARRAY_TASK_ID - 1]}"
R1_FILENAME=$(basename "$R1")
SAMPLE=$(basename "$R1_FILENAME" | sed 's/_1.fastq.gz//' | sed 's/_2.fastq.gz//')
filename="$SAMPLE"

#Output sample name for checking
echo "Processing sample: $SAMPLE"
echo "R1 file: $R1"
echo "R2 file: $R2"

# Run fastp trimming inside container
  apptainer exec \
    --bind "$input_fastq_dir":"$input_fastq_dir" \
    --bind "$output_fastp_dir":"$output_fastp_dir" \
    "$SIF_PATH" fastp \
      -i "$R1" \
      -I "$R2" \
      -o "$output_fastp_dir/${filename}_1_trimmed.fastq.gz" \
      -O "$output_fastp_dir/${filename}_2_trimmed.fastq.gz" \
      --html "$output_fastp_dir/${filename}_fastp.html" \
      --json "$output_fastp_dir/${filename}_fastp.json" \
      --detect_adapter_for_pe \
      --thread 1

