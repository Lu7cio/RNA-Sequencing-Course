#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/hisat2_GRCm39/mapping/trimmed/hisat2_mapping_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/hisat2_GRCm39/mapping/trimmed/hisat2_mapping_error_%j.e
#SBATCH --array=1-15

#Define raw data paths, continers path and output paths and raw data pattern
input_fastq_trimmed_dir="/data/users/mkummer/RNA-Sequencing-Course/results/fastp_trimming"
input_index_reference="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/indexing/GRCm39"
RESULTS_DIR="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/mapping"
SIF_PATH="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
FASTQ_PATTERN_1="*_1_trimmed.fastq.gz"
FASTQ_PATTERN_2="*_2_trimmed.fastq.gz"

# Ensure results directory exists
mkdir -p "$RESULTS_DIR"

# Get list of R1 files and select one for this array task
files_R1=("$input_fastq_trimmed_dir"/$FASTQ_PATTERN_1)
R1="${files_R1[$SLURM_ARRAY_TASK_ID - 1]}"
files_R2=("$input_fastq_trimmed_dir"/$FASTQ_PATTERN_2)
R2="${files_R2[$SLURM_ARRAY_TASK_ID - 1]}"
R1_FILENAME=$(basename "$R1")
SAMPLE=$(basename "$R1_FILENAME" | sed 's/_1_trimmed.fastq.gz//' | sed 's/_2_trimmed.fastq.gz//')

#Output sample name for checking
echo "Processing sample: $SAMPLE"
echo "R1 file: $R1"
echo "R2 file: $R2"


# Run HISAT2 inside container
apptainer exec \
  --bind /data/users/mkummer/RNA-Sequencing-Course:/data/users/mkummer/RNA-Sequencing-Course \
  "$SIF_PATH" \
  hisat2 \
    -x "$input_index_reference" \
    -1 "$R1" \
    -2 "$R2" \
    --rna-strandness RF \
    --dta \
    -p 4 \
    -S "$RESULTS_DIR/${SAMPLE}.sam"





