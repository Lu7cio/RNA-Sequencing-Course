#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/hisat2_GRCm39/samtools/sort/samtools_sort_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/hisat2_GRCm39/samtools/sort/samtools_sort_error_%j.e
#SBATCH --array=1-15

#Define raw data paths, continers path and output paths and raw data pattern
input_bam_files="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/samtools/converted_bam"
RESULTS_DIR="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/samtools/sorted_bam"
SIF_PATH="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
SAM_PATTERN_1="*.bam"

# Ensure results directory exists
mkdir -p "$RESULTS_DIR"

# Get list of files and select one for this array task
files=("$input_bam_files"/$SAM_PATTERN_1)
sample="${files[$SLURM_ARRAY_TASK_ID - 1]}"

#Output sample name for checking
echo "Processing sample: $sample"

# Run samtools inside container for sorting bam
apptainer exec \
  --bind /data/users/mkummer/RNA-Sequencing-Course:/data/users/mkummer/RNA-Sequencing-Course \
  "$SIF_PATH" \
  samtools sort -@ 4 $sample -o $RESULTS_DIR/$(basename "$sample" .bam)_sorted.bam











