#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/samtools/convert/samtools_convert_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/samtools/convert/samtools_convert_error_%j.e
#SBATCH --array=1-15

#Define raw data paths, continers path and output paths and raw data pattern
input_sam_files="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/mapping/"
RESULTS_DIR="/data/users/mkummer/RNA-Sequencing-Course/results/samtools/converted_bam"
SIF_PATH="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
SAM_PATTERN_1="*.sam"


# Ensure results directory exists
mkdir -p "$RESULTS_DIR"

# Get list of R1 files and select one for this array task
files=("$input_sam_files"/$SAM_PATTERN_1)
sample="${files[$SLURM_ARRAY_TASK_ID - 1]}"

#Output sample name for checking
echo "Processing sample: $sample"

# Run samtools inside container for conversion sam to bam
apptainer exec \
  --bind /data/users/mkummer/RNA-Sequencing-Course:/data/users/mkummer/RNA-Sequencing-Course \
  "$SIF_PATH" \
  samtools view -bS $sample > $RESULTS_DIR/$(basename "$sample" .sam).bam