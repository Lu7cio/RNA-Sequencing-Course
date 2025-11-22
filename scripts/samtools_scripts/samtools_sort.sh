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

#Define bam file path, continers path and output paths and bam pattern
input_bam_files="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/samtools/converted_bam"
RESULTS_DIR="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/samtools/sorted_bam"
SIF_PATH="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
BAM_PATTERN_1="*.bam"

# Ensure results directory exists
mkdir -p "$RESULTS_DIR"

# Get list of files and select one for this array task
files=("$input_bam_files"/$BAM_PATTERN_1)
sample="${files[$SLURM_ARRAY_TASK_ID - 1]}"

#Output sample name for control
echo "Processing sample: $sample"

# Run samtools inside container for sorting bam
apptainer exec \
  --bind /data/users/mkummer/RNA-Sequencing-Course:/data/users/mkummer/RNA-Sequencing-Course \
  "$SIF_PATH" \
  samtools sort -@ 4 $sample -o $RESULTS_DIR/$(basename "$sample" .bam)_sorted.bam

#-@ 4 : use 4 threads
#samtools sort : sort bam file
#--bind : bind mount directories from host to container
#"$sample" : input bam file to sort
#"$SIF_PATH" : path to the Apptainer container with samtools installed
#"$RESULTS_DIR/$(basename "$sample" .bam)_sorted.bam : output
#$(basename "$sample" .bam) : extract sample name without .bam suffix












