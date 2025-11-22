#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/hisat2_GRCm39/samtools/index/samtools_index_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/hisat2_GRCm39/samtools/index/samtools_index_error_%j.e
#SBATCH --array=1-15

#Define sorted bam file path, containers path and output paths and raw data pattern
input_bam_files="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/samtools/sorted_bam"
RESULTS_DIR="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/samtools/indexed_bam"
SIF_PATH="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
SAM_PATTERN_1="*.bam"

# Ensure results directory exists
mkdir -p "$RESULTS_DIR"

# Get list of files and select one for this array task
files=("$input_bam_files"/$SAM_PATTERN_1)
sample="${files[$SLURM_ARRAY_TASK_ID - 1]}"

#Output sample name for control
echo "Processing sample: $sample"

# Run samtools inside container for indexing bam
apptainer exec \
  --bind /data/users/mkummer/RNA-Sequencing-Course:/data/users/mkummer/RNA-Sequencing-Course \
  "$SIF_PATH" \
  samtools index -@ 1 $sample -o $RESULTS_DIR/$(basename "$sample" _sorted.bam)_indexed_sorted.bam.bai

#-@ 1 : use 1 thread
#samtools index : index bam file
#--bind : bind mount directories from host to container
#"$sample" : input sorted bam file to index
#"$SIF_PATH" : path to the Apptainer container with samtools installed
#"$RESULTS_DIR/$(basename "$sample" _sorted.bam)_indexed_sorted.bam.bai : output indexed bam file path   
#$(basename "$sample" _sorted.bam) : extract sample name without _sorted.bam suffix 










