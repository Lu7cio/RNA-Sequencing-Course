#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/featureCounts/featureCounts_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/featureCounts/featureCounts_error_%j.e
#SBATCH --array=1-15

#Define raw data paths, continers path and output paths and raw data pattern
input_bam_files="/data/users/mkummer/RNA-Sequencing-Course/results/samtools/sorted_bam"
input_annotation_file="/data/users/mkummer/RNA-Sequencing-Course/reference_genome/Mus_musculus.GRCm39.115.gtf"
RESULTS_DIR="/data/users/mkummer/RNA-Sequencing-Course/results/featureCounts"
SIF_PATH="/containers/apptainer/subread_2.0.6.sif"
BAM_PATTERN_1="*_sorted.bam"


# Ensure results directory exists
mkdir -p "$RESULTS_DIR"

# Get list of R1 files and select one for this array task
files=("$input_bam_files"/$BAM_PATTERN_1)
sample="${files[$SLURM_ARRAY_TASK_ID - 1]}"

#Output sample name for checking
echo "Processing sample: $sample"

# Run featureCounts inside container for counting number of reads per gene
apptainer exec \
  --bind /data/users/mkummer/RNA-Sequencing-Course:/data/users/mkummer/RNA-Sequencing-Course \
  "$SIF_PATH" \
  featureCounts \
    -p -B -C \
    -T 4 \
    -s 2 \
    -a "$input_annotation_file" \
    -o "$RESULTS_DIR"/$(basename "$sample" .bam).txt \
    "$sample"