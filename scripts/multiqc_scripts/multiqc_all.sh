#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_analysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/logs/all_multiqc/all_multiqc_analysis__featureCounts_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/logs/all_multiqc/all_multiqc_analysis__featureCounts_error_%j.e

#Define data paths, containers path and output paths
input_fastqc_dir="/data/users/mkummer/RNA-Sequencing-Course/results/fastqc"
input_fastp_trimmed_fastqc_dir="/data/users/mkummer/RNA-Sequencing-Course/results/fastqc_fastp_trimmed"
input_hisat2_dir="/data/users/mkummer/RNA-Sequencing-Course/results/hisat2_GRCm39/samtools/converted_bam"
input_featureCounts_dir="/data/users/mkummer/RNA-Sequencing-Course/results/featureCounts"
output_all_multiqc_dir="/data/users/mkummer/RNA-Sequencing-Course/results/all_workflow_multiqc"
SIF_PATH="/containers/apptainer/multiqc-1.19.sif"

#Ensure results directory exists
mkdir -p "$output_all_multiqc_dir"

#Print input directories for control
echo "$input_fastqc_dir"
echo "$input_fastp_trimmed_fastqc_dir"
echo "$input_hisat2_dir"
echo "$input_featureCounts_dir"
echo "$output_all_multiqc_dir"

# Run MultiQC on featureCounts results 

apptainer exec \
  --bind "$input_fastqc_dir":"$input_fastqc_dir" \
  --bind "$input_fastp_trimmed_fastqc_dir":"$input_fastp_trimmed_fastqc_dir" \
  --bind "$input_hisat2_dir":"$input_hisat2_dir" \
  --bind "$input_featureCounts_dir":"$input_featureCounts_dir" \
  --bind "$output_all_multiqc_dir":"$output_all_multiqc_dir" \
  --bind /data/users/mkummer/RNA-Sequencing-Course/multiqc_config_yaml:/data/users/mkummer/RNA-Sequencing-Course/multiqc_config_yaml \
  "$SIF_PATH" multiqc \
    -c "/data/users/mkummer/RNA-Sequencing-Course/multiqc_config_yaml/multiqc_config.yaml" \
    -o "$output_all_multiqc_dir" \
    "$input_fastqc_dir" "$input_fastp_trimmed_fastqc_dir" "$input_hisat2_dir" "$input_featureCounts_dir"

