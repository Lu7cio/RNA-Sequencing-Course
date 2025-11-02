#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cwpu=1000M
#SBATCH --time=00:01:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_anlaysis
#SBATCH --mail-user=mario.kummer@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/mkummer/RNA-Sequencing-Course/results/fastqc_analysis_quality_check_%j.o
#SBATCH --error=/data/users/mkummer/RNA-Sequencing-Course/results/logs/fastqc_analysis_error_%j.e
