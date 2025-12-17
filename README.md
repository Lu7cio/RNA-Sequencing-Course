# RNA-Sequencing Project

## Author
Mario Kummer

## Course & University
RNA-Sequencing Course Autumn Semester 2025

# University
University of Bern

## Description
This repository contains a complete RNA-seq analysis pipeline for the study of lung immune responses to Toxoplasma gondii infection in wild-type (WT) and double-knockout (Ifnar−/− × Ifngr−/−; DKO) mice. 
The data analyzed is subset from the publication Singhania et al. 2019 (https://doi.org/10.1038/s41467-019-10601-6). 
More precise the analysis is based on the data of 4 different experimental groups (All samples from lung tissues): 
- WT Case (5x): SRR7821918–SRR7821922  
- DKO Case (4x): SRR7821923–SRR7821927  
- WT Control (3x): SRR7821937–SRR7821939  
- DKO Control (3x): SRR7821940–SRR7821942 

The goal was to produce lists of genes that are differentially expressed between two experimental groups, and identify gene ontology (GO) terms enriched for DE genes.


## Software & Tool versions
All analyses were performed using fixed software versions and a reproducible
R environment managed with renv. See below for full version details.

#### Bash-based analysis
| Tool / Software            | Version  |
|----------------------------|---------|
| FastQC                     | 0.12.1  |
| MultiQC                    | 1.32    |
| fastp                       | 24.1    |
| HISAT2                     | 2.2.1   |
| SAMtools                   | 1.20    |
| Subread (featureCounts)    | 2.0.6   |


#### R-based analysis:

- R: version 4.5.1 (2025-06-13) 
- Reproducing the R analysis (to get R environment renv):
  From the project root directory in R:
  ```r
  install.packages("renv")   # if not already installed
  renv::restore()
  ```
  
## Workflow of analysis:

### Workflow Summary steps:
    The Workflow was split into to major parts, the first part (Step 1 -5) was done in bash scripts for submitting it to SLURM on a HPC-Cluster and the second part (Step 6 - 8) was done in a R script locally:
    1. FastQC quality control on raw lung reads (Bash)
    2. FastP trimming of the raw lung reads (Bash)
    3. FastQC quality control on the trimmed lung reads (Bash)
    4. Hisat2 maping reads to reference genome / SamTools sam/bam conversion and indexing (Bash)
    5. FeatureCounts count reads per gen (Bash)
    6. Exploratory data analysis (R)
    7. Differential expression analysis (R)
    8. Overrepresentation analysis (R)

### Detailed Workflow
1. FastQC quality control on raw lung reads:
    This was performed to asses the quality of the raw lung reads, which means we look at the base quality along the read or if there is evidence for adapter sequences. This report serves as decicsion foundation for trimming of the reads.

    The FastQC was run on the *fastq* files of the raw lung data (For both reads for each sample). 

    --> Used container/apptainer: fastqc-0.12.1.sif (FastQC version 0.12.1)

    --> Used script: fastqc_analysis.sh    (In folder fastqc_scripts)



2. FastP trimming of the raw lung reads:
    This was performed to asses the quality of the raw lung reads. Primarly we want to get ride of the adapter seqeunces at the end of the read, and for this we used the option --detect_adapter_for_pe on the FastP trimm.

    The FastP was run on the *fastq* files of the raw lung data (For both reads for each sample). 

    --> Used container/apptainer: fastp_0.24.1.sif (FastP version 0.24.1)

    --> Used script: fastp_trimming.sh    (In folder fastp_scripts)



3. FastQC quality control on the trimmed lung reads:
    This was performed to asses the quality of the trimmed lung reads, which means we look at the base quality along the read or if there is still  evidence for adapter sequences. This report serves as decicsion foundation if the trimming was secssful and if now the quality of the lung reads is good enough for further downstream analysis.

    The FastQC was run on the *fastq* files of the trimmed lung data (For both reads for each sample). 

    --> Used container/apptainer: fastqc-0.12.1.sif (FastQC version 0.12.1)

    --> Used script: fastqc_fastp_trimmed.sh    (In folder fastqc_scripts)


3. FastQC quality control on the trimmed lung reads:
    This was performed to asses the quality of the trimmed lung reads, which means we look at the base quality along the read or if there is still  evidence for adapter sequences. This report serves as decicsion foundation if the trimming was secssful and if now the quality of the lung reads is good enough for further downstream analysis.

    The FastQC was run on the *fastq* files of the trimmed lung data (For both reads for each sample). 

    --> Used container/apptainer: fastqc-0.12.1.sif (FastQC version 0.12.1)

    --> Used script: fastqc_fastp_trimmed.sh    (In folder fastqc_scripts)


4. Hisat2 maping reads to reference genome / SamTools sam/bam conversion and indexing:
    ???


5. FeatureCounts count reads per gen:


    The FeatureCounts was run on the sorted *.bam files of the sorted mapped reads to the referece Genome (Mus_musculus.GRCm39). 

    --> Used container/apptainer: subread_2.0.6.sif (Subread version 0.12.1)

    --> Used script: feature_counts.sh    (In folder featureCounts_scripts)






### Scripts folder structure
```bash
scripts/
├── DESeq2_scripts/                         # Folder for all DESeq2 analysis (R)
│   ├── DESeq2_analysis.R                   # Script for DESeq2 analysis --> Workflow 5.- 7.
├── fastp_scripts/                          # Folder for all FastP scripts (Bash)
│   ├── fastp_trimming.sh                   # Script for running FastP trimming on raw lung reads
│ 
├── fastqc_scripts/                         # Folder for all FastQC scripts (Bash)
│   ├── fastqc_analysis.sh                  # Script for running FastQC on raw lung reads
│   └── fastqc_fastp_trimmed.sh             # Script for running FastQC on trimmed lung reads
│ 
├── featureCounts_scripts/                  # Folder for all FeatureCounts scripts (Bash)
│   └── feature_counts.sh                   # Script for running FeatureCounts on sorted mapped bam files for each sample
│   
├── hisat2_scripts/                         # Folder for all Hisat2 scripts (Bash)
│   ├── hisat2_indexing.sh                  # Script for running Hisat2 indexing on reference genome file (*.fa / *.gtf)
│   └── hisat2_mapping_trimmed.sh           # Script for running Hisat2 mapping of trimmed lung reads to reference genome
│ 
├── multiqc_scripts/                        # Folder for all Multiqc scripts (Bash)
│   ├── multiqc_all.sh                      # Script for running Multiqc on all FastQC reports (Fastqc, Hisat2 featureCounts)
│   ├── multiqc_analysis_fastp_trimmed.sh   # Script for running Multiqc on Fastp trimmed FastQC reports
│   ├── multiqc_analysis.sh                 # Script for running FastQC on raw lung reads FastQC reports
│   ├── multiqc_featureCounts.sh            # Script for running FastQC on FeatureCounts FastQC reports
│   └── multiqc_hisat2_mapping_trimmed.sh   # Script for running FastQC on Hisat2 FastQC reports
│
├── old_scripts/                            # Folder for all old scripts, which where not used in final analysis (Bash)
│   ├── hisat2_mapping.sh                   # Old script for running hisat2 mapping but without direct conversion sam to bam
│   └── samtools_convert.sh                 # Old script for running samtools conversion sam to bam
│
└── samtools_scripts/                       # Folder for all Samtools scripts (Bash)
    ├── sammtools_index.sh                  # Script for running Samtools indexing mapped, sorted bam files
    └── samtools_sort.sh                    # Script for running Samtools sorting for mapped bam files


```


### Config folder structure
```bash
├── config
 └── multiqc_config_yaml
     └── multiqc_config.yaml                # yaml file for custom config for the mmultiqc report for Workflow steps 1-6


```


