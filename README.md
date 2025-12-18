# RNA-Sequencing Project

## Author
Mario Kummer

## Course
RNA-Sequencing Course, Autumn Semester 2025

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
  
#### Download of reference genome from Ensembl FTP (https://www.ensembl.org/info/data/ftp/index.html)
The reference genome Mus musculus GRCm39 version 115 was downloaded from the Ensemble
ftp website. From the site the file Mus musculus.GRCm39.dna.primary assembly.fa.gz (Un-
der section DNA(FASTA)) and for the annotation file named Mus musculus.GRCm39.115.gtf.gz
(Under section gene sets) was downloaded. 
  
## Workflow of analysis:

### Workflow steps Summary:
The Workflow was split into to major parts, the first part (Step 1 -5) was done in bash scripts for submitting it to SLURM on a HPC-Cluster and the second part (Step 6 - 8) was done in a R script locally. For running the bash scripts on the HPC cluster (With SLURM) the command ```bash sbatch my_script.sh```  followed by the name of the desired script was used.
    1. FastQC quality control on raw lung reads (Bash)
    2. FastP trimming of the raw lung reads (Bash)
    3. FastQC quality control on the trimmed lung reads (Bash)
    4. Hisat2 maping reads to reference genome / SAMtools sam/bam conversion and indexing (Bash)
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

4. Hisat2 mapping reads to reference genome & SAMtools sam/bam conversion and sorting/indexing:
    For this part we need to download the latest reference genome sequence and associated annotation for Mus musculus GRCm39 (Download instruction, see section above):
    
    
    - First we have to index the reference genome with Hisat2. The Hisat2 indexing was run on the \*.fa and \*.gtf file:
    
        --> Used container/apptainer: hisat2_samtools_408dfd02f175cd88.sif (Hisat2 version 2.2.1)
    
        --> Used script: hisat2_indexing.sh    (In folder hisat2_scripts)
    
    
    - Second we map for each sample separately,the fastp trimmed reads to the reference genome using Hisat2. 
      In the same step we converted the resulting *.sam files from the mapping to readable *.bam files with SAMtools:
    
        --> Used container/apptainer: hisat2_samtools_408dfd02f175cd88.sif (Hisat2 version 2.2.1 & SAMtools version 1.20)
    
        --> Used script: hisat2_mapping_trimmed.sh    (In folder hisat2_scripts)
    
    
    - Third we sorted the *.bam files by genomic coordinates using Samtools:
    
        --> Used container/apptainer: hisat2_samtools_408dfd02f175cd88.sif (SAMtools version 1.20)
    
        --> Used script: samtools_index.sh    (In folder samtools_scripts)
      
    
    - Forth we index the coordinate sorted bam files using SAMtools:
    
        --> Used container/apptainer: hisat2_samtools_408dfd02f175cd88.sif (SAMtools version 1.20)
    
        --> Used script: samtools_sort.sh    (In folder samtools_scripts)
    
5. FeatureCounts count reads per gen:

    The FeatureCounts was run on the sorted *.bam files of the sorted mapped reads to the referece Genome (Mus_musculus.GRCm39). 

      --> Used container/apptainer: subread_2.0.6.sif (Subread version 0.12.1)

      --> Used script: feature_counts.sh    (In folder featureCounts_scripts)

6.- 8. DESeq2 Differential expression analysis and clusterProfiler Go terms enrichment:

  This part was all done in one R script. For each of the steps 6. - 8. an section was created in the script. 
  When you want only one part execute of thes 3 step execute onyl this part in the R script (Pay attention, 
  some variables or calculations depends on previous step!). It is adviced to run the whole script in once, 
  so the correctness is ensured and you have all the necessary steps performed. 
  The steps contains:
    
  6. Exploratory data analysis (How the samples cluster based on their gene expression profiles --> PCA Plot).
  7. Differential expression analysis (Differential Expression for one pairwise contrast, i.e. the comparison between two experimental groups).
  8. Overrepresentation analysis (Identification Gene Ontology terms that contain more differentially expressed genes than expected by chance for the choosen pairwise comparison)


### Scripts folder structure
```bash
scripts/
├── DESeq2_scripts/                         # Folder for all DESeq2 analysis (R)
│   ├── DESeq2_analysis.R                   # Script for DESeq2 analysis --> Workflow 5.- 7.
│ 
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


