# Analysis and Read Trimming with Makefile

## Introduction
This report explains the usage of the Makefile developed to automate the process of genome download, read simulation, SRA data download, read trimming, and quality control analysis.


The Makefile contains the following key tasks:
1. **Genome Download**: Downloads the genome of *Saccharomyces cerevisiae*.
2. **Read Simulation**: Simulates FASTQ reads from the downloaded genome.
3. **SRA Data Download**: Downloads reads from the SRA database for accession number **SRR387901**.
4. **Read Trimming**: Trims the downloaded reads using **Trimmomatic**.
5. **Quality Control (FastQC)**: Performs quality checks on the reads before and after trimming.

## Setup and Usage
### Variables
The following variables are defined in the Makefile to ensure flexibility:
- `SRA_ID`: The accession number of the dataset from SRA.
- `GENOME_URL`: The URL for downloading the yeast genome from Ensembl.
- `TRIMMOMATIC_JAR`: The path to the Trimmomatic JAR file, downloaded if not present.
- `FASTQC_CMD`: Command to run FastQC for quality control.
- `QC_DIR`: Directory for storing FastQC reports.
- `TRIMMED_FILE`: Output file for trimmed reads.
- `IMPROVED_TRIMMED_FILE`: Output file for aggressively trimmed reads.

### How to Run the Makefile
The following targets were run to perform specific tasks in the workflow:

1. **Usage**
   To view the list of available targets in the Makefile:
````
make usage
````

**To download and unzip the Saccharomyces cerevisiae genome:**

````
make genome
````
**To simulate FASTQ reads from the downloaded genome using wgsim:**

````
make simulate
````

**To download reads from the SRA for accession number SRR387901 and run FastQC on the raw data:**

````
make download
````

**To trim the downloaded reads with Trimmomatic, generate FastQC reports before and after trimming:**

````
make trim
````
**To run FastQC on all the trimmed reads and generate reports:**

````
make fastqc
````
**Workflow Steps**
1. Downloading the Genome
The genome is downloaded from the Ensembl FTP server. It is a compressed .gz file that is unzipped for further analysis. The downloaded file is placed in the genome_data directory, and its size and contents are reported.

2. Simulating Reads
Simulated reads are generated from the downloaded genome using wgsim. Two FASTQ files, one for each read pair, are created, and then these files are compressed for efficient storage.

3. Downloading SRA Reads
Reads for the SRA accession number SRR387901 are downloaded using the fastq-dump tool from the SRA Toolkit. After downloading, a FastQC report is generated to evaluate the quality of the raw reads.

4. Trimming the Reads
Trimming is performed in two stages using Trimmomatic:

Initial trimming: Uses moderate parameters to remove low-quality bases and adapters.
Aggressive trimming: Uses stricter parameters to further clean up the reads.
5. Quality Control (FastQC)
FastQC is run on both the raw and trimmed reads to assess the quality at different stages. This provides a visual report of base quality, GC content, adapter contamination, and overrepresented sequences.

## Conclusion
The Makefile provides a streamlined and reproducible workflow for processing sequencing data. It automates the download, trimming, and quality control steps, reducing manual intervention and ensuring consistency. The use of variables in the Makefile makes it easy to adjust parameters such as the genome URL or SRA accession number for different datasets

**Make file Link**

````
https://github.com/rmahib/Mahib852/blob/main/Assignment_7/makefile
````