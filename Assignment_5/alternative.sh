#!/bin/bash

# Set error handling and trace
set -uex

# ----- DEFINITIONS -----

# The URL of the genome file
GENOME_URL="ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"

# The name of the genome file after download
GENOME_FILE_GZ="Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
GENOME_FILE="Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"

# Simulated FASTQ file names
SIM_READS_1="simulated_reads_1.fastq"
SIM_READS_2="simulated_reads_2.fastq"
SIM_READS_1_GZ="${SIM_READS_1}.gz"
SIM_READS_2_GZ="${SIM_READS_2}.gz"

# ----- ACTIONS -----

# Step 1: Download the genome file
echo "Downloading genome..."
wget ${GENOME_URL}

# Step 2: Unzip the genome file
echo "Unzipping genome..."
gunzip ${GENOME_FILE_GZ}

# Step 3: Report the size of the genome file
echo "Reporting file size..."
du -sh ${GENOME_FILE}

# Step 4: Calculate and report the total genome size (in bases)
genome_size=$(grep -v ">" ${GENOME_FILE} | wc -c)
echo "Total genome size: $genome_size bases"

# Step 5: Report the number of chromosomes (sequences)
chromosomes=$(grep ">" ${GENOME_FILE} | wc -l)
echo "Number of chromosomes: $chromosomes"

# Step 6: Report the name (ID) and length of each chromosome
echo "Chromosome lengths:"
awk '/^>/{if (seqlen){print seqlen}; printf substr($0, 2, 30) " "; seqlen=0; next}{seqlen += length($0)}END{print seqlen}' ${GENOME_FILE}

# Step 7: Simulate FASTQ reads using wgsim
echo "Simulating FASTQ reads..."
wgsim -N 240000 -1 150 -2 150 ${GENOME_FILE} ${SIM_READS_1} ${SIM_READS_2}

# Step 8: Report the size of the simulated FASTQ files
echo "Reporting FASTQ file sizes..."
du -sh ${SIM_READS_1} ${SIM_READS_2}

# Step 9: Compress the FASTQ files
echo "Compressing FASTQ files..."
gzip ${SIM_READS_1} ${SIM_READS_2}

# Step 10: Report the size of the compressed FASTQ files
echo "Reporting compressed file sizes..."
du -sh ${SIM_READS_1_GZ} ${SIM_READS_2_GZ}
