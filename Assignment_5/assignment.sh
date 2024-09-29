#!/bin/bash

# Set error handling and trace
set -uex

# ----- DEFINITIONS -----

# URL for downloading the Saccharomyces cerevisiae genome
GENOME_URL="ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"

# Genome file name after download
GENOME_FILE="Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"

# Unzipped genome file name
UNZIPPED_GENOME_FILE="Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"

# Simulated FASTQ file names
SIM_READS_1="simulated_reads_1.fastq"
SIM_READS_2="simulated_reads_2.fastq"
SIM_READS_1_GZ="${SIM_READS_1}.gz"
SIM_READS_2_GZ="${SIM_READS_2}.gz"

# Number of reads and read lengths for wgsim
NUM_READS=240000
READ_LENGTH_1=150
READ_LENGTH_2=150

# Path to wgsim binary after build
WGSIM_PATH="./wgsim"

# ----- ACTIONS -----

# Install dependencies
sudo apt-get install -y gcc zlib1g-dev

# Clone wgsim repository
git clone https://github.com/lh3/wgsim.git

# Change directory to wgsim
cd wgsim

# Build wgsim binary
make

# Download Saccharomyces cerevisiae genome
wget ${GENOME_URL}

# Unzip the genome file
gunzip ${GENOME_FILE}

# Check file size of the unzipped genome
du -sh ${UNZIPPED_GENOME_FILE}

# Count the number of characters in the genome sequence (excluding headers)
grep -v ">" ${UNZIPPED_GENOME_FILE} | wc -c

# Count the number of sequences (headers)
grep ">" ${UNZIPPED_GENOME_FILE} | wc -l

# Extract sequence names and lengths
awk '/^>/{if (seqlen){print seqlen}; printf substr($0, 2, 30) " "; seqlen=0; next}{seqlen += length($0)}END{print seqlen}' ${UNZIPPED_GENOME_FILE}

# Generate simulated reads using wgsim
${WGSIM_PATH} -N ${NUM_READS} -1 ${READ_LENGTH_1} -2 ${READ_LENGTH_2} ${UNZIPPED_GENOME_FILE} ${SIM_READS_1} ${SIM_READS_2}

# Count the number of reads in the first simulated FASTQ file
grep -c "^@" ${SIM_READS_1}

# Calculate average read length for the first FASTQ file (compressed)
zcat ${SIM_READS_1_GZ} | awk 'NR%4==2 {sum+=length($0); count++} END {print "Average read length for " FILENAME ": " sum/count " bp"}'

# Calculate average read length for the second FASTQ file (compressed)
zcat ${SIM_READS_2_GZ} | awk 'NR%4==2 {sum+=length($0); count++} END {print "Average read length for " FILENAME ": " sum/count " bp"}'

# Check file size of simulated reads before compression
du -sh ${SIM_READS_1} ${SIM_READS_2}

# Compress the simulated reads
gzip ${SIM_READS_1} ${SIM_READS_2}

# Check file size of compressed simulated reads
du -sh ${SIM_READS_1_GZ} ${SIM_READS_2_GZ}


