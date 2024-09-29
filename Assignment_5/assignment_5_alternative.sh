#!/bin/bash

# Download and extract genome
wget ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

# Report file size
echo "File size: $(du -sh Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa)"

# Total genome size
genome_size=$(grep -v ">" Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa | wc -c)
echo "Total genome size: $genome_size bases"

# Number of chromosomes
chromosomes=$(grep ">" Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa | wc -l)
echo "Number of chromosomes: $chromosomes"

# Chromosome lengths
awk '/^>/{if (seqlen){print seqlen}; printf substr($0, 2, 30) " "; seqlen=0; next}{seqlen += length($0)}END{print seqlen}' Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

# Simulate FASTQ with wgsim (assuming installed)
wgsim -N 240000 -1 150 -2 150 Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa simulated_reads_1.fastq simulated_reads_2.fastq

# Report FASTQ size and compress
du -sh simulated_reads_1.fastq simulated_reads_2.fastq
gzip simulated_reads_1.fastq simulated_reads_2.fastq
du -sh simulated_reads_1.fastq.gz simulated_reads_2.fastq.gz