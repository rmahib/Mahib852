#!/bin/bash

# Set the error handling and trace
set -uex

# Define all variables at the top.

# URL to download files
GFF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/771/035/GCA_009771035.1_CNRS_Lalb_1.0/GCA_009771035.1_CNRS_Lalb_1.0_genomic.gff.gz"
FNA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/771/035/GCA_009771035.1_CNRS_Lalb_1.0/GCA_009771035.1_CNRS_Lalb_1.0_genomic.fna.gz"

# Filenames
GFF_FILE="GCA_009771035.1_CNRS_Lalb_1.0_genomic.gff"
FNA_FILE="GCA_009771035.1_CNRS_Lalb_1.0_genomic.fna"
UNZIPPED_GFF="GCA_009771035.1_CNRS_Lalb_1.0_genomic.gff"
UNZIPPED_FNA="GCA_009771035.1_CNRS_Lalb_1.0_genomic.fna"

# Download data into files
wget ${GFF_URL}
wget ${FNA_URL}

# Unzip the files
gunzip ${GFF_FILE}.gz
gunzip ${FNA_FILE}.gz

# Make new GFF files with specific features
grep -w 'exon' ${UNZIPPED_GFF} > lupinus_exons.gff
grep -w 'mRNA' ${UNZIPPED_GFF} > lupinus_mRNA.gff
grep -w 'gene' ${UNZIPPED_GFF} > lupinus_genes.gff
grep -w 'CDS' ${UNZIPPED_GFF} > lupinus_CDS.gff
grep -w 'biological_region' ${UNZIPPED_GFF} > lupinus_biological_region.gff
grep -w 'transcript' ${UNZIPPED_GFF} > lupinus_transcripts.gff
grep -w 'rRNA' ${UNZIPPED_GFF} > lupinus_rRNA.gff
grep -w 'snRNA' ${UNZIPPED_GFF} > lupinus_snRNA.gff
grep -w 'ncRNA_gene' ${UNZIPPED_GFF} > lupinus_ncRNA_gene.gff

# Print the number of each feature
echo "Number of exons:"
grep -w 'exon' ${UNZIPPED_GFF} | wc -l

echo "Number of mRNA:"
grep -w 'mRNA' ${UNZIPPED_GFF} | wc -l

echo "Number of genes:"
grep -w 'gene' ${UNZIPPED_GFF} | wc -l

echo "Number of CDS:"
grep -w 'CDS' ${UNZIPPED_GFF} | wc -l

echo "Number of biological regions:"
grep -w 'biological_region' ${UNZIPPED_GFF} | wc -l

echo "Number of transcripts:"
grep -w 'transcript' ${UNZIPPED_GFF} | wc -l

echo "Number of rRNA:"
grep -w 'rRNA' ${UNZIPPED_GFF} | wc -l

echo "Number of snRNA:"
grep -w 'snRNA' ${UNZIPPED_GFF} | wc -l

echo "Number of ncRNA_gene:"
grep -w 'ncRNA_gene' ${UNZIPPED_GFF} | wc -l

# Print the number of protein coding genes
echo "Number of protein-coding genes:"
awk '$3 == "gene" && $9 ~ /biotype=protein_coding/ {print}' ${UNZIPPED_GFF} | wc -l
