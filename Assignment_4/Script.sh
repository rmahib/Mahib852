#!/bin/bash

# Set the error handling and trace
set -uex

# Define all variables at the top

# URLs for downloading the files
URL_FNA="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz"
URL_GFF="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz"

# Names of the files
FNA_FILE="GCF_000027325.1_ASM2732v1_genomic.fna.gz"
GFF_FILE="GCF_000027325.1_ASM2732v1_genomic.gff.gz"
UNZIPPED_FNA_FILE="GCF_000027325.1_ASM2732v1_genomic.fna"
UNZIPPED_GFF_FILE="GCF_000027325.1_ASM2732v1_genomic.gff"

# Output files for extracted features
GENES_FILE="mycoplasma_genes.gff"
CDS_FILE="mycoplasma_CDS.gff"
MRNA_FILE="mycoplasma_mRNA.gff"
EXONS_FILE="mycoplasma_exons.gff"

# ------ ALL THE ACTIONS FOLLOW ------

# Download the genomic files
wget ${URL_FNA}
wget ${URL_GFF}

# Unzip the files
gunzip ${FNA_FILE}
gunzip ${GFF_FILE}

# Extract genes from the GFF file
cat ${UNZIPPED_GFF_FILE} | awk '$3 == "gene"' > ${GENES_FILE}

# Extract CDS from the GFF file
cat ${UNZIPPED_GFF_FILE} | awk '$3 == "CDS"' > ${CDS_FILE}

# Extract mRNA from the GFF file
cat ${UNZIPPED_GFF_FILE} | awk '$3 == "mRNA"' > ${MRNA_FILE}

# Extract exons from the GFF file
cat ${UNZIPPED_GFF_FILE} | awk '$3 == "exon"' > ${EXONS_FILE}

# Print the number of genes
echo "Number of genes:"
cat ${GENES_FILE} | wc -l

# Print the number of CDS
echo "Number of CDS:"
cat ${CDS_FILE} | wc -l

# Print the number of mRNA
echo "Number of mRNA:"
cat ${MRNA_FILE} | wc -l

# Print the number of exons
echo "Number of exons:"
cat ${EXONS_FILE} | wc -l

# Optional: Show the first few lines of each output file
echo "First few lines of genes file:"
head ${GENES_FILE}

echo "First few lines of CDS file:"
head ${CDS_FILE}

echo "First few lines of mRNA file:"
head ${MRNA_FILE}

echo "First few lines of exons file:"
head ${EXONS_FILE}
