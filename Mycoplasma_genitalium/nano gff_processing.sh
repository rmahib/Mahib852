##!/bin/bash
 #
 ## Define variables for input files and outputs
 #GFF_FILE="GCF_000027325.1_ASM2732v1_genomic.gff"
 #GENE_OUTPUT="mycoplasma_genes.gff"
 #CDS_OUTPUT="mycoplasma_CDS.gff"
 #INTERGENIC_OUTPUT="intergenic_regions.gff"
 #
 ## Step 1: Extract gene features
 #echo "Extracting gene features..."
 #awk '$3 == "gene"' $GFF_FILE > $GENE_OUTPUT
 #
 ## Step 2: Extract CDS features
 #echo "Extracting CDS features..."
 #awk '$3 == "CDS"' $GFF_FILE > $CDS_OUTPUT
 #
 ## Step 3: Identify intergenic regions
 #echo "Generating intergenic regions..."
 #awk '{print $1"\t.\tintergenic_region\t"$2"\t"$3"\t.\t.\t.\tID=intergenic_"NR";Name=intergenic_region_"NR}' gene_intervals.bed > $INTERGENIC_OUTPUT
 #
 ## Step 4: Summarize gene details
 #echo "Summarizing gene details..."
 #awk '$3 == "gene" && $9 ~ /biotype=protein_coding/ {match($9, /ID=[^;]*/, a); match($9, /Name=[^;]*/, b); print $1 ": " a[0] (b[0] ? "; " b[0] : "") "; biotype=protein_coding"}' $GENE_OUTPUT > protein_coding_genes_summary.txt
 #
 ## Step 5: Count number of gene and CDS entries
 #GENE_COUNT=$(awk '$3 == "gene"' $GFF_FILE | wc -l)
 #CDS_COUNT=$(awk '$3 == "CDS"' $GFF_FILE | wc -l)
 #
 #echo "Number of genes: $GENE_COUNT" > genes_summary.txt
 #echo "Number of CDS: $CDS_COUNT" >> genes_summary.txt
 #
 ## Step 6: Output results
 #echo "Script completed. Gene and CDS counts are saved in genes_summary.txt"
 #echo "Intergenic regions saved in $INTERGENIC_OUTPUT"
