# RNA-Seq Analysis Report

## Overview
This report documents the RNA-Seq analysis performed on the Saccharomyces cerevisiae dataset using three replicates: `SRR19383416`, `SRR19383417`, and `SRR19383418`. The analysis includes data acquisition, alignment, quantification, and visualization to assess gene expression levels.

---

## Workflow Explanation
1. **Data Acquisition**: Raw RNA-Seq data for three replicates was downloaded from SRA.
2. **Genome Preparation**: The reference genome and corresponding annotation (GTF) were obtained from Ensembl.
3. **Alignment**: Reads were aligned to the reference genome using HISAT2, and SAM files were converted to sorted BAM files using samtools.
4. **Quantification**: FeatureCounts was used to generate a count matrix that summarizes read counts for each gene in each replicate.
5. **Visualization**: IGV was employed to visually confirm the alignment of reads to the genome.
6. **Data Analysis**: R was used for data normalization, identification of consistently expressed genes, and exploratory data analysis.

---

## Methods

### 1. Data Acquisition
Three replicates from Saccharomyces cerevisiae were downloaded using `fasterq-dump`. A CSV file containing metadata about the samples was generated:
```bash
echo -e "sample,accession\nreplicate_1,SRR19383416\nreplicate_2,SRR19383417\nreplicate_3,SRR19383418" > design.csv
```

2. Genome Reference
The reference genome and GTF file for Saccharomyces cerevisiae were downloaded from Ensembl:

````
wget ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-105/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.105.gtf.gz
gunzip *.gz
````

3. Alignment
Reads were aligned to the reference genome using HISAT2:

````
hisat2 -x genome/scerevisiae_index -1 fastq/SRR19383416_1.fastq -2 fastq/SRR19383416_2.fastq -S alignments/SRR19383416.sam
````

This step was repeated for all replicates. SAM files were converted to sorted BAM files:

````
samtools sort -o alignments/sorted/SRR19383416_sorted.bam alignments/SRR19383416.sam
````

4. Quantification
Gene-level counts were obtained using FeatureCounts:

````
featureCounts -a genome/Saccharomyces_cerevisiae.R64-1-1.105.gtf -o counts/gene_counts.txt alignments/sorted/*.bam
````

5. Visualization

Aligned reads were visualized in IGV. Several regions of the genome were inspected to confirm the accuracy of RNA-Seq alignments. IGV screenshots are included in the results.

6. Data Analysis

````
# Load count data
count_data <- read.table("counts/gene_counts.txt", header = TRUE, row.names = 1)
count_data <- count_data[, -c(1:5)] # Remove unnecessary columns

# Normalize data
normalized_data <- t(t(count_data) / colSums(count_data)) * 1e6

# Identify consistently expressed genes
consistent_data <- normalized_data[rowSums(normalized_data > 1) == ncol(normalized_data), ]

# Save consistent genes
write.table(consistent_data, "counts/consistent_genes.txt", sep = "\t", quote = FALSE)
````

Results
1. Count Matrix Summary
The count matrix was generated successfully and contains read counts for all samples. A subset of the data is shown below:

````
## Count Matrix Results

GeneID     SRR19383416 SRR19383417 SRR19383418
YDL246C    2           8           4
YDL243C    139         191         182
YDR387C    581         670         643
YDL094C    137         188         162
YDR438W    273         366         327
YDR523C    60          84          28
YDR542W    0           0           4
YDR492W    964         1328        1033
````

2. Visualization
IGV confirmed the alignment of reads to the reference genome. Below are the screenshots of RNA-Seq data visualized on IGV:

Chromosome VIII Region

Chromosome XI Region

Chromosome VII Region

<img width="842" alt="Assignment_13-9" src="https://github.com/user-attachments/assets/599e6987-968a-4c42-9dbf-097404038b5b">

<img width="851" alt="Assignment 13_10" src="https://github.com/user-attachments/assets/02f9785b-1188-4507-bd70-862b88370398">

<img width="850" alt="Assignment_13_11" src="https://github.com/user-attachments/assets/27702be2-6fb4-41f3-9264-02ebba9b4a00">


3. Consistently Expressed Genes

Consistently expressed genes across all samples were identified:
````
GeneID     SRR19383416 SRR19383417 SRR19383418
YDR438W    273         366         327
YDR523C     60          84          28
````


*Critical Observations*

## Coverage: 

The IGV screenshots confirm robust coverage for exonic regions, validating the RNA-Seq dataset's quality.

## Gene Expression: 

The count matrix shows variability in gene expression across samples, with some genes exhibiting high expression levels, indicating active transcription.

## Consistency: 

Consistently expressed genes were observed across replicates, providing confidence in the experimental reproducibility.


## Discussion

This analysis demonstrates successful RNA-Seq data processing, from alignment to feature quantification and visualization. The results highlight gene expression patterns in Saccharomyces cerevisiae and provide insights into transcriptional activity.
