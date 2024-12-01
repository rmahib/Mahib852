# Variant Calling and Annotation Report

## Overview
This report details the workflow and results of a variant calling and annotation pipeline performed on Saccharomyces cerevisiae. The goal was to identify, merge, and annotate genetic variants across multiple samples from the SRA database.

---

## Workflow Description

### 1. **Reference Genome and GFF Annotation**
- **Objective**: Download and prepare the reference genome and gene annotation files.
- **Steps**:
  - Download the FASTA file for the Saccharomyces cerevisiae reference genome.
  - Index the FASTA file using `samtools` for efficient querying during alignment.
  - Download the GFF3 file for annotation, filter out irrelevant chromosome names, and index it using `tabix`.

### 2. **Sample Data Acquisition**
- **Objective**: Download SRA data for the provided samples and convert them to FASTQ format.
- **Steps**:
  - A `design.csv` file lists the SRA accessions to process.
  - The `fasterq-dump` tool fetches and converts the data into FASTQ files.

### 3. **Read Alignment**
- **Objective**: Align the FASTQ reads to the reference genome.
- **Steps**:
  - Use `bwa mem` for alignment and `samtools` to sort and generate BAM files.
  - Each BAM file corresponds to an aligned sample.

### 4. **Variant Calling**
- **Objective**: Call genetic variants for each aligned sample.
- **Steps**:
  - Use `bcftools mpileup` to identify genomic variants.
  - Variants are output in VCF format and compressed with `bgzip`.

### 5. **VCF File Merging**
- **Objective**: Merge individual VCF files into a single VCF file.
- **Steps**:
  - Use `bcftools merge` to combine VCF files across samples.
  - The merged VCF file represents the combined genetic variants across the dataset.

### 6. **Variant Annotation**
- **Objective**: Annotate the merged VCF file to determine the functional impact of variants.
- **Steps**:
  - Use the Ensembl VEP tool with the GFF3 file and reference genome to annotate variants.
  - The output includes details on variant consequences, affected genes, and variant statistics.

---

## Results and Observations

### 1. **Variant Statistics**
- **Total Variants**: 213 variants processed.
- **Variant Classes**:
  - **Single Nucleotide Variants (SNVs)**: 210 (98.6%)
  - **Insertions**: 3 (1.4%)

### 2. **Variants by Chromosome**
- Variants are distributed across multiple chromosomes. The most variants were found on:
  - Chromosome V (49 variants)
  - Chromosome XII (45 variants)

| Chromosome | Variant Count |
|------------|---------------|
| IX         | 2             |
| I          | 4             |
| II         | 10            |
| IV         | 20            |
| V          | 49            |
| VI         | 6             |
| VII        | 3             |
| VIII       | 5             |
| X          | 11            |
| XV         | 21            |

### 3. **Variant Consequences**
- **Synonymous Variants**: 69.8%
- **Missense Variants**: 28.8%
- **Inframe Insertions**: 1.4%

#### Pie Chart of Coding Consequences
- The largest proportion of coding consequences are synonymous variants, followed by missense variants.

### 4. **Most Severe Consequences**
- Based on severity, the majority of variants are synonymous or missense, with minor effects on gene function.

| Consequence Type      | Count |
|-----------------------|-------|
| Synonymous Variant    | 142   |
| Missense Variant      | 58    |
| Inframe Insertion     | 3     |
| Upstream Gene Variant | 5     |
| Downstream Gene Variant | 4   |

---

## Key Findings
- **Variant Distribution**: Chromosomes V and XII had the most variants, indicating potential areas of interest for further analysis.
- **Functional Impact**: Most variants are synonymous, suggesting minimal functional impact. However, the missense variants and inframe insertions could potentially alter protein function.
- **Variant Class**: SNVs overwhelmingly dominate the dataset.

---

## Expanded Workflow Explanation

### Makefile
The Makefile orchestrates the entire pipeline. Key targets include:
- **`genome_data`**: Downloads and prepares reference genome and GFF annotation.
- **`fastq`**: Fetches and converts SRA data to FASTQ format.
- **`alignments`**: Aligns reads to the reference genome.
- **`variants`**: Calls variants for each sample.
- **`merge`**: Merges all sample VCF files into a single file.
- **`vep`**: Annotates the merged VCF file using VEP.

### Tools Used
- **`bwa`**: For sequence alignment.
- **`samtools`**: For BAM file processing.
- **`bcftools`**: For variant calling and merging.
- **`vep`**: For variant annotation.

---

## Screenshots
The following screenshots illustrate the analysis results:
1. **Variants by Chromosome**:


2. **Coding Consequences**:


3. **Consequence Types**:
   

4. **Variant Classes**:
   

5. **Summary Statistics**:
   

---

## Discussion
The results highlight the utility of combining multiple samples to identify genetic variants and annotate their effects. While most variants are synonymous, the small subset of missense and inframe insertion variants are worth exploring for potential functional implications. Chromosomes V and XII show a higher density of variants, which could guide future studies focusing on these genomic regions.

---

## Recommendations for Future Work
1. **Deeper Annotation**: Include additional VEP plugins or databases for more detailed variant characterization.
2. **Functional Validation**: Prioritize missense variants and inframe insertions for experimental validation.
3. **Data Expansion**: Incorporate more SRA samples to increase variant detection sensitivity.

---

## Running the Pipeline
1. Clone the repository containing the `Makefile`, `design.csv`, and reference data.
2. Run the following command to execute the pipeline:
   ```bash
   make
