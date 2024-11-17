# Variant Effect Prediction using VEP

## Overview

This report summarizes the process and results of predicting the effects of genetic variants in **Saccharomyces cerevisiae** using the Variant Effect Predictor (VEP). The workflow incorporated generating a VCF file through variant calling and performing downstream analysis with VEP to identify variant impacts on genes, transcripts, and proteins.

---

## Methods

### 1. Data Preparation
- **Reference Genome:** The reference genome for *Saccharomyces cerevisiae* was downloaded and indexed using BWA.
- **SRA Reads:** Single-end reads were obtained from the SRA database using `SRR5078495`.
- **Alignment:** Reads were aligned to the reference genome using BWA, followed by sorting and indexing with `samtools`.
- **Variant Calling:** Variants were called with `bcftools mpileup` and filtered to produce a compressed VCF file.

### 2. VEP Setup
- **GFF3 File Preparation:** The annotation file was sorted, compressed, and indexed using `bgzip` and `tabix`.
- **VEP Execution:** VEP was run using the prepared GFF3 and reference genome FASTA file to annotate the variants.

---

## Results

### 1. Variant Categories
The following types of variants were identified:
- **Synonymous Variants:** Changes in the nucleotide sequence that do not alter the amino acid sequence of the encoded protein.
- **Missense Variants:** Variants causing single amino acid changes in the protein sequence.
- **Frameshift Variants:** Insertions or deletions that disrupt the reading frame of the gene.
- **Start/Stop Codon Variants:** Variants that affect the start or stop codon, potentially truncating or elongating the protein.
- **Upstream/Downstream Variants:** Variants located upstream or downstream of genes that may affect gene regulation.

### 2. Detailed Analysis
| **Variant ID**   | **Location** | **Consequence**          | **Gene**    | **Impact**   | **Distance** | **Codon Change** |
|-------------------|--------------|--------------------------|-------------|--------------|--------------|-------------------|
| `I_96809_C/T`     | I:96809      | Synonymous Variant       | `YAL026C`   | Low          | -            | `caG/caA`         |
| `IV_654214_A/C`   | IV:654214    | Missense Variant         | `YDR099W`   | Moderate     | -            | `gAt/gCt`         |
| `IV_116422_G/T`   | IV:116422    | Missense Variant         | `YDL192W`   | Moderate     | -            | `ttG/ttT`         |
| `II_163695_T/C`   | II:163695    | Synonymous Variant       | `YBL030C`   | Low          | -            | `caA/caG`         |
| `VII_451122_G/-`  | VII:451122   | Frameshift Variant       | `YHR012W`   | High         | -            | `g/-`             |
| `IX_920877_AT/A`  | IX:920877    | Frameshift Deletion      | `YIL054C`   | High         | -            | `atg/a`           |
| `XIV_1267313_G/A` | XIV:1267313  | Stop Gained              | `YNL009W`   | High         | -            | `gaT/gaA`         |
| `XII_31988_C/T`   | XII:31988    | Start Lost               | `YLR097C`   | High         | -            | `aTg/aCg`         |

---

## Conclusion
The VEP analysis identified a range of variant effects across coding and regulatory regions of the genome in *Saccharomyces cerevisiae*. The majority of the variants were synonymous or missense, while some frameshift, stop-gained, and start-lost variants were also detected, suggesting potentially significant impacts on protein structure and function. These results provide a foundation for further experimental validation to explore their biological relevance.

---
