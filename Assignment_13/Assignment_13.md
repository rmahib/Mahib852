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

## Results
## Gene Expression Levels

Each row corresponds to a gene.  
Each column represents the read counts from a specific replicate:  
**SRR19383416_sorted.bam**, **SRR19383417_sorted.bam**, **SRR19383418_sorted.bam**  

| **Gene ID**       | **Replicate 1 (Counts)** | **Replicate 2 (Counts)** | **Replicate 3 (Counts)** |
|--------------------|--------------------------|--------------------------|--------------------------|
| **Gene 1 (5836)** | 796,955                  | 842,071                  | 921,303                  |
| **Gene 2 (1909)** | 569,941                  | 690,080                  | 647,602                  |
| **Gene 3 (6161)** | 442,479                  | 557,638                  | 513,660                  |
| **Gene 4 (2421)** | 415,645                  | 509,213                  | 473,398                  |
| **Gene 5 (2529)** | 277,252                  | 353,533                  | 305,020                  |

---

<img width="500" alt="13-12" src="https://github.com/user-attachments/assets/9081c1f1-b992-4bba-8c5c-b177cbe18c7d">


## Observations

## Gene 1 (ID: 5836):
- Highest expression across all replicates.
- Counts increase across replicates, suggesting a consistent trend.

## Gene 2 (ID: 1909):
- Second-highest expression.
- Counts are slightly more variable across replicates compared to Gene 1.

## Gene 3 to Gene 5 (IDs: 6161, 2421, 2529):
- Lower expression levels compared to Gene 1 and Gene 2.
- Counts remain consistent across replicates, suggesting stable expression.

---

## Insights
- Genes with consistently high counts (e.g., 5836 and 1909) could represent highly expressed or functionally important genes in your dataset.
- Variation between replicates appears minimal, indicating good experimental reproducibility.
- The counts follow a descending pattern from Gene 1 to Gene 5, reflecting their relative expression levels.

## 1. Count Matrix Summary
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

## 2. Visualization
IGV confirmed the alignment of reads to the reference genome. Below are the screenshots of RNA-Seq data visualized on IGV:

Chromosome VIII Region

Chromosome XI Region

Chromosome VII Region

<img width="842" alt="Assignment_13-9" src="https://github.com/user-attachments/assets/599e6987-968a-4c42-9dbf-097404038b5b">

<img width="851" alt="Assignment 13_10" src="https://github.com/user-attachments/assets/02f9785b-1188-4507-bd70-862b88370398">

<img width="850" alt="Assignment_13_11" src="https://github.com/user-attachments/assets/27702be2-6fb4-41f3-9264-02ebba9b4a00">


## 3. Consistently Expressed Genes

Consistently expressed genes across all samples were identified:
````
GeneID     SRR19383416 SRR19383417 SRR19383418
YDR438W    273         366         327
YDR523C     60          84          28
````


### Critical Observations

## Coverage: 

The IGV screenshots confirm robust coverage for exonic regions, validating the RNA-Seq dataset's quality.

## Gene Expression: 

The count matrix shows variability in gene expression across samples, with some genes exhibiting high expression levels, indicating active transcription.

## Consistency: 

Consistently expressed genes were observed across replicates, providing confidence in the experimental reproducibility.


## Discussion

This analysis demonstrates successful RNA-Seq data processing, from alignment to feature quantification and visualization. The results highlight gene expression patterns in Saccharomyces cerevisiae and provide insights into transcriptional activity.
