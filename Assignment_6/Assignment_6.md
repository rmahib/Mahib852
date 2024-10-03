## Quality Control Analysis

The data analyzed in this report comes from the Sequence Read Archive (SRA) under the accession number **SRR387901**. This dataset corresponds to a single-end sequencing run of *Saccharomyces cerevisiae* (yeast), generated as part of the study revealing pervasive translational control in meiosis and helps to illuminate the molecular basis of the broad restructuring of meiotic cells.

The study, titled **"High-resolution view of the yeast meiotic program revealed by ribosome profiling"**, was conducted by **Brar GA, Yassour M, Friedman N, Regev A, Ingolia NT, Weissman JS** and published in **Science**.

## Results of the Analysis

### 1. Initial Data Quality (Pre-Trimming)
After downloading the data, an initial quality control check using FastQC revealed the following issues:
- **Low per-base sequence quality**: The quality scores dropped significantly toward the end of the reads.
- **Adapter contamination**: Adapter sequences were detected in a large number of reads.
- **Overrepresented sequences**: A significant proportion of reads were overrepresented.

### 2. Improvement after Trimming
The dataset was trimmed using Trimmomatic to remove low-quality bases and adapter sequences. After trimming, the quality improved:
- **Per-base sequence quality** improved, with fewer low-quality bases.
- **Adapter sequences** were successfully removed.
- **Overrepresented sequences** were reduced to acceptable levels.

A detailed comparison of the quality metrics before and after trimming is provided below:

| Metric                          | Before Trimming | After Initial Trimming | After Aggressive Trimming |
|----------------------------------|-----------------|------------------------|---------------------------|
| Per Base Sequence Quality        | Poor            | Improved               | Good                      |
| Adapter Content                  | Present         | Removed                | Removed                   |
| Overrepresented Sequences        | Present         | Reduced                | None                      |

**Script Link**

````
https://github.com/rmahib/Mahib852/blob/main/Assignment_6/Script.sh
````