# Variant Calling

## Objective
The objective of this assignment was to perform variant calling on SRA data aligned to a reference genome. We used the tools `bcftools` and `freebayes` for variant calling and validated the variants visually in IGV. This report includes:
- Description of the variant calling workflow.
- Answers to specific questions about alignment and variant data.
- Validation results using IGV screenshots.

## Workflow Summary

### Tools and Data
- **Reference Genome**: *Saccharomyces cerevisiae* (downloaded and indexed with BWA).
- **SRA Data**: Sample SRR5078495 from NCBI was chosen for alignment.
- **Alignment**: BWA was used to align reads to the reference genome.
- **Variant Calling**: 
  - `bcftools` and `freebayes` were used to call variants on the aligned BAM file.
  - VCF files were generated for both tools and analyzed for accuracy.
- **Visualization**: Variants were visualized in IGV to confirm the calls and check for false positives/negatives.

### Makefile
The Makefile automates the workflow with the following key targets:
- **genome**: Download and index the reference genome.
- **sra_data**: Download and convert SRA data to FASTQ.
- **align**: Align SRA reads to the reference genome.
- **call_bcftools**: Call variants using `bcftools`.
- **call_freebayes**: Call variants using `freebayes`.
- **compare**: Compare the variants identified by both tools.

## Analysis and Results

### Alignment Information
The alignment step produced the following:
- BAM file for the aligned reads.
- Coverage analysis to confirm that the sample had sufficient coverage for variant calling.

### Variant Calling Results
The variant calling with `bcftools` and `freebayes` produced the following VCF files:
- `srr5078495_bcftools.vcf`
- `srr5078495_freebayes.vcf`

These files were analyzed to identify differences between the two tools.

### Answers to Assignment Questions

1. **How many reads did not align with the reference genome?**
   - From the alignment stats, 10.99% of reads aligned, indicating that the remaining 89.01% of reads did not align with the reference.

2. **How many primary, secondary, and supplementary alignments are in the BAM file?**
   - Primary alignments: 2022016
   - Secondary alignments: 0
   - Supplementary alignments: 0

3. **How many properly paired alignments on the reverse strand are formed by reads contained in the first pair?**
   - Properly paired reverse strand alignments: 0 (no proper pairs found in the BAM file).

4. **Make a new BAM file that contains only the properly paired primary alignments with a mapping quality of over 10.**
   - A filtered BAM file was created using the following command:
     ```bash
     samtools view -h -b -f 2 -q 10 alignments/srr5078495_sorted.bam > alignments/srr5078495_filtered.bam
     ```
   - This filtered BAM file only includes reads with high confidence alignments.

5. **Compare the flagstats for the original and filtered BAM files.**
   - **Original BAM file**: 
     - Total reads: 18,402,980
     - Mapped reads: 10.99%
   - **Filtered BAM file**:
     - Total reads: 0 (after filtering for quality and pairing criteria).

### Variant Validation in IGV

The following images illustrate the results of variant validation in IGV:
1. **Chromosome V (Position 1-576,874)**

<img width="854" alt="Chromosome 5" src="https://github.com/user-attachments/assets/7b94422e-e534-4453-a55e-e2d215d436b1">


Observations: Variants identified by both `bcftools` and `freebayes` align with the reference, with some discrepancies visible.

2. **Chromosome X (Position 1-745,751)**

<img width="848" alt="Chromosome 10" src="https://github.com/user-attachments/assets/5f03e776-e7a7-4341-909b-52b86ee3dad0">


Observations: Both tools identified variants, but `freebayes` shows additional calls

3. **Chromosome VII (Position 1-1,090,940)**

<img width="851" alt="Chromosome 7" src="https://github.com/user-attachments/assets/703b334b-2da7-409f-8a3f-cb3518bd72eb">


Observations: Some variants identified by `freebayes` do not match the consensus from `bcftools`, suggesting potential false positives

4. **Chromosome IV (Position 1-1,531,933)**

<img width="852" alt="Chromosome 4" src="https://github.com/user-attachments/assets/86abd419-48fe-4ffb-8ee8-455390a99a9e">

Observations: Variants detected, but `freebayes` seems to have additional calls not corroborated by `bcftools`.

5. **Chromosome III (Position 1-316,620)**

<img width="843" alt="Chromosome 3" src="https://github.com/user-attachments/assets/4492365d-01a5-4d93-a06e-d957dbce4cd5">

Observations: Minimal variants were found, indicating a high concordance in this region

### Variant Caller Performance
The IGV analysis highlighted several instances of variant calling differences:
- **False Positives**: Freebayes occasionally detected variants not supported by bcftools.
- **False Negatives**: In regions with high coverage, certain true variants were missed by bcftools.
- **Concordant Variants**: Many variants were consistently identified by both tools, showing high reliability.

## Conclusion
This work involved a complete workflow from alignment to variant calling and validation. The comparison of `bcftools` and `freebayes` demonstrated differences in variant detection capabilities. Through IGV validation, we identified examples of both concordant and discordant variant calls. This analysis provides insights into the strengths and limitations of each variant caller, particularly in regions with high or low coverage.

## Some of the variants Analyzed

### BCFtools
- **Number of Variants**: SNPs identified on chromosomes `I` and `II`.
- **Example Variant**:
  - **Chromosome**: `I`
  - **Position**: `96809`
  - **Reference Allele**: `C`
  - **Alternate Allele**: `T`
  - **Genotype**: `1/1` (Homozygous alternate)
  - **Quality Score**: `67.41`
- **Additional Notes**:
  - Moderate sequencing depth (`DP` values).
  - Focused on high-confidence SNP detection.


### FreeBayes
- **Number of Variants**: SNPs and indels identified.
- **Example Variant**:
  - **Chromosome**: `I`
  - **Position**: `96809`
  - **Reference Allele**: `CTGT`
  - **Alternate Allele**: `TTGA` (Indel)
  - **Genotype**: `1/1` (Homozygous alternate)
  - **Quality Score**: `75.60`
- **Additional Notes**:
  - Detected more complex variants (e.g., indels).
  - Provides detailed allele frequency (`AF`), depth of coverage (`DP`), and alternate allele observations (`AO`).


## Comparison of Tools
- **Shared Findings**:
  - Both tools identified SNPs at overlapping positions (e.g., `96809` and `163683`).
- **Unique Findings**:
  - FreeBayes detected additional complex indels not reported by BCFtools.
- **Confidence**:
  - High confidence for SNPs from both tools.
  - Lower confidence in complex indels detected by FreeBayes.


## References
- Tools: `bcftools`, `freebayes`, `samtools`, IGV
- Reference Genome: *Saccharomyces cerevisiae* (Ensembl)
