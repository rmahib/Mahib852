# BAM File Analysis and Filtering Report

## Introduction
In this analysis, we aligned sequencing data from the SRA (Sequence Read Archive) to the **Saccharomyces cerevisiae** reference genome and performed various filtering and statistical analyses on the resulting BAM files.

## Commands and Answers

### 1. How many reads did not align with the reference genome?
To find the number of reads that did not align, we used the `samtools flagstat` command on the original BAM file. This command provides the total reads and the number of mapped reads, allowing us to calculate the unmapped reads.

**Command**:

````
samtools flagstat alignments/sra_reads_sorted.bam
````

**Output**
````
18402980 + 0 in total (QC-passed reads + QC-failed reads)
2022016 + 0 mapped (10.99% : N/A)
````
There were 16,180,964 reads that did not align with the reference genome (calculated as 18,402,980 total reads minus 2,022,016 mapped reads).

**How many primary, secondary, and supplementary alignments are in the BAM file?**

The samtools flagstat command also provides information on primary, secondary, and supplementary alignments

**Command**

````
samtools flagstat alignments/sra_reads_sorted.bam
````

**Output**

````
18402980 + 0 primary
0 + 0 secondary
0 + 0 supplementary
````
**Answer:**

Primary alignments: 18,402,980
Secondary alignments: 0
Supplementary alignments: 0

**How many properly-paired alignments on the reverse strand are formed by reads contained in the first pair?**

**Command**

````
samtools view -f 1 -F 4 -F 8 -F 256 -F 2048 alignments/sra_reads_sorted.bam | awk '$2 ~ /16/' | wc -l
````

**Output**
````
0
````
**Answer:** 
There were 0 properly paired alignments on the reverse strand formed by reads in the first pair.

**Make a new BAM file that contains only the properly paired primary alignments with a mapping quality of over 10**

**Command**

````
samtools view -h -f 2 -q 10 -b alignments/sra_reads_sorted.bam > alignments/sra_reads_filtered.bam
samtools index alignments/sra_reads_filtered.bam
````
This command created the file alignments/sra_reads_filtered.bam containing the filtered alignments.

**Compare the flagstats for your original and filtered BAM files**

**Command**

````
samtools flagstat alignments/sra_reads_sorted.bam
````
**Output**

````
18402980 + 0 in total
2022016 + 0 mapped (10.99% : N/A)
````
**Command for Filtered BAM:**

````
samtools flagstat alignments/sra_reads_filtered.bam
````
**Output**

````
0 + 0 in total
0 + 0 mapped
````

**Answer:**

Original BAM: 18,402,980 total reads with 2,022,016 mapped reads (10.99% mapped).
Filtered BAM: 0 reads remained after filtering, indicating that no reads met the criteria for proper pairing and mapping quality >10.

**Conclusion**

This analysis examined alignment statistics and applied quality and pairing filters to the BAM file. The results revealed that most reads in the original BAM file did not align or failed to meet the quality criteria after filtering.

**Makefile Link**

https://github.com/rmahib/Mahib852/blob/main/Assignment_9/makefile
