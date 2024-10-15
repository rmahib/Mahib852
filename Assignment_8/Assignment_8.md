## Generate a BAM alignment file

````
# Alignment of Simulated and SRA Reads Using BWA and Visualization in IGV

## Introduction
In this analysis, we aligned both simulated reads and real sequencing data from the SRA (Sequence Read Archive) to the **Saccharomyces cerevisiae** reference genome using the BWA aligner. The goal was to compare the quality of alignments between simulated and real sequencing data. The alignments were visualized using **IGV (Integrative Genomics Viewer)**, and alignment statistics were generated using **Samtools**.

Additionally, we created a **Makefile** to streamline the analysis by automating tasks such as downloading the genome, simulating reads, trimming reads, indexing the reference genome, aligning the reads, and generating sorted BAM files.

## Makefile Overview
The **Makefile** was designed to automate the steps of aligning reads to the reference genome and managing the preprocessing steps. Two key targets, `index` and `align`, were added to the Makefile to handle the indexing of the reference genome and the alignment of both simulated and real reads.

### Key Targets in the Makefile

- **genome**: Downloads and extracts the reference genome.
- **simulate**: Simulates reads from the reference genome.
- **download**: Downloads the real sequencing data from SRA.
- **trim**: Trims the SRA reads using Trimmomatic.
- **index**: Indexes the reference genome using BWA.
- **align**: Aligns both simulated and SRA reads to the indexed reference genome and produces sorted BAM files.

### Makefile Code
Here is the core section of the Makefile that includes the `index` and `align` targets:

```makefile
# Index the reference genome
index: genome
	@echo "Indexing the reference genome with BWA..."
	$(BWA) index genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

# Align the reads to the reference genome and sort the BAM files
align: index trim simulate
	mkdir -p alignments
	@echo "Aligning simulated reads to the reference genome..."
	$(BWA) mem genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa simulated_reads_1.fastq simulated_reads_2.fastq | $(SAMTOOLS) view -Sb - > alignments/simulated_reads.bam
	$(SAMTOOLS) sort alignments/simulated_reads.bam -o alignments/simulated_reads_sorted.bam
	$(SAMTOOLS) index alignments/simulated_reads_sorted.bam
	@echo "Aligning SRA reads to the reference genome..."
	$(BWA) mem genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa sra_data/SRR387901_1.fastq | $(SAMTOOLS) view -Sb - > alignments/sra_reads.bam
	$(SAMTOOLS) sort alignments/sra_reads.bam -o alignments/sra_reads_sorted.bam
	$(SAMTOOLS) index alignments/sra_reads_sorted.bam
````
**Explanation of Targets**
1. index: This target uses BWA to index the reference genome. Indexing is a necessary step before aligning any reads.
2. align: This target aligns both the simulated reads and the SRA reads to the indexed reference genome using BWA, and the resulting BAM files are sorted and indexed using Samtools

**Materials and Methods**

1.Reference Genome Download

2.Simulating Reads

3.Downloading SRA Data

4.Trimming the SRA Reads

5.Aligning the Reads to the Reference Genome

6.Visualization in IGV

7.Alignment Statistics

**Code Used**

````
wget ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
````
````
wgsim -N 240000 -1 150 -2 150 Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa simulated_reads_1.fastq simulated_reads_2.fastq
````
````
fastq-dump --split-files --outdir sra_data SRR387901
````
````
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 sra_data/SRR387901_1.fastq sra_trimmed.fastq \
  ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
````
````
# Simulated reads
bwa mem Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa simulated_reads_1.fastq simulated_reads_2.fastq | samtools view -Sb - > simulated_reads.bam
samtools sort simulated_reads.bam -o simulated_reads_sorted.bam
samtools index simulated_reads_sorted.bam

# SRA reads
bwa mem Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa sra_trimmed.fastq | samtools view -Sb - > sra_reads.bam
samtools sort sra_reads.bam -o sra_reads_sorted.bam
samtools index sra_reads_sorted.bam
````
````
samtools flagstat simulated_reads_sorted.bam
````

````
480000 + 0 in total (QC-passed reads + QC-failed reads)
480000 + 0 mapped (100.00% : N/A)
480000 + 0 properly paired (100.00% : N/A)
````
````
samtools flagstat sra_reads_sorted.bam
````
````
18402980 + 0 in total (QC-passed reads + QC-failed reads)
2022016 + 0 mapped (10.99% : N/A)
0 + 0 properly paired (N/A : N/A)
````
## Results
**Visualization in IGV**

The BAM files were visualized in IGV:

**Simulated Reads** 

The alignments were perfect, with 100% of the reads mapping to the reference genome. The coverage was uniform across the entire genome, as expected from simulated data.

**SRA Reads** 

The real sequencing data showed much lower alignment rates (~10.99% mapped). This could indicate that the reads either had low quality or came from a strain that differs from the reference genome. There were also significant gaps in the coverage.
Differences Between Simulated and SRA Datasets

**Simulated Reads** 

These reads aligned perfectly to the reference genome (100% alignment), which is expected since the reads were generated directly from the reference genome sequence.

**SRA Reads** 

Only 10.99% of the SRA reads mapped to the reference genome. This low mapping rate suggests that the SRA reads may come from a different strain of Saccharomyces cerevisiae or could be of lower quality compared to the simulated reads. Additionally, there were no properly paired reads, which could indicate issues with sequencing or read quality.

**Conclusion**

The alignment of simulated reads served as a baseline for perfect mapping, as expected from simulated data. However, the SRA dataset exhibited much lower alignment efficiency, which highlights the challenges of working with real-world sequencing data, such as quality issues or genome differences. The Makefile streamlined the process of indexing, aligning, and generating BAM

**Makefile link**

https://github.com/rmahib/Mahib852/blob/main/Assignment_8/makefile