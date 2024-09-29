# Simulating FASTQ files

**Size of the file**

12M

**The total size of the genome**

12359733

**The number of chromosomes in the genome**

17

**The name (id) and length of each chromosome in the genome**
````
I dna:chromosome chromosome:R6 230218
II dna:chromosome chromosome:R 813184
III dna:chromosome chromosome: 316620
IV dna:chromosome chromosome:R 1531933
V dna:chromosome chromosome:R6 576874
VI dna:chromosome chromosome:R 270161
VII dna:chromosome chromosome: 1090940
VIII dna:chromosome chromosome 562643
IX dna:chromosome chromosome:R 439888
X dna:chromosome chromosome:R6 745751
XI dna:chromosome chromosome:R 666816
XII dna:chromosome chromosome: 1078177
XIII dna:chromosome chromosome 924431
XIV dna:chromosome chromosome: 784333
XV dna:chromosome chromosome:R 1091291
XVI dna:chromosome chromosome: 948066
Mito dna:chromosome chromosome 85779
````

**How many reads have you generated?**

240000

**What is the average read length?**

150 BP

How big are the FASTQ files?

79M

**After Compression**

16M (63M saved)

**Discuss whether you could get the same coverage with different parameter settings (read length vs. read number)**

Yes, I can achieve the same coverage with different combinations of read length and number of reads as long as the total number of bases sequenced (i.e., read length multiplied by number of reads) remains the same.

Example:
If I generate shorter reads, I'll need more reads to cover the same genome size.

If I generate longer reads, I'll need fewer reads to cover the same genome size

**Estimate Coverage for Other Genomes**


**Yeast (12 Mb genome)**

30x coverage = 12 Mb × 30 = 360 Mb total bases.

If we assume 150 bp reads, the number of reads required would be:
360,000,000 bp / 150 bp = 2.4 million reads

**Drosophila (180 Mb genome)**

30x coverage = 180 Mb × 30 = 5.4 Gb total bases.

Number of reads = 5,400,000,000 / 150 bp = 36 million reads.

**Human (3.2 Gb genome)**

30x coverage = 3.2 Gb × 30 = 96 Gb total bases.

Number of reads = 96,000,000,000 / 150 bp = 640 million reads

**Estimated file Size Uncompressed**

Yeast: 180MB

Drosophila: 2700 MB

Human: 48000 MB

**Script Link**

````
https://github.com/rmahib/Mahib852/blob/main/Assignment_5/alternative.sh

````

````
https://github.com/rmahib/Mahib852/blob/main/Assignment_5/assignment.sh
````