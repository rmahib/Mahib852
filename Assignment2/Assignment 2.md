# Genome Feature Analysis of *Nannospalax galili* (Upper Galilee Blind Mole Rat)

**Organism Information**
*Nannospalax galili*
** Tell us bit about the organism**  
*Nannospalax galili*, also known as the Upper Galilee blind mole rat, belongs to the Spalacidae family. This species is well-known for its ability to survive in low oxygen environments. It is blind and uses touch and vibrations to navigate.

**Total Number of Features**

The total number of features in this GFF3 file is 988,215.

**Number of Sequence Regions**

The organism contains 154,975 sequence regions

**Total Number of Genes**

There are 57,497 genes listed for the organism

**Feature Breakdown**

Below is a breakdown of the main annotated features for Nannospalax galili:

-298,499 exons
-282,550 CDS (coding sequences)
-162,234 regions
-26,872 mRNA
-17,719 five_prime_UTR
-14,741 three_prime_UTR
-5,623 ncRNA_gene


**Observation**

The large number of Exons and CDS annotations suggests that the coding region of the genome is well characterized. However, it has also have a number of ncRNA which indicates that it might be under annotated. It also does have some general biological regions. Considering everything it could be said the organism is well annotated but might under annotated in some regions.

**Pseudogenes**

The organism contains 649 pseudogenes

**Command Used for downloading the file**

````
wget https://ftp.ensembl.org/pub/current_gff3/nannospalax_galili/Nannospalax_galili.S.galili_v1.0.112.gff3.gz
````
**Question 1**

````
zgrep -v "^#" Nannospalax_galili.S.galili_v1.0.112.gff3.gz | wc -l
````
**Question 2**

````zgrep "##sequence-region" Nannospalax_galili.S.galili_v1.0.112.gff3.gz | wc -l
````
**Question 3**

````
zgrep -w "gene" Nannospalax_galili.S.galili_v1.0.112.gff3.gz | wc -l
````
**Question4**


````
zcat Nannospalax_galili.S.galili_v1.0.112.gff3.gz | awk '$3 != "" {print $3}' | sort | uniq -c | sort -nr | head -10
````
**Psudogene**

````
zgrep -w "pseudogene" Nannospalax_galili.S.galili_v1.0.112.gff3.gz" | wc -l
````


