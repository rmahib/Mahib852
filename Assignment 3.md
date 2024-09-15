** Visualizing GFF file of my choice**

First I strated with the following approach. I decide to go with *Mycoplasma_genitalium*

**Command Used for Downloading the file**

`mkdir Mycoplasma_genitalium`

`cd Mycoplasma_genitalium`

** Command used for downloading the corresponding gff and fna**

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz`

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz`

`gunzip GCF_000027325.1_ASM2732v1_genomic.fna.gz`

`gunzip GCF_000027325.1_ASM2732v1_genomic.gff.gz`

#Corresponding gff and fna file were downloaded in directory. The files were loaded in IGV and relevant screenshots were taken

**Command used for getting gene.gff and cds.gff**


`cat GCF_000027325.1_ASM2732v1_genomic.gff | awk '$3 == "gene"' > mycoplasma_genes.gff`


`cat GCF_000027325.1_ASM2732v1_genomic.gff | awk '$3 == "CDS"' > mycoplasma_CDS.gff`


`head mycoplasma_genes.gff`

`head mycoplasma_CDS.gff`


**Itergenic Region**

`awk '{print $1"\t.\tintergenic_region\t"$2"\t"$3"\t.\t.\t.\tID=intergenic_"NR";Name=intergenic_region_"NR}' gene_intervals.bed > intergenic_regions.gff`

`cat intergenic_regions.gff`

NC_000908.2     .       intergenic_region       0       686     .       .       .       ID=intergenic_1;Name=intergenic_region_1

NC_000908.2     .       intergenic_region       2760    2845    .       .       .       ID=intergenic_2;Name=intergenic_region_2

NC_000908.2     .       intergenic_region       4797    4812    .       .       .       ID=intergenic_3;Name=intergenic_region_3

NC_000908.2     .       intergenic_region       8547    8551    .       .       .       ID=intergenic_4;Name=intergenic_region_4

NC_000908.2     .       intergenic_region       9920    9923    .       .       .       ID=intergenic_5;Name=intergenic_region_5

NC_000908.2     .       intergenic_region       12039   12068   .       .       .       ID=intergenic_6;Name=intergenic_region_6


**Genes protein_coding**

`awk '$3 == "gene" && $9 ~ /biotype=protein_coding/ {
    split($9, a, ";");
    id="";
    name="";
    for (i in a) {
        if (a[i] ~ /^ID=/) id=a[i];
        if (a[i] ~ /^Name=/) name=a[i];
    }
    if (name) {
        print $1": "id"; "name"; biotype=protein_coding";
    } else {
        print $1": "id"; biotype=protein_coding";
    }
}' mycoplasma_genes.gff3`

NC_000908.2: ID=gene-MG_RS00005; Name=dnaN; biotype=protein_coding

NC_000908.2: ID=gene-MG_RS00010; Name=MG_RS00010; biotype=protein_coding

NC_000908.2: ID=gene-MG_RS00015; Name=gyrB; biotype=protein_coding

NC_000908.2: ID=gene-MG_RS00020; Name=gyrA; biotype=protein_coding

NC_000908.2: ID=gene-MG_RS00025; Name=serS; biotype=protein_coding

NC_000908.2: ID=gene-MG_RS00030; Name=tmk; biotype=protein_coding

NC_000908.2: ID=gene-MG_RS00035; Name=MG_RS00035; biotype=protein_coding

NC_000908.2: ID=gene-MG_RS00040; Name=mnmE; biotype=protein_coding

NC_000908.2: ID=gene-MG_RS00045; Name=MG_RS00045; biotype=protein_coding

NC_000908.2: ID=gene-MG_RS00050; Name=MG_RS00050; biotype=protein_coding


**Number of genes**

--439

`awk '$3 == "gene"' GCF_000027325.1_ASM2732v1_genomic.gff | wc -l`

**Number of CDS**

--10

`awk '$3 == "CDS"' GCF_000027325.1_ASM2732v1_genomic.gff | wc -l`

**Manual gff**

A manual gff was created and uploaded as a seperate track in IGV

**Gene only gff**

`awk '$3 == "gene"' mycoplasma_genes.gff3 | head -n 20`

NC_000908.2     RefSeq  gene    686     1828    .       +       .       ID=gene-MG_RS00005;Dbxref=GeneID:88282116;Name=dnaN;gbkey=Gene;gene=dnaN;gene_biotype=protein_coding;locus_tag=MG_RS00005;old_locus_tag=MG_001

NC_000908.2     RefSeq  gene    1828    2760    .       +       .       ID=gene-MG_RS00010;Dbxref=GeneID:88282117;Name=MG_RS00010;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MG_RS00010;old_locus_tag=MG_002

NC_000908.2     RefSeq  gene    2845    4797    .       +       .       ID=gene-MG_RS00015;Dbxref=GeneID:88282118;Name=gyrB;gbkey=Gene;gene=gyrB;gene_biotype=protein_coding;locus_tag=MG_RS00015;old_locus_tag=MG_003

NC_000908.2     RefSeq  gene    4812    7322    .       +       .       ID=gene-MG_RS00020;Dbxref=GeneID:88282119;Name=gyrA;gbkey=Gene;gene=gyrA;gene_biotype=protein_coding;locus_tag=MG_RS00020;old_locus_tag=MG_004

NC_000908.2     RefSeq  gene    7294    8547    .       +       .       ID=gene-MG_RS00025;Dbxref=GeneID:88282120;Name=serS;gbkey=Gene;gene=serS;gene_biotype=protein_coding;locus_tag=MG_RS00025;old_locus_tag=MG_005

NC_000908.2     RefSeq  gene    8551    9183    .       +       .       ID=gene-MG_RS00030;Dbxref=GeneID:88282121;Name=tmk;gbkey=Gene;gene=tmk;gene_biotype=protein_coding;locus_tag=MG_RS00030;old_locus_tag=MG_006

NC_000908.2     RefSeq  gene    9156    9920    .       +       .       ID=gene-MG_RS00035;Dbxref=GeneID:88282122;Name=MG_RS00035;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MG_RS00035;old_locus_tag=MG_007

NC_000908.2     RefSeq  gene    9923    11251   .       +       .       ID=gene-MG_RS00040;Dbxref=GeneID:88282123;Name=mnmE;gbkey=Gene;gene=mnmE;gene_biotype=protein_coding;locus_tag=MG_RS00040;old_locus_tag=MG_008

NC_000908.2     RefSeq  gene    11251   12039   .       +       .       ID=gene-MG_RS00045;Dbxref=GeneID:88282124;Name=MG_RS00045;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MG_RS00045;old_locus_tag=MG_009

NC_000908.2     RefSeq  gene    12068   12724   .       +       .       ID=gene-MG_RS00050;Dbxref=GeneID:88282125;Name=MG_RS00050;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MG_RS00050;old_locus_tag=MG_010



##The total number of genes were way more than the number of CDS.

## In total genes were identified with gene gene_intervals

#Upon Visualizing in IGV genes were found to be merged with the sequence. Hence starting with a start codon and ending with a stop codon

#I've loaded the gff made manually but could not able to align it with the sequence. I could not find what went wrong

**Importance of this report**

#It tells us Where the gene is located on the chromosome

#Whether the gene is on the forward (+) or reverse (-) strand

#By having the genomic coordinates and biotype (protein-coding), researchers can study the structural organization of the genome, track gene expression, and understand how these genes function within the organism’s biology


By using IGV I could see the gene was overlapping with the cds region