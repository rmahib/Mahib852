# Makefile for VEP analysis pipeline

# Variables
REF=genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
GFF=genome_data/Saccharomyces_cerevisiae.R64-1-1.105.gff3
SRA=SRR5078495
VCF=variants/sample_variants.vcf

all: usage download_genome sra_data index align call_variants vep

usage:
	@echo "Makefile targets:"
	@echo "  usage          - Show this help message"
	@echo "  download_genome - Download and unzip the reference genome and GFF"
	@echo "  sra_data       - Download SRA data and convert to FASTQ"
	@echo "  index          - Index the reference genome with BWA"
	@echo "  align          - Align SRA reads to the reference genome and sort BAM files"
	@echo "  call_variants  - Call variants using bcftools"
	@echo "  vep            - Predict variant effects using VEP"
	@echo "  clean          - Remove all downloaded and generated files"

download_genome:
	mkdir -p genome_data
	wget -P genome_data ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
	wget -P genome_data ftp://ftp.ensembl.org/pub/release-105/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.105.gff3.gz
	gunzip -f genome_data/*.gz

sra_data:
	mkdir -p sra_data
	fastq-dump --split-files --outdir sra_data $(SRA)

index: download_genome
	bwa index $(REF)

align: sra_data index
	mkdir -p alignments
	bwa mem $(REF) sra_data/$(SRA)_1.fastq | samtools view -Sb - | samtools sort -o alignments/$(SRA)_sorted.bam
	samtools index alignments/$(SRA)_sorted.bam

call_variants: align
	mkdir -p variants
	bcftools mpileup -Ou -f $(REF) alignments/$(SRA)_sorted.bam | bcftools call -mv -Oz -o $(VCF).gz
	bcftools index $(VCF).gz
	bcftools view $(VCF).gz > $(VCF)

vep: call_variants
	mkdir -p results
	bgzip -c $(GFF) > $(GFF).gz
	tabix -p gff $(GFF).gz
	micromamba run -n vep \
		~/src/ensembl-vep/vep \
		-i $(VCF) \
		-o results/vep_output.txt \
		--gff $(GFF).gz \
		--fasta $(REF) \
		--force_overwrite

clean:
	rm -rf genome_data sra_data alignments variants results
