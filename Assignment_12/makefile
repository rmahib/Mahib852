# Makefile for variant calling and annotation pipeline

SAMPLES := $(shell awk -F',' 'NR>1 {print $$1}' design.csv)
GENOME_URL := ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
GFF_URL := ftp://ftp.ensembl.org/pub/release-105/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.105.gff3.gz

.PHONY: all clean

all: vep_results

# Download reference genome
genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz:
	mkdir -p genome_data
	curl -o $@ $(GENOME_URL)

# Download GFF file and index it
genome_data/Saccharomyces_cerevisiae.R64-1-1.105_filtered.gff3.gz:
	mkdir -p genome_data
	curl -o $@ $(GFF_URL)
	tabix -p gff $@

# Decompress genome and index
genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa: genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
	gunzip -c $< > $@
	samtools faidx $@

# Download SRA data and convert to fastq
fastq/%.fastq: design.csv
	mkdir -p fastq
	fasterq-dump $*

# Align reads to reference genome
alignments/%.bam: fastq/%.fastq genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
	mkdir -p alignments
	bwa mem genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa $< | samtools sort -o $@

# Variant calling
variants/%.vcf.gz: alignments/%.bam
	mkdir -p variants
	bcftools mpileup -f genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa $< | \
	bcftools call -mv -Oz -o $@

# Merge all VCF files
variants/merged_variants.vcf.gz: $(addprefix variants/,$(addsuffix .vcf.gz, $(SAMPLES)))
	bcftools merge $^ -Oz -o $@
	tabix -p vcf $@

# Annotate variants using VEP
vep_results: variants/merged_variants.vcf.gz genome_data/Saccharomyces_cerevisiae.R64-1-1.105_filtered.gff3.gz
	mkdir -p results
	vep -i $< -o results/vep_output.txt --gff genome_data/Saccharomyces_cerevisiae.R64-1-1.105_filtered.gff3.gz \
		--fasta genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa --force_overwrite

# Clean up intermediate files
clean:
	rm -rf genome_data fastq alignments variants results


