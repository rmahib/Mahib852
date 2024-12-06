# Variables
SAMPLES := $(shell awk -F',' 'NR>1 {print $$1}' design.csv)
GENOME_URL := ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
GFF_URL := ftp://ftp.ensembl.org/pub/release-105/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.105.gff3.gz
GENOME_FA := genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
GFF_FILE := genome_data/Saccharomyces_cerevisiae.R64-1-1.105_filtered.gff3.gz

# Targets
.PHONY: all clean

all: vep_results

# Download and preprocess genome
$(GENOME_FA): genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
	mkdir -p genome_data
	curl -o $@.gz $(GENOME_URL)
	gunzip -c $@.gz > $@
	samtools faidx $@

# Download and preprocess GFF file
$(GFF_FILE): genome_data/Saccharomyces_cerevisiae.R64-1-1.105.gff3.gz
	mkdir -p genome_data
	curl -o genome_data/Saccharomyces_cerevisiae.R64-1-1.105.gff3.gz $(GFF_URL)
	zcat genome_data/Saccharomyces_cerevisiae.R64-1-1.105.gff3.gz | \
	awk '!/Mito/' > genome_data/Saccharomyces_cerevisiae.R64-1-1.105_filtered.gff3
	bgzip -c genome_data/Saccharomyces_cerevisiae.R64-1-1.105_filtered.gff3 > $@
	tabix -p gff $@

# Download SRA data
fastq/%.fastq: design.csv
	mkdir -p fastq
	fasterq-dump $* --outdir fastq || { echo "Error downloading $*"; exit 1; }

# Align reads to reference genome
alignments/%.bam: fastq/%.fastq $(GENOME_FA)
	mkdir -p alignments
	bwa mem $(GENOME_FA) $< | samtools sort -o $@

# Call variants
variants/%.vcf.gz: alignments/%.bam $(GENOME_FA)
	mkdir -p variants
	bcftools mpileup -f $(GENOME_FA) $< | bcftools call -mv -Oz -o $@
	tabix -p vcf $@

# Merge VCF files
variants/merged_variants.vcf.gz: $(addprefix variants/,$(addsuffix .vcf.gz, $(SAMPLES)))
	bcftools merge $^ -Oz -o $@
	tabix -p vcf $@

# Annotate variants with VEP
vep_results: variants/merged_variants.vcf.gz $(GFF_FILE)
	mkdir -p results
	~/src/ensembl-vep/vep \
		-i $< -o results/vep_output.txt \
		--gff $(GFF_FILE) \
		--fasta $(GENOME_FA) \
		--force_overwrite

# Clean up
clean:
	rm -rf genome_data fastq alignments variants results
