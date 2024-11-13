# Variables
GENOME_URL := ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
SRA_ID := SRR5078495
GENOME_FILE := genome_data/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
BAM_FILE := alignments/srr5078495_sorted.bam
BCFTOOLS_VCF := variants/srr5078495_bcftools.vcf
FREEBAYES_VCF := variants/srr5078495_freebayes.vcf

# Tools
BWA := bwa
SAMTOOLS := samtools
FASTQ_DUMP := fastq-dump
BCFTOOLS := bcftools
FREEBAYES := freebayes

# Targets
.PHONY: usage genome sra_data index align call_bcftools call_freebayes compare clean

# Help
usage:
	@echo "Makefile targets:"
	@echo "  genome         - Download and unzip the reference genome"
	@echo "  sra_data       - Download SRA data for analysis"
	@echo "  index          - Index the reference genome"
	@echo "  align          - Align SRA reads to the reference genome (single-end)"
	@echo "  call_bcftools  - Call variants using bcftools"
	@echo "  call_freebayes - Call variants using freebayes"
	@echo "  compare        - Compare results between bcftools and freebayes"
	@echo "  clean          - Remove all generated files to start fresh"

# Download and prepare the genome
genome:
	mkdir -p genome_data
	@echo "Downloading reference genome..."
	wget -P genome_data $(GENOME_URL)
	@echo "Unzipping reference genome..."
	gunzip -f genome_data/*.fa.gz

# Download the SRA data
sra_data:
	mkdir -p sra_data
	@echo "Downloading SRA data..."
	$(FASTQ_DUMP) --split-files --outdir sra_data $(SRA_ID)

# Index the reference genome
index: genome
	@echo "Indexing the reference genome with BWA..."
	$(BWA) index $(GENOME_FILE)

# Align the SRA reads to the reference genome and sort the BAM files (single-end)
align: index sra_data
	mkdir -p alignments
	@echo "Aligning SRA reads to the reference genome..."
	$(BWA) mem $(GENOME_FILE) sra_data/$(SRA_ID)_1.fastq | $(SAMTOOLS) view -Sb - > alignments/srr5078495.bam
	$(SAMTOOLS) sort alignments/srr5078495.bam -o $(BAM_FILE)
	$(SAMTOOLS) index $(BAM_FILE)

# Call variants using bcftools
call_bcftools: align
	mkdir -p variants
	@echo "Calling variants with bcftools..."
	$(BCFTOOLS) mpileup -f $(GENOME_FILE) $(BAM_FILE) | $(BCFTOOLS) call -mv -Ov -o $(BCFTOOLS_VCF)

# Call variants using freebayes
call_freebayes: align
	mkdir -p variants
	@echo "Calling variants with freebayes..."
	$(FREEBAYES) -f $(GENOME_FILE) $(BAM_FILE) > $(FREEBAYES_VCF)

# Compare the VCF files
compare: call_bcftools call_freebayes
	@echo "Comparing the VCF files from bcftools and freebayes..."
	diff $(BCFTOOLS_VCF) $(FREEBAYES_VCF) || true

# Clean up all generated files
clean:
	rm -rf genome_data sra_data alignments variants




