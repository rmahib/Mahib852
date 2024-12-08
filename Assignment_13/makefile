# Variables
GENOME_URL = ftp://ftp.ensembl.org/pub/release-105/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
GTF_URL = ftp://ftp.ensembl.org/pub/release-105/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.105.gtf.gz
SAMPLES_CSV = sample_info.csv
FASTQ_DIR = fastq
ALIGNMENTS_DIR = alignments/sorted
COUNTS_DIR = counts
GENOME_DIR = genome
INDEX_PREFIX = $(GENOME_DIR)/scerevisiae_index
GENOME_FA = $(GENOME_DIR)/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
GTF_FILE = $(GENOME_DIR)/Saccharomyces_cerevisiae.R64-1-1.105.gtf
R_SCRIPT = analyze_counts.R
R_OUTPUT_DIR = r_analysis

# Programs
FASTQ_DUMP = fasterq-dump
HISAT2 = hisat2
HISAT2_BUILD = hisat2-build
SAMTOOLS = samtools
FEATURECOUNTS = featureCounts
RSCRIPT = Rscript

# Default Target
.PHONY: all
all: setup design_csv download_genome build_index align count_matrix analyze

# Step 1: Set up directories
.PHONY: setup
setup:
	mkdir -p $(FASTQ_DIR) $(ALIGNMENTS_DIR) $(COUNTS_DIR) $(GENOME_DIR) $(R_OUTPUT_DIR)

# Step 2: Generate the design file if it doesn't exist
.PHONY: design_csv
design_csv: $(SAMPLES_CSV)

$(SAMPLES_CSV):
	@echo "Sample_Name,SRA_Accession,Condition,Replicate" > $@
	@echo "T0_rep1,SRR19383418,T0,1" >> $@
	@echo "T0_rep2,SRR19383417,T0,2" >> $@
	@echo "T0_rep3,SRR19383416,T0,3" >> $@

# Step 3: Download genome and GTF files
.PHONY: download_genome
download_genome: $(GENOME_FA) $(GTF_FILE)

$(GENOME_FA):
	wget -O $@ $(GENOME_URL)

$(GTF_FILE):
	wget -O $@ $(GTF_URL)

# Step 4: Build HISAT2 index
.PHONY: build_index
build_index: $(INDEX_PREFIX).1.ht2

$(INDEX_PREFIX).1.ht2: $(GENOME_FA)
	$(HISAT2_BUILD) $< $(INDEX_PREFIX)

# Step 5: Download FASTQ files
.PHONY: download_fastq
download_fastq: setup design_csv
	cat $(SAMPLES_CSV) | tail -n +2 | cut -d ',' -f2 | while read SRA; do \
		$(FASTQ_DUMP) --split-files -O $(FASTQ_DIR) $$SRA; \
	done

# Step 6: Align reads
.PHONY: align
align: $(SAMPLES_CSV) build_index
	cat $(SAMPLES_CSV) | tail -n +2 | cut -d ',' -f2 | while read SRA; do \
		$(HISAT2) -x $(INDEX_PREFIX) \
		-1 $(FASTQ_DIR)/$$SRA\_1.fastq \
		-2 $(FASTQ_DIR)/$$SRA\_2.fastq \
		-S $(ALIGNMENTS_DIR)/$$SRA.sam; \
		$(SAMTOOLS) view -S -b $(ALIGNMENTS_DIR)/$$SRA.sam > $(ALIGNMENTS_DIR)/$$SRA.bam; \
		$(SAMTOOLS) sort $(ALIGNMENTS_DIR)/$$SRA.bam -o $(ALIGNMENTS_DIR)/$$SRA\_sorted.bam; \
		$(SAMTOOLS) index $(ALIGNMENTS_DIR)/$$SRA\_sorted.bam; \
		rm $(ALIGNMENTS_DIR)/$$SRA.sam $(ALIGNMENTS_DIR)/$$SRA.bam; \
	done

# Step 7: Generate the count matrix
.PHONY: count_matrix
count_matrix: $(COUNTS_DIR)/gene_counts.txt

$(COUNTS_DIR)/gene_counts.txt: $(GTF_FILE)
	mkdir -p $(COUNTS_DIR)
	$(FEATURECOUNTS) -a $< -o $@ $(ALIGNMENTS_DIR)/*_sorted.bam

# Step 8: Perform R Analysis
.PHONY: analyze
analyze: $(COUNTS_DIR)/gene_counts.txt
	$(RSCRIPT) $(R_SCRIPT) $(COUNTS_DIR)/gene_counts.txt $(R_OUTPUT_DIR)

# Cleanup intermediate files
.PHONY: clean
clean:
	rm -rf $(FASTQ_DIR) $(ALIGNMENTS_DIR) $(COUNTS_DIR) $(GENOME_DIR)
	rm -f $(SAMPLES_CSV)
	rm -rf $(R_OUTPUT_DIR)
