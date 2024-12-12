# Makefile for RNA-Seq Differential Expression Analysis

# Variables
SAMPLES_CSV = sample_info_simulated.csv
COUNTS_MATRIX = counts/simulated_gene_counts.txt
RESULTS_DIR = results
R_SCRIPT = deseq2_analysis.R
RSCRIPT = Rscript

.PHONY: all clean setup analyze

# Default target
all: setup analyze

# Step 1: Setup directories and files
setup:
	mkdir -p counts $(RESULTS_DIR)
	cp $(SAMPLES_CSV) design.csv

# Step 2: Differential expression analysis using DESeq2
analyze: $(COUNTS_MATRIX) $(R_SCRIPT)
	$(RSCRIPT) $(R_SCRIPT) $(COUNTS_MATRIX) $(RESULTS_DIR)

# Clean up intermediate files
clean:
	rm -rf counts $(RESULTS_DIR) design.csv

