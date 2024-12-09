args <- commandArgs(trailingOnly = TRUE)

# Input count matrix and output directory
input_file <- args[1]
output_dir <- args[2]

# Load count data
count_data <- read.table(input_file, header = TRUE, row.names = 1, check.names = FALSE)
count_data <- count_data[, -c(1:5)]  # Remove unnecessary columns

# Normalize data
normalized_data <- t(t(count_data) / colSums(count_data)) * 1e6

# Identify consistently expressed genes
consistent_genes <- normalized_data[rowSums(normalized_data > 1) == ncol(normalized_data), ]

# Save consistent genes
write.table(consistent_genes, file.path(output_dir, "consistent_genes.txt"), sep = "\t", quote = FALSE)

# Generate a boxplot for gene expression levels across samples
png(filename = file.path(output_dir, "gene_expression_boxplot.png"), width = 800, height = 800)
boxplot(as.data.frame(count_data),
        main = "Gene Expression Distribution Across Samples",
        xlab = "Samples",
        ylab = "Read Counts",
        col = rainbow(ncol(count_data)))
dev.off()

