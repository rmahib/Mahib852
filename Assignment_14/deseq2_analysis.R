args <- commandArgs(trailingOnly = TRUE)

# Input count matrix and output directory
input_file <- args[1]
output_dir <- args[2]

# Load count data
count_data <- read.table(input_file, header = TRUE, row.names = 1, check.names = FALSE)

# Prepare sample information
sample_info <- data.frame(
  Sample = colnames(count_data),
  Condition = c("Control", "Control", "Treated", "Treated")
)

# Differential expression analysis
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = ~Condition)
dds <- DESeq(dds)
res <- results(dds)

# Save significant genes
sig_genes <- res[which(res$padj < 0.05), ]
write.table(as.data.frame(sig_genes), file.path(output_dir, "significant_genes.txt"), sep = "\t", quote = FALSE)

# Save Volcano Plot
library(ggplot2)
res$significant <- ifelse(res$padj < 0.05, "Yes", "No")
volcano_plot <- ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value")
ggsave(file.path(output_dir, "Volcano_plot.png"), volcano_plot)

# Save PCA Plot
normalized_counts <- vst(dds)
pca_data <- plotPCA(normalized_counts, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
ggsave(file.path(output_dir, "PCA_plot.png"), pca_plot)

# Save Heatmap
library(pheatmap)
top_genes <- head(order(rowSums(count_data), decreasing = TRUE), 50)
top_gene_data <- count_data[top_genes, ]
pheatmap(as.matrix(top_gene_data),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Heatmap of Top Differentially Expressed Genes",
         fontsize_row = 8,
         fontsize_col = 8,
         filename = file.path(output_dir, "Heatmap.png"))
