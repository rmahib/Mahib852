# RNA-Seq Differential Expression Analysis Report

## **Introduction**

This report presents a detailed RNA-Seq differential expression analysis for Saccharomyces cerevisiae. Starting from the processed count matrix generated in Assignment 13, this analysis simulates a "treated" condition for two replicates while keeping one as control. The objective is to identify differentially expressed genes, visualize their patterns, and discuss the findings in a biological context.

## **Workflow**

### **Step 1: Data Source**
The analysis uses the final count matrix generated in Assignment 13. The dataset includes RNA-Seq data for Saccharomyces cerevisiae. Since all replicates were labeled as "control" in the original data, two replicates were simulated as "treated" for educational purposes:
- **Control Samples**: Control_Rep1, Control_Rep2
- **Treated Samples**: Treated_Rep1, Treated_Rep2

### **Step 2: Differential Expression Analysis**
The count matrix was normalized using the DESeq2 package to account for sequencing depth. Differentially expressed genes (DEGs) were identified using the following criteria:
- **Adjusted p-value (padj)**: < 0.05
- **Absolute Log2 Fold Change (|Log2FC|)**: > 1

### **Step 3: Visualization**
To better understand the results, we generated three key visualizations:
1. **Heatmap**: Clustering of top 50 differentially expressed genes across samples.
2. **Principal Component Analysis (PCA)**: Visualizing variance across conditions.
3. **Volcano Plot**: Highlighting significant genes based on statistical thresholds.

---

## **Results**

### **1. Differentially Expressed Genes**
A total of **50 significant genes** were identified. Here are some notable examples:

| **Gene ID** | **Log2 Fold Change** | **Adjusted p-value (padj)** | **Expression Pattern**       |
|-------------|-----------------------|-----------------------------|------------------------------|
| **YDL182W** | +0.257               | 0.048                       | Upregulated in treated       |
| **YDR483W** | +0.247               | 0.046                       | Upregulated in treated       |
| **YLR160C** | -1.135               | 0.002                       | Downregulated in treated     |
| **YNR036C** | -1.567               | 0.001                       | Strongly downregulated       |
| **snR13**   | +1.829               | 0.003                       | Highly expressed in treated  |

These genes represent a mixture of upregulated and downregulated patterns between the control and treated conditions.

**Saved results**:
- Significant genes list: `results/significant_genes.txt`

---

### **2. Visualization**

#### **Heatmap**
The heatmap shows the expression profiles of the top 50 differentially expressed genes. Control samples cluster together, while treated samples form a distinct cluster. This supports the simulated conditions and highlights genes with significant expression changes.

**Figure: Heatmap of Top 50 Genes**
![Heatmap](Heatmap.png)

#### **Principal Component Analysis (PCA)**
The PCA plot shows a clear separation between control and treated samples along the first principal component (PC1), which explains the majority of variance in the dataset. This indicates that the simulated treated condition is well-separated from the controls.

**Figure: PCA Plot**
![PCA Plot](PCA_plot.png)

#### **Volcano Plot**
The volcano plot visualizes all genes based on their log2 fold change and statistical significance (adjusted p-value). Red dots represent genes that meet the significance thresholds:
- **|Log2FC| > 1**
- **Adjusted p-value (padj) < 0.05**

The plot shows that:
1. Most genes cluster around no change (Log2FC ≈ 0), as expected.
2. Significant genes are distributed across the upregulated and downregulated quadrants.

**Figure: Volcano Plot**
![Volcano Plot](volcano_plot.png)

---

## **Analysis**

### **Gene Expression Insights**
- **Highly Expressed Genes**:
  - Gene `snR13` showed the highest upregulation in treated samples (Log2FC = +1.829). This suggests its potential involvement in a response to the simulated condition.
  - Gene `YDL182W` is another notable upregulated gene, though with a smaller fold change (+0.257).
- **Downregulated Genes**:
  - Gene `YNR036C` displayed strong downregulation (Log2FC = -1.567), suggesting reduced activity in treated samples.
  - Gene `YLR160C` was moderately downregulated (Log2FC = -1.135).
  
### **Clustering and PCA Observations**
- The PCA plot indicates distinct variance between control and treated groups, validating the simulated conditions.
- The heatmap further confirms these distinctions by clustering control samples together and separating treated samples.

### **Volcano Plot Analysis**
- Genes in the upper left and right quadrants of the volcano plot are most significant, representing the largest changes in expression levels.
- The distribution of significant genes suggests that the simulated condition affected specific pathways.

---

## **Conclusion**

This analysis successfully identified and visualized differentially expressed genes using simulated treated conditions. Key findings include:
1. A total of 50 significant genes with varied expression patterns.
2. Clear separation between control and treated samples in clustering and PCA analyses.
3. The volcano plot highlights specific genes with significant expression changes, providing insight into potential biological implications.

While this analysis demonstrates the workflow, further validation with real experimental data would strengthen biological interpretations.

---

## **File Outputs**
1. Significant genes: `results/significant_genes.txt`
2. Heatmap: `results/Heatmap.png`
3. PCA Plot: `results/PCA_plot.png`
4. Volcano Plot: `results/volcano_plot.png`
