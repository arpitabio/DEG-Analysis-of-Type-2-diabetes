# Differential Gene Expression Analysis using DESeq2
This project performs differential gene expression (DGE) analysis using the DESeq2 package in R. It identifies genes that are significantly upregulated or downregulated between control and diabetic samples.
## Overview
This project uses RNA-Seq count data to analyze differentially expressed genes (DEGs) using DESeq2 in R. The results are visualized with volcano plots, heatmaps, and PCA plots.
## Tools & Packages Used
- **DESeq2** – Differential expression analysis
- **ggplot2** – Data visualization
- **EnhancedVolcano** – Volcano plot generation
- **pheatmap** – Heatmap visualization
- **RColorBrewer** – Color palettes
## Installation
To install the required R packages, run the following command:

```r
install.packages(c("ggplot2", "pheatmap", "RColorBrewer"))
BiocManager::install(c("DESeq2", "EnhancedVolcano"))
```
## How to Run the Code
1. Clone the repository:
   ```bash
   git clone https://github.com/arpitabio/DEG-Analysis-of-Type-2-diabetes.git
2. Navigate into the project directory:
 `cd repo-name`
3. Open RStudio and run the script diabetes diabetes_R.R

## Workflow
1. Load Count Data
- Read count data matrix
- Assign gene symbols
- Remove unwanted columns

2. Create DESeq2 Object
- Define experimental conditions
- Construct DESeqDataSetFromMatrix()
- Filter low-expression genes

3. Differential Expression Analysis
- Run DESeq() to normalize counts and estimate dispersions
- Extract significant DEGs

4. Visualization
- Volcano plot: Log2 fold change vs. p-value
- Heatmap: Top 20 DEGs clustered
- PCA plot: Sample clustering
- MA plot: Average expression vs. log2 fold change

5. Functional Enrichment Analysis
- DAVID analysis for GO terms and KEGG pathways
- Generate pathway bar plots

## Files in This Repository

| File Name            | Description |
|----------------------|------------|
| `DEG_results.csv`    | Full differential expression results |
| `Upregulated_genes.csv` | List of significantly upregulated genes |
| `Volcano_Plot.png`   | Visualization of DEGs |
| `Heatmap.png`        | Heatmap of top 20 DEGs |

## Notes
- The dataset should contain raw count values for RNA-seq analysis.
- Pre-filtering is applied to remove lowly expressed genes.
- Adjust sample conditions as per the experimental design.

## References
DESeq2 Documentation: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
EnhancedVolcano: https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html
DAVID Functional Annotation Tool: https://david.ncifcrf.gov

## Contact
For any questions, feel free to reach out via [gangulyarpita561@gmail.com].

