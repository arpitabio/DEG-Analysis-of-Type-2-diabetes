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


