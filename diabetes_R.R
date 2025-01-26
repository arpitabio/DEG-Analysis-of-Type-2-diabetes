# Set your working directory to the folder where the dataset is located
setwd("D:/Diabetes_Profiling")

# Load the required libraries
library(DESeq2)
library(dplyr)

# 1. Load the dataset (assuming it's a tab-separated file)
data <- read.csv("D:/Diabetes_Profiling/raw_counts.csv",sep = ",", header = TRUE, row.names = 1)
data_plot <- read.csv("D:/Diabetes_Profiling/raw_counts.csv",sep = ",", header = TRUE)

# Check the first few rows to ensure it loaded correctly
head(data)
ncol(data)

# Separate the symbol column
symbols <- data$Gene_Name  # Store the symbol column separately
symbols
# Remove the symbol column for DESeq2 input
data <- data[, -1]  # Remove the first column (symbol column)
data

# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)

# 2. Create a condition vector (manually assign based on your experimental design)
# Example: If there are 3 control and 3 treatment samples
conditions <- factor(c(rep("control", 4), rep("diabetes", 4)))  # Adjust to your actual data
length(conditions)

# 3. Create DESeq2 object
summary(data)
data <- round(data)  # Round to nearest integer
str(data)  # This should show integer type for count data

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = DataFrame(condition = conditions),
                              design = ~ condition)

# 4. Pre-filtering: Remove genes with low counts (optional)
dds <- dds[rowSums(counts(dds)) > 1, ]

# 5. Perform DESeq2 differential expression analysis
dds <- DESeq(dds)

# 6. Extract results
res <- results(dds)

# Display summary of results
summary(res)

# Convert DESeq2 results to a data frame
res_df <- as.data.frame(res)

colnames(res_df)
colnames(data_plot)

# Add Ensembl IDs as a column for merging
res_df$ENSEMBL <- rownames(res_df)
# Merge res_df with data_plot (assuming data_plot contains ENSEMBL and SYMBOL columns)
res_df <- merge(res_df, data_plot[, c("ENSEMBL", "Gene_Name")], by = "ENSEMBL", all.x = TRUE)



# Check the updated results
head(res_df)
res_df

# 7. Volcano Plot (Visualizing p-values and log2 fold changes)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano Plot of DEGs')

# Volcano Plot with Gene Symbols
EnhancedVolcano(res_df,
                lab = res_df$Gene_Name,  # Use Gene_name column from res_df
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano Plot of DEGs')


# 8. Heatmap of Top 20 DEGs (based on smallest p-value)
# Get results from DESeq2 analysis (e.g., for condition comparison)
res <- results(dds)

# Sort results by p-value or log-fold change, and select the top genes
top_genes <- order(res_df$padj)[1:20]  # Example: top 100 genes by adjusted p-value

# Subset normalized counts for top genes
data_top_genes <- counts(dds, normalized = TRUE)[top_genes, ]

# Check for missing or infinite values and replace them with 0
data_top_genes[is.na(data_top_genes)] <- 0
data_top_genes[is.infinite(data_top_genes)] <- 0

# Ensure data is in matrix form
data_top_genes <- as.matrix(data_top_genes)

# Load RColorBrewer package
library(RColorBrewer)

# Define a color palette for better visualization
color_palette <- colorRampPalette(brewer.pal(11, "RdBu"))(100)

# Generate heatmap with a different color palette
pheatmap(data_top_genes,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = FALSE,
         color = color_palette)

# Add annotations and generate heatmap again
annotation_col <- data.frame(condition = factor(conditions))
rownames(annotation_col) <- colnames(data_top_genes)

pheatmap(data_top_genes,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         scale = "row",
         show_rownames = FALSE,
         color = color_palette)


# Load RColorBrewer package
library(RColorBrewer)

# Check if color palette is generated correctly
color_palette <- colorRampPalette(brewer.pal(9, "Blues"))(100)
length(conditions)
ncol(data_top_genes)


# Generate heatmap
pheatmap(data_top_genes,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = FALSE,
         color = color_palette)
annotation_col <- data.frame(condition = factor(conditions))
rownames(annotation_col) <- colnames(data_top_genes)

pheatmap(data_top_genes,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         scale = "row",
         show_rownames = FALSE,
         color = color_palette)



# 9. PCA Plot: Visualizing sample clustering based on gene expression
# Perform variance stabilizing transformation
# Perform rlog transformation
rld <- rlog(dds, blind = FALSE)

# Perform PCA on rlog-transformed data
pca_rlog <- prcomp(t(assay(rld)))

# Get percentage of variance explained by each principal component
percent_variance_rlog <- (pca_rlog$sdev^2 / sum(pca_rlog$sdev^2)) * 100

vsd <- vst(dds, blind = FALSE)

# Perform PCA using prcomp
pca <- prcomp(t(assay(vsd)))

# Get percentage of variance explained by each principal component
percent_variance <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

# Prepare PCA data for plotting
pca_data <- as.data.frame(pca$x)
pca_data$condition <- colData(vsd)$condition  # Add condition information

# Plot PCA with % variance in labels
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(
    title = "PCA of Samples",
    x = paste0("PC1 (", round(percent_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_variance[2], 1), "%)")
  ) +
  theme_minimal()
#pca increasing components

# Get all percentage of variance explained
percent_variance_all <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

# Prepare PCA data for plotting
pca_data <- as.data.frame(pca$x)
pca_data$condition <- colData(vsd)$condition  # Add condition information

# Create a cumulative variance plot
ggplot(data.frame(PC = 1:length(percent_variance_all), Variance = cumsum(percent_variance_all)),
       aes(x = PC, y = Variance)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Cumulative Variance Explained by PCA",
    x = "Principal Component",
    y = "Cumulative Variance (%)"
  ) +
  theme_minimal()

# Optionally, you could use more components for plotting (PC1 to PC5 or more)
ggplot(pca_data, aes(x = PC1, y = PC3, color = condition)) +  # PC3 instead of PC2
  geom_point(size = 3) +
  labs(
    title = "PCA of Samples",
    x = paste0("PC1 (", round(percent_variance[1], 1), "%)"),
    y = paste0("PC3 (", round(percent_variance[3], 1), "%)")
  ) +
  theme_minimal()
#pca
library(plotly)

# Prepare PCA data
pca_data <- as.data.frame(pca$x)
pca_data$condition <- colData(vsd)$condition  # Add condition information

# Create 3D PCA plot
fig <- plot_ly(
  pca_data,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~condition,
  colors = "Set1",
  type = "scatter3d",
  mode = "markers"
) %>%
  layout(
    title = "3D PCA of Samples",
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(percent_variance[1], 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(percent_variance[2], 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(percent_variance[3], 1), "%)"))
    )
  )

fig


# Create a volcano plot
library(ggplot2)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(padj)") +
  theme(legend.position = "none")

# Obtain DESeq2 results
res <- results(dds)
res
# Filter upregulated genes (logFC > 1, padj < 0.05)
upregulated_genes <- res_df[which(res_df$log2FoldChange > 1 & res_df$padj < 0.05), ]

# Filter downregulated genes (logFC < -1, padj < 0.05)
downregulated_genes <- res_df[which(res_df$log2FoldChange < -1 & res_df$padj < 0.05), ]
# View upregulated genes
upregulated_genes

# View downregulated genes
downregulated_genes

str(upregulated_genes)
str(downregulated_genes)

#MA PLOT
#Load necessary libraries
library(ggplot2)

# Example data: You can replace this with your own results from DESeq2 or similar
# Sample data with log2FoldChange and padj values
res <- data.frame(
  log2FoldChange = rnorm(1000),  # Example log2 fold change values
  padj = runif(1000)  # Example adjusted p-values
)

# Calculate the average expression (example: just a random value for demonstration)
# In real cases, this would be the average expression of each gene across samples
res$avg_expression <- runif(1000, 1, 100)  # Replace with real average expression values

# Create the MA plot
ggplot(res, aes(x = avg_expression, y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "MA Plot", x = "Average Expression", y = "Log2 Fold Change") +
  theme(legend.position = "none")

# 10. Save DEG results to CSV
write.csv(as.data.frame(res_df), "DEG_results.csv")
write.csv(as.data.frame(upregulated_genes), "Upregulated_genes.csv")
write.csv(as.data.frame(downregulated_genes), "Downregulated_genes.csv")

#GO Pathway DAVID
# Replace the path with your downloaded file
data <- read.delim("https://davidbioinformatics.nih.gov/data/download/chart_75B1DEC91F701737674267483.txt",
                   header = TRUE, sep = "\t")

# Filter significant terms
significant_terms <- data[data$PValue < 0.05, ]

#Bar plot 
library(ggplot2)

# Plot Fold Enrichment for significant terms
# Subset top 20 pathways based on fold enrichment
top_20_terms <- significant_terms[order(-significant_terms$Fold.Enrichment),][1:20,]

# Plot the bar plot
ggplot(top_20_terms, aes(x = reorder(Term, Fold.Enrichment), y = Fold.Enrichment, fill = -log10(PValue))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(x = "Pathways", y = "Fold Enrichment", title = "Top 20 Enriched Pathways (p < 0.05)", fill = "-log10(PValue)") +
  theme(axis.text.y = element_text(size = 12, hjust = 0), # Adjust text size
        plot.title = element_text(size = 10))



#DAVID KEGG Pathway
# URLs for the files
chart_url <- "https://davidbioinformatics.nih.gov/data/download/chart_57DAF60854C21737874240485.txt"

# Download files locally
download.file(chart_url, destfile = "chart.txt", method = "curl")

# Read the downloaded files
chart_data <- read.delim("chart.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Inspect the first few rows
head(chart_data)

# Filter for KEGG pathways
kegg_data <- subset(chart_data, Category == "KEGG_PATHWAY")

# Add -log10(p-value) for visualization
kegg_data$logP <- -log10(kegg_data$PValue)

# Select relevant columns
kegg_data <- kegg_data[, c("Term", "logP", "Count", "Fold.Enrichment")]

# View the processed data
head(kegg_data)

library(ggplot2)

# Example: Bar plot for top KEGG pathways
library(ggplot2)

# Bar plot with color gradient based on -log10(p-value)
ggplot(kegg_data, aes(x = reorder(Term, -logP), y = logP, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip the plot for better readability
  scale_fill_gradient(low = "red", high = "blue") +  # Gradient color scale (red to blue)
  labs(title = "Top Enriched KEGG Pathways",
       x = "Pathway", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))

#DAVID GO_BP
# Filter for GO biological processes
go_bp_data <- subset(chart_data, Category == "GOTERM_BP_DIRECT")

# Add -log10(p-value) for visualization
go_bp_data$logP <- -log10(go_bp_data$PValue)

# Select relevant columns
go_bp_data <- go_bp_data[, c("Term", "logP", "Count", "Fold.Enrichment")]

# View the processed data
head(go_bp_data)


library(ggplot2)
# Scatter plot of Fold Enrichment vs. -log10(p-value)
ggplot(go_bp_data, aes(x = Fold.Enrichment, y = logP, label = Term)) +
  geom_point(color = "tomato", size = 3) +
  geom_text(aes(label = ifelse(logP > 5, as.character(Term), "")), 
            hjust = 0, vjust = 0, size = 3, color = "black") +  # Annotate top terms
  labs(title = "GO Biological Process: Fold Enrichment vs. -log10(p-value)",
       x = "Fold Enrichment", 
       y = "-log10(p-value)") +
  theme_minimal()

library(ggplot2)

# Select the top 20 pathways based on -log10(p-value)
top_20_go_bp_data <- go_bp_data[order(go_bp_data$logP, decreasing = TRUE),][1:20, ]

# Bar plot with color gradient based on -log10(p-value)
ggplot(top_20_go_bp_data, aes(x = reorder(Term, -logP), y = logP, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip the plot for better readability
  scale_fill_gradient(low = "green", high = "blue") +  # Gradient color scale (red to blue)
  labs(title = "Top 20 GO Biological Process Enrichment",
       x = "GO Term", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
















