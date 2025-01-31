# -------------------------------------------------------------------------
# 5 Exploratory data analysis
# -------------------------------------------------------------------------

library(tidyverse) 

# Data preparation --------------------------------------------------------

gene_counts <- 'gene_counts.txt'

# Reformatting of the table to correspond to the format expected by DESeq2 
count_data <- read.table(gene_counts, header=TRUE, stringsAsFactors=FALSE)
count_data <- count_data[, -c(2:6)]

count_data <- as.data.frame(count_data)
colnames(count_data) <- gsub("^X.*\\.mapping\\.(SRR[0-9]+)_sorted\\.bam$", "\\1", colnames(count_data))

rownames(count_data) <- NULL
rownames(count_data) <- count_data$Geneid
count_data <- count_data[, -which(names(count_data) == "Geneid")]

# Formatting a meta data table
meta_data <- data.frame(read.csv2('metadata.csv', header=TRUE))
meta_data <- meta_data %>% remove_rownames %>% column_to_rownames(var="Sample")
meta_data$Condition <- as.factor(meta_data$Condition)
meta_data$Organ <- as.factor(meta_data$Organ)

# Checking if the name of the columns in the counts matrix is the same as the names in the meta data file
all(colnames(count_data) %in% rownames(meta_data))
all(colnames(count_data) == rownames(meta_data))

# -------------------------------------------------------------------------
# 5.1 Differential expression analysis
# -------------------------------------------------------------------------

library(DESeq2)

# Set up DESeq objects ----------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = count_data, 
                              colData = meta_data, 
                              design = ~ Condition * Organ) 
# The design includes Condition * Organ, which models both the main effects 
# (Condition, Organ) and their interaction.

dds
estimateSizeFactors(dds)

# -------------------------------------------------------------------------
# 5.2 Visualising results
# -------------------------------------------------------------------------

library(pheatmap)

# PCA plot ----------------------------------------------------------------

# Removing the dependence of the variance on the mean
vst <- vst(dds, blind = TRUE)
plotPCA(vst, intgroup=c("Condition", "Organ"))

pcaData <- plotPCA(vst, intgroup = c("Condition", "Organ"), returnData = TRUE)

# Extraction of the percentage of variance for each PC
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Creation of a ggplot with Condition mapped to color, Organ to shape
ggplot(pcaData, aes(x = PC1, y = PC2, color = Condition, shape = Organ)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_minimal() +
  ggtitle("PCA based on the samples gene expression profile")

# Alternatively
# rld <- rlog(dds, blind = TRUE)
# plotPCA(rld, intgroup=c("Condition", "Organ"))

# Heatmap -----------------------------------------------------------------

Realisation of a Heatmap graph
sample_dists <- dist(t(assay(vst)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- paste(meta_data$Condition, meta_data$Organ, sep = "_")
colnames(sample_dist_matrix) <- rownames(sample_dist_matrix)

pheatmap(sample_dist_matrix, clustering_distance_rows = sample_dists, 
         clustering_distance_cols = sample_dists)

# -------------------------------------------------------------------------
# 6 Differential expression analysis
# -------------------------------------------------------------------------

dds <- DESeq(dds)
res <- results(dds)

summary(res)

# Exact number of adjusted p-values that are <0.05 for the all dataset
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res$padj < 0.05, na.rm=TRUE)

plotMA(res)

# Comparison of infected vs. control within lung samples ------------------

dds_lung <- dds[, meta_data$Organ == "lung"]

dds_lung <- DESeqDataSetFromMatrix(
  countData = counts(dds_lung),
  colData = colData(dds_lung),
  design = ~ Condition
)

dds_lung <- DESeq(dds_lung)

res_lung <- results(dds_lung, contrast = c("Condition", "infected", "control"))
res_lung
res_interaction <- results(dds, name = "Conditioninfected.Organlung")

summary(res_lung)

# Number of significant genes at padj < 0.05
sum(res_lung$padj < 0.05, na.rm = TRUE)

# Separate upregulated and downregulated
upregulated <- sum(res_lung$padj < 0.05 & res_lung$log2FoldChange > 0, na.rm = TRUE)
downregulated <- sum(res_lung$padj < 0.05 & res_lung$log2FoldChange < 0, na.rm = TRUE)

cat("Upregulated:", upregulated, "\nDownregulated:", downregulated)

res_lung_filtered <- res_lung[!is.na(res_lung$padj), ]

# Select significant genes with padj < 0.05
significant_genes <- res_lung_filtered[res_lung_filtered$padj < 0.05, ]

#####

library(AnnotationDbi)
library(org.Mm.eg.db)

# Convert Ensembl IDs to gene symbols
gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = rownames(significant_genes),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Add gene symbols to the result
significant_genes$Gene_Symbol <- gene_symbols[rownames(significant_genes)]

# Check if it works
head(significant_genes)

#####

# Analyze the significant genes

head(significant_genes[order(significant_genes$log2FoldChange, decreasing = TRUE), ]) # Top upregulated
head(significant_genes[order(significant_genes$log2FoldChange, decreasing = FALSE), ]) # Top downregulated

# Volcano Plot

res_lung_df <- as.data.frame(res_lung)
res_lung_df <- res_lung_df[!is.na(res_lung_df$padj), ]
res_lung_df$Significance <- "Not Significant"
res_lung_df$Significance[res_lung_df$padj < 0.05 & res_lung_df$log2FoldChange > 1] <- "Upregulated" # expression has doubled
res_lung_df$Significance[res_lung_df$padj < 0.05 & res_lung_df$log2FoldChange < -1] <- "Downregulated" # expression has halved

ggplot(res_lung_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot - Lung",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Regulation"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")

# Plot counts

top5_upregulated <- rownames(res_lung[order(res_lung$log2FoldChange, decreasing = TRUE), ])[1:5]
top5_downregulated <- rownames(res_lung[order(res_lung$log2FoldChange), ])[1:5]

# Combine into one vector
top_genes <- c(top5_upregulated, top5_downregulated)

# View the selected top genes
print(top_genes)

library(gridExtra)

# Create a list to store the plots
plot_list <- list()

for (gene in top_genes) {
  # Create a count plot for the current gene
  plot <- plotCounts(dds_lung, gene = gene, intgroup = "Condition", returnData = TRUE)
  
  # Generate the ggplot2 object
  p <- ggplot(plot, aes(x = Condition, y = count)) + 
    geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) + 
    scale_y_log10(breaks = c(25, 100, 400)) + 
    ggtitle(paste("Gene Expression for", gene))
  
  # Store the plot in the list
  plot_list[[gene]] <- p
}

# Arrange the plots in a 2x5 grid
grid.arrange(grobs = plot_list, nrow = 2, ncol = 5)

# Selection of the genes of interest
selected_genes <- c("Tgtp1", "Tgtp2", "Batf2", "Gbp5")

selected_ensembl_ids <- rownames(res_lung)[res_lung$gene_symbol %in% selected_genes]

plot_list <- list()

for (gene in selected_ensembl_ids) {
  # Map Ensembl ID to gene symbol
  gene_symbol <- res_lung$gene_symbol[rownames(res_lung) == gene]
  
  # Create a count plot for the current gene
  plot <- plotCounts(dds_lung, gene = gene, intgroup = "Condition", returnData = TRUE)
  
  # Generate the plot
  p <- ggplot(plot, aes(x = Condition, y = count)) + 
    geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) + 
    scale_y_log10(breaks = c(25, 100, 400, 1000, 5000, 10000, 50000)) + 
    ggtitle(paste("Expression of", gene_symbol))
  
  # Store the plot in the list
  plot_list[[gene_symbol]] <- p
}

# Arrange the plots in a single row (or adjust as needed)
grid.arrange(grobs = plot_list, nrow = 1, ncol = length(plot_list))

# Originally, other comparisions such as infected vs. control within blood samples, lung vs. blood within infected samples,
# lung vs. blood within control samples and a global comparison: Infected vs. Control (ignoring tissue) were also made but
# not analysed due to time limit and word requirement of the final report.

# -------------------------------------------------------------------------
# 7 Overrepresentation analysis
# -------------------------------------------------------------------------

library(clusterProfiler)
library(org.Mm.eg.db) 
library(AnnotationDbi)

# Comparison: Infected vs. Control in Lung --------------------------------

# Filter for significant genes (padj < 0.05)
de_genes_lung <- rownames(res_lung)[!is.na(res_lung$padj) & res_lung$padj < 0.05]

# Prepare the universe (all genes measured in the analysis)
all_genes_lung <- rownames(res_lung)

# Run GO enrichment analysis
ego_lung <- enrichGO(
  gene          = de_genes_lung,     # Significant DE genes
  universe      = all_genes_lung,    # All measured genes
  OrgDb         = org.Mm.eg.db,      # Mouse annotation database
  keyType       = "ENSEMBL",         # Format of gene IDs
  ont           = "BP",              # GO subontology: "BP", "MF", "CC", or "ALL"
  pAdjustMethod = "BH",              # Adjust p-values for multiple testing
  pvalueCutoff  = 0.05,              # p-value threshold
  qvalueCutoff  = 0.2,               # q-value threshold
  readable      = TRUE               # Convert gene IDs to gene symbols for readability
)

# View summary of results
head(ego_lung)

# Bar plot of top 15 enriched GO terms
barplot(ego_lung, showCategory = 15, title = "Top 15 Enriched GO Terms")

# Dot plot of top 15 enriched GO terms
dotplot(ego_lung, showCategory = 15, title = "Top 15 Enriched GO Terms")

# Convert results to a data frame
ego_lung_df <- as.data.frame(ego_lung)

# Save results to a CSV file
write.csv(ego_lung_df, file = "GO_enrichment_results_lung.csv", row.names = FALSE)