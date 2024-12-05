# -------------------------------------------------------------------------
# 5 Exploratory data analysis
# -------------------------------------------------------------------------

# With the help of the DESeq2 Tutorial from Ashley Valentina Schwartz
# https://ashleyschwartz.com/posts/2023/05/deseq2-tutorial

library(DESeq2)
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

# Checking if the name of the columns in the counts matrix is the same as the names in the meta data file
all(colnames(count_data) %in% rownames(meta_data))
all(colnames(count_data) == rownames(meta_data))

# Set up DESeq objects ----------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = count_data, 
                              colData = meta_data, 
                              design = ~Condition)
dds
estimateSizeFactors(dds)

# Differential expression -------------------------------------------------

dds <- DESeq(dds)
res <- results(dds)
res

summary(res)
# Exact number of adjusted p-values that are <0.1
sum(res$padj < 0.1, na.rm=TRUE) 

# Exact number of adjusted p-values that are <0.05
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res$padj < 0.05, na.rm=TRUE)

# -------------------------------------------------------------------------
# 5.1 Visualising results
# -------------------------------------------------------------------------

library(ggrepel)
library(pheatmap)

# conversion of the data to a basic data frame
data <- data.frame(res)
head(data)

# PCA plot ----------------------------------------------------------------

rld <- rlog(dds, blind = TRUE)
plotPCA(rld, intgroup=c("Condition", "Organ"))

# Volcano plots -----------------------------------------------------------

# Add an additional column that identifies a gene as upregulated, downregulated, or unchanged
data <- data %>%
  mutate(
    Expression = case_when(log2FoldChange >= log(1) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(1) & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(data)

# Getting the top 10 most up- or downregulated genes
top <- 10

top_genes <- bind_rows(
  data %>%
    filter(Expression == 'Up-regulated') %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top),
  data %>%
    filter(Expression == 'Down-regulated') %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top)
)

# Dataframe with the top 10 genes
Top_Hits = head(arrange(data,pvalue),10)
Top_Hits

data$label = if_else(rownames(data) %in% rownames(Top_Hits), rownames(data), "")

# Basic plot
p1 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point( size = 2/5) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  xlim(-4.5, 4.5)
p1

# Basic plot with line + red for p < 0.05
p2 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point(aes(color = Expression), size = 2/5) +
  #geom_hline(yintercept= -log10(0.05), linetype="dashed", linewidth = .4) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  scale_color_manual(values = c("firebrick3", "black", "firebrick3")) +
  xlim(-4.5, 4.5) +
  theme(legend.position = "none")
p2

# With labels for top 10 sig overall
p3 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point(aes(color = Expression), size = 2/5) +
  # geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth = .4) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  scale_color_manual(values = c("firebrick3", "black", "firebrick3")) +
  xlim(-4.5, 4.5) +
  theme(legend.position = "none") +
  geom_text_repel(aes(label = label), size = 2.5)
p3

# Plot with up/down
p4 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point(aes(color = Expression), size = 2/5) +
  #geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth = .4) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
  xlim(-4.5, 4.5) 
  #ylim(0, 50)
p4

p5 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point(aes(color = Expression), size = 2/5) +
  # geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth = .4) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
  xlim(-4.5, 4.5) +
  geom_text_repel(aes(label = label), size = 2.5)
p5

# Heatmap -----------------------------------------------------------------

select <- order(rowMeans(count(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Condition","Organ")])
pheatmap(assay(dds)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# -------------------------------------------------------------------------
# 6 Differential expression analysis
# -------------------------------------------------------------------------