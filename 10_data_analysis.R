require(DESeq2)
require(tidyverse)

# data preparation
gene_counts <- 'gene_counts.txt'

count_data <- read.table(gene_counts, header=TRUE, stringsAsFactors=FALSE)
count_data <- count_data[, -c(2:6)]

count_data <- as.data.frame(count_data)
colnames(count_data) <- gsub("^X.*\\.mapping\\.(SRR[0-9]+)_sorted\\.bam$", "\\1", colnames(count_data))

rownames(count_data) <- NULL
rownames(count_data) <- count_data$Geneid
count_data <- count_data[, -which(names(count_data) == "Geneid")]

meta_data <- data.frame(read.csv2('metadata.csv', header=TRUE))

# samples_old <- colnames(count_data)
# samples_new <- setNames(meta_data$Condition, meta_data$Sample)
# colnames(count_data) <- samples_new[samples_old]

meta_data <- meta_data %>% remove_rownames %>% column_to_rownames(var="Sample")

all(colnames(count_data) %in% rownames(meta_data))
all(colnames(count_data) == rownames(meta_data))

# set up DESeq objects
dds <- DESeqDataSetFromMatrix(countData = count_data, 
                              colData = meta_data, 
                              design = ~Condition)
dds

# pre-filtering
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
# dds

# differential expression
dds <- DESeq(dds)
res <- results(dds)
res

summary(res)
sum(res$padj < 0.1, na.rm=TRUE) #exact number of adjusted p-values that are <0.01

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res$padj < 0.05, na.rm=TRUE)

# results visualisation
# convert results data to basic dataframe
data <- data.frame(res)
head(data)

# PCA
rld <- rlog(dds, blind = TRUE)
plotPCA(rld, intgroup=c("Condition", "Organ"))

# volcano plot
data <- data %>%
  mutate(
    Expression = case_when(log2FoldChange >= log(1) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(1) & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(data)

top <- 10
# we are getting the top 10 up and down regulated genes by filtering the column Up-regulated and Down-regulated and sorting by the adjusted p-value. 
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
# create a datframe just holding the top 10 genes
Top_Hits = head(arrange(data,pvalue),10)
Top_Hits

data$label = if_else(rownames(data) %in% rownames(Top_Hits), rownames(data), "")
# basic plot
p1 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point( size = 2/5) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  xlim(-4.5, 4.5)
p1

# basic plot with line + red for p < 0.05
p2 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point(aes(color = Expression), size = 2/5) +
  #geom_hline(yintercept= -log10(0.05), linetype="dashed", linewidth = .4) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  scale_color_manual(values = c("firebrick3", "black", "firebrick3")) +
  xlim(-4.5, 4.5) +
  theme(legend.position = "none")
p2

# with labels for top 10 sig overall
library(ggrepel)
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

# plot with up/down
p4 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point(aes(color = Expression), size = 2/5) +
  #geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth = .4) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
  xlim(-4.5, 4.5) 
  #ylim(0, 50)
p4

require(ggrepel)
p5 <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point(aes(color = Expression), size = 2/5) +
  # geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth = .4) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"P-Value")) +
  scale_color_manual(values = c("dodgerblue3", "black", "firebrick3")) +
  xlim(-4.5, 4.5) +
  geom_text_repel(aes(label = label), size = 2.5)
p5

# Heatmap
library("pheatmap")
select <- order(rowMeans(count(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Condition","Organ")])
pheatmap(assay(dds)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)