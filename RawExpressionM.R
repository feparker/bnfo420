library(DESeq2)
library(tidyverse)
library(ggplot2)

#read in raw counts 
counts <- read.delim(
  "M_annotated_raw_counts.tsv",
  header = TRUE,
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_ids <- counts[, 1]
counts <- counts[, -1]
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"
rownames(counts) <- gene_ids

#read in patient data 
meta <- read.delim(
  "M_covid_status_2.tsv",
  header = TRUE,
  stringsAsFactors = FALSE
)

rownames(meta) <- meta$PatientIDs
meta$PatientIDs <- NULL

# Clean up the trailing space in CovidStatus
meta$CovidStatus <- trimws(meta$CovidStatus)

# Check cleaned values
print("CovidStatus values after cleaning:")
print(unique(meta$CovidStatus))
print(table(meta$CovidStatus))

#create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ CovidStatus
)

# Set "non-long covid" as the reference (baseline)
# This means positive log2FC = higher in "long covid"
# Negative log2FC = higher in "non-long covid"
dds$CovidStatus <- relevel(dds$CovidStatus, ref = "non-long covid")

# Verify the reference level
print("Factor levels (first is reference):")
print(levels(dds$CovidStatus))

#pre-filter low counts 
dds <- dds[rowSums(counts(dds)) > 1, ]

#run DESeq
dds <- DESeq(dds)

# Get results - now comparing "long covid" vs "non-long covid"
res_overall <- results(dds)
summary(res_overall)

# Check significant genes in results
res_overall_df <- as.data.frame(res_overall)
res_overall_df$gene <- rownames(res_overall_df)
sig_genes_actual <- res_overall_df %>%
  filter(padj < 0.1 & !is.na(padj)) %>%
  arrange(padj)

print("Actually significant genes:")
print(sig_genes_actual[, c("gene", "log2FoldChange", "padj")])

# ============================================================================
# PLOT RAW COUNTS FOR SIGNIFICANT GENES
# ============================================================================

# List of significant genes from male analysis
sig_genes <- c("MIR8071-2")

# Get RAW counts from the DESeq2 object
raw_counts <- counts(dds, normalized = FALSE)

# Extract counts for significant genes (keep as matrix even for single gene)
sig_gene_counts <- raw_counts[sig_genes, , drop = FALSE]

# Convert to long format for plotting
plot_data <- as.data.frame(t(sig_gene_counts))
plot_data$Sample <- rownames(plot_data)
plot_data$CovidStatus <- meta[rownames(plot_data), "CovidStatus"]

# Reshape to long format
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = all_of(sig_genes),
    names_to = "Gene",
    values_to = "RawCounts"
  )

# Create individual plots for each gene
ggplot(plot_data_long, aes(x = CovidStatus, y = RawCounts, color = CovidStatus)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  facet_wrap(~Gene, scales = "free_y", ncol = 2) +
  scale_color_manual(
    values = c("non-long covid" = "#4393C3", "long covid" = "#D6604D"),
    name = "COVID Status"
  ) +
  labs(
    title = "Raw Expression of Significant Genes in Males",
    subtitle = "Long COVID vs Non-Long COVID",
    x = "COVID Status",
    y = "Raw Counts"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  )

ggsave("male_significant_genes_raw_counts.png", width = 10, height = 8, dpi = 300)

# Alternative: Create separate plots for each gene
for (gene in sig_genes) {
  gene_data <- plot_data_long %>% filter(Gene == gene)
  
  p <- ggplot(gene_data, aes(x = CovidStatus, y = RawCounts, color = CovidStatus)) +
    geom_jitter(width = 0.2, size = 4, alpha = 0.7) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.5) +
    scale_color_manual(
      values = c("non-long covid" = "#4393C3", "long covid" = "#D6604D"),
      name = "COVID Status"
    ) +
    labs(
      title = paste("Raw Expression of", gene, "in Males"),
      subtitle = "Long COVID vs Non-Long COVID",
      x = "COVID Status",
      y = "Raw Counts"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
      legend.position = "bottom"
    )
  
  ggsave(paste0("male_", gene, "_raw_counts.png"), plot = p, width = 7, height = 6, dpi = 300)
}

# Print summary statistics
print("Summary statistics for each gene:")
for (gene in sig_genes) {
  cat("\n", gene, ":\n")
  gene_summary <- plot_data_long %>%
    filter(Gene == gene) %>%
    group_by(CovidStatus) %>%
    summarise(
      Mean = mean(RawCounts),
      Median = median(RawCounts),
      SD = sd(RawCounts),
      N = n()
    )
  print(gene_summary)
}