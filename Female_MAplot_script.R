library(DESeq2)
library(tidyverse)
library(ggplot2)

#read in raw counts 
counts <- read.delim(
  "F_annotated_raw_counts.tsv",
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
  "F_covid_status_2.tsv",
  header = TRUE,
  stringsAsFactors = FALSE
)

rownames(meta) <- meta$PatientIDs
meta$PatientIDs <- NULL

# Clean up the trailing space in CovidStatus (optional but recommended)
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

# Create results dataframe
res_overall_df <- as.data.frame(res_overall)
res_overall_df$gene <- rownames(res_overall_df)
res_overall_df$significant <- res_overall_df$padj < 0.05 & !is.na(res_overall_df$padj)

sig_genes <- res_overall_df %>%
  filter(significant == TRUE)

print("Significant genes (positive = upregulated in long covid):")
print(sig_genes[, c("gene", "log2FoldChange", "padj", "baseMean")])

# Create MA plot
ggplot(res_overall_df, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = significant), alpha = 0.5, size = 1) +
  scale_color_manual(
    values = c("grey50", "red"), 
    labels = c("Not significant", "FDR < 0.05"),
    name = "Significance"
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue", linewidth = 0.5) +
  geom_text(
    data = sig_genes,
    aes(label = gene),
    size = 3.5,
    fontface = "bold",
    vjust = -0.5,
    hjust = 0.5,
    check_overlap = TRUE
  ) +
  labs(
    title = "MA Plot: Long COVID vs Non-Long COVID",
    subtitle = paste0(sum(res_overall_df$significant, na.rm = TRUE), 
                      " significant genes (FDR < 0.05)"),
    x = "Mean of raw counts (log10 scale)",
    y = "Log2 Fold Change"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )