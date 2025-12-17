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
res_sig <- subset(res_overall_df, padj < 0.05)

# Reorder columns for IPA
res_ipa <- res_sig[, c("log2FoldChange", "pvalue", "padj", "gene")]

write.csv(res_ipa, "DESeq2_for_IPA.csv", row.names = FALSE)



