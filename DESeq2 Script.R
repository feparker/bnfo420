###############################################################
# 1. Load Libraries
###############################################################
library(DESeq2)
library(readxl)
library(dplyr)
library(tibble)
library(EnhancedVolcano)
library(pheatmap)


###############################################################
# 2. Detect and Load Annotated Count File
###############################################################
files <- list.files("C:/DESeq2", full.names = TRUE)

annot_file <- files[grep("\\.tsv$", files)]
annot_file <- annot_file[grep("annot|raw|count", annot_file, ignore.case = TRUE)]

if (length(annot_file) == 0) {
  stop("No annotated .tsv file found in C:/DESeq2")
}

cat("Using annotated file:", annot_file, "\n")


###############################################################
# 3. Load Annotated Counts
###############################################################
annot <- read.table(
  annot_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)


###############################################################
# 4. Prepare Count Matrix
###############################################################
count_cols <- grep("^GSM", colnames(annot), value = TRUE)

counts <- annot[, c("Symbol", count_cols)]
counts <- counts[!is.na(counts$Symbol), ]
counts <- counts[!duplicated(counts$Symbol), ]

rownames(counts) <- counts$Symbol
counts$Symbol <- NULL


###############################################################
# 5. Load Metadata
###############################################################
meta <- read_excel("C:/DESeq2/Copy of covid_status.xlsx")

meta$Covid_status <- trimws(meta$`Covid status`)
meta$Covid_status <- factor(meta$Covid_status,
                            levels = c("non-long covid", "long covid"))
meta$Gender <- factor(meta$Gender)

rownames(meta) <- meta$`Patient identifier`


###############################################################
# 6. Align Counts and Metadata
###############################################################
counts <- counts[, rownames(meta)]


###############################################################
# 7. Build DESeq2 Object and Run DESeq
###############################################################
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ Gender + Covid_status
)

dds <- DESeq(dds)


###############################################################
# 8. Male-Only Differential Expression
###############################################################
dds_male <- dds[, dds$Gender == "M"]
dds_male$Gender <- droplevels(dds_male$Gender)
design(dds_male) <- ~ Covid_status
dds_male <- DESeq(dds_male)

res_male <- results(dds_male,
                    contrast = c("Covid_status", "long covid", "non-long covid"))

write.csv(as.data.frame(res_male),
          "annotated_DEG_male_LC_vs_nonLC.csv")


###############################################################
# 9. Female-Only Differential Expression
###############################################################
dds_female <- dds[, dds$Gender == "F"]
dds_female$Gender <- droplevels(dds_female$Gender)
design(dds_female) <- ~ Covid_status
dds_female <- DESeq(dds_female)

res_female <- results(dds_female,
                      contrast = c("Covid_status", "long covid", "non-long covid"))

write.csv(as.data.frame(res_female),
          "annotated_DEG_female_LC_vs_nonLC.csv")


###############################################################
# 10. Volcano Plots (Male + Female)
###############################################################
# --- Larger plots ---
par(mfrow = c(1, 1))
options(repr.plot.width = 12, repr.plot.height = 10)

# --- Volcano: Male ---
EnhancedVolcano(
  res_male,
  lab = rownames(res_male),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Volcano: Male',
  pCutoff = 0.05,
  FCcutoff = 1,
  labSize = 5,
  titleLabSize = 20,
  axisLabSize = 16,
  legendLabSize = 14,
  colAlpha = 0.8,
  pointSize = 3.5,
  drawConnectors = TRUE
)

# --- Volcano: Female ---
EnhancedVolcano(
  res_female,
  lab = rownames(res_female),
  x = 'log2FoldChange',
  y = 'padj',
  title = 'Volcano: Female',
  pCutoff = 0.05,
  FCcutoff = 1,
  labSize = 5,
  titleLabSize = 20,
  axisLabSize = 16,
  legendLabSize = 14,
  colAlpha = 0.8,
  pointSize = 3.5,
  drawConnectors = TRUE
)

###############################################################
# 11. MA Plots (Male + Female)
###############################################################
plotMA(
  res_male,
  ylim = c(-5, 5),
  main = "MA Plot: Male",
  alpha = 0.05
)

plotMA(
  res_female,
  ylim = c(-5, 5),
  main = "MA Plot: Female",
  alpha = 0.05
)


###############################################################
# 12. Heatmaps (Male + Female)
###############################################################
# Reset plot window (heatmap needs full window)
par(mfrow = c(1, 1))

norm <- read.table(
  "C:/DESeq2/GSE224615_gene_names (1).tsv",
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

# Remove duplicates BEFORE setting rownames
norm <- norm[!duplicated(norm$Gene), ]

rownames(norm) <- norm$Gene
norm$Gene <- NULL


###############################################################
# 12A. Build Annotation Table (LC vs R + Sex)
###############################################################
annotation_col <- data.frame(
  Condition = meta$Covid_status,
  Sex = meta$Gender
)
rownames(annotation_col) <- rownames(meta)

ann_colors <- list(
  Condition = c("long covid" = "#FFFACD", "non-long covid" = "#C1E1C1"),
  Sex = c("F" = "#FF66C4", "M" = "#6A5ACD")
)

# Color palette (blue = low, red = high)
library(RColorBrewer)
heat_colors <- colorRampPalette(rev(brewer.pal(n = 9, "RdBu")))(255)


###############################################################
# Select Top 40 Genes with lowest p-value
###############################################################
top_male <- rownames(res_male)[order(res_male$padj)][1:40]
top_male <- top_male[!is.na(top_male)]

top_female <- rownames(res_female)[order(res_female$padj)][1:40]
top_female <- top_female[!is.na(top_female)]


###############################################################
# 12B. Heatmap: Male ONLY
###############################################################
pheatmap(
  norm[top_male, rownames(meta)[meta$Gender == "M"]],
  scale = "row",
  annotation_col = annotation_col[meta$Gender == "M", ],
  annotation_colors = ann_colors,
  color = heat_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize = 12,
  fontsize_row = 8,
  fontsize_col = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  main = "Heatmap: 40 Genes with the Lowest Adjusted p-values (Male)"
)


###############################################################
# 12C. Heatmap: Female ONLY
###############################################################
pheatmap(
  norm[top_female, rownames(meta)[meta$Gender == "F"]],
  scale = "row",
  annotation_col = annotation_col[meta$Gender == "F", ],
  annotation_colors = ann_colors,
  color = heat_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize = 12,
  fontsize_row = 8,
  fontsize_col = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  main = "Heatmap: 40 Genes with the Lowest Adjusted p-values (Female)"
)
