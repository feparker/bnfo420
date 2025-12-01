import pandas as pd

#load tsv files
raw_counts_file = "/Users/faithparker/Downloads/raw_counts_NCBI.tsv"
annotation_file = "/Users/faithparker/Downloads/human_annotation_genes.tsv"

#read tsv files
raw_counts = pd.read_csv(raw_counts_file, sep="\t")
annotation = pd.read_csv(annotation_file, sep="\t", low_memory=False)

#verify that gene ID columns exist
required_cols_raw = {"GeneID"}
required_cols_annotation = {"GeneID"}

if not required_cols_raw.issubset(raw_counts.columns):
    raise ValueError(f"Raw counts file is missing required columns: {required_cols_raw}")

if not required_cols_annotation.issubset(annotation.columns):
    raise ValueError(f"Annotation file is missing required columns: {required_cols_annotation}")

#ensuring that gene IDs overlap between the two files
raw_gene_ids = set(raw_counts["GeneID"])
annotation_gene_ids = set(annotation["GeneID"])

matching_gene_ids = raw_gene_ids.intersection(annotation_gene_ids)
missing_in_annotation = raw_gene_ids - annotation_gene_ids
missing_in_raw = annotation_gene_ids - raw_gene_ids

print(f" Total matching genes: {len(matching_gene_ids)}")
print(f"Gene IDs missing in annotation file: {len(missing_in_annotation)}")
print(f"Gene IDs missing in raw file: {len(missing_in_raw)}")

#merging raw counts with annotation file
merged = raw_counts.merge(annotation, on= "GeneID", how= "left")

#saved merged file
merged_output = "annotated_raw_counts.tsv"
merged.to_csv(merged_output, sep="\t", index=False)

print(f"Merged file saves as: {merged_output}")

