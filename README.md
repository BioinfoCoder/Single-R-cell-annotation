# Efficient Single-Cell Annotation with SingleR

SingleR is a powerful tool for annotating single-cell RNA sequencing (scRNA-seq) data by comparing it to reference datasets. It's efficient, adaptable, and highly accurate for cell type discovery.

## What is SingleR?

**SingleR** is an R package used for annotating single-cell RNA sequencing (scRNA-seq) data. It matches gene expression profiles from your single-cell dataset with a reference dataset, assigning cell type labels based on the most similar reference cells. It's an efficient and effective way to discover cell types in large, complex single-cell datasets.

## How does SingleR work?

### Input
- **Test dataset**: Your single-cell RNA-seq data (e.g., `norm_counts`).
- **Reference dataset**: A curated dataset with labeled cell types (e.g., `ref`).

### Process
- **Expression comparison**: SingleR compares the gene expression of each test cell to all cells in the reference dataset.
- **Cell annotation**: The algorithm assigns a cell type label to each test cell based on similarity scores derived from statistical methods like the Wilcoxon test.

## Advantages of SingleR
- **Fast and efficient**: Annotate large single-cell datasets quickly.
- **Versatile**: Works with any reference dataset, including your own or public repositories.
- **Accurate**: Helps identify rare or novel cell types using a robust similarity-based approach.

## Example Code: SingleR Usage

Here's an example of how to use **SingleR** for cell annotation:

```r
# Install and load necessary packages
install.packages("BiocManager")
BiocManager::install("SingleR")

library(SingleR)

# Running SingleR for cell type annotation
ct_ann <- SingleR(test = norm_counts, 
                  ref = ref, 
                  labels = ref$label.main, 
                  de.method = "wilcox")
