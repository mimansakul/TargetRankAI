library(AnnotationDbi)
library(org.Hs.eg.db)

# Read DESeq2 results
res <- read.csv(
    "results/deseq2/deseq2_results.csv",
    row.names = 1
)

# Remove Ensembl version numbers
ensembl_ids <- sub("\\..*$", "", rownames(res))

# Convert to gene symbols
symbols <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
)

# Convert to Entrez IDs (needed later for GO/KEGG)
entrez <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
)

# Add annotation columns
res$ENSEMBL <- ensembl_ids
res$SYMBOL <- symbols
res$ENTREZID <- entrez

# Reorder columns
res <- res[, c(
    "ENSEMBL",
    "SYMBOL",
    "ENTREZID",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "pvalue",
    "padj"
)]

# Save annotated table
write.csv(
    res,
    "results/deseq2/deseq2_results_annotated.csv",
    row.names = FALSE
)

cat("\nAnnotation Complete!\n")
cat("Genes:", nrow(res), "\n")
cat("Mapped Symbols:", sum(!is.na(res$SYMBOL)), "\n")
cat("Mapped Entrez IDs:", sum(!is.na(res$ENTREZID)), "\n")
