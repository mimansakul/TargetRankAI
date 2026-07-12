# Load annotated results
res <- read.csv("results/deseq2/deseq2_results_annotated.csv")

# Remove genes without symbols
res <- subset(res, !is.na(SYMBOL))

# Remove genes without adjusted p-values
res <- subset(res, !is.na(padj))

# Order by significance
res <- res[order(res$padj), ]

write.csv(
    res,
    "results/deseq2/deseq2_results_clean.csv",
    row.names = FALSE
)

# Significant genes
sig <- subset(
    res,
    padj < 0.05 &
    abs(log2FoldChange) > 1
)

write.csv(
    sig,
    "results/deseq2/significant_genes_annotated.csv",
    row.names = FALSE
)

cat("Clean genes:", nrow(res), "\n")
cat("Significant annotated genes:", nrow(sig), "\n")
