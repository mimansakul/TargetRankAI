library(gprofiler2)

# Read annotated significant genes
genes <- read.csv(
  "results/deseq2/significant_genes_annotated.csv",
  stringsAsFactors = FALSE
)

# Keep only genes with symbols
genes <- genes[!is.na(genes$SYMBOL), ]

gene_list <- unique(genes$SYMBOL)

cat("Genes submitted:", length(gene_list), "\n")

# GO + KEGG enrichment
gost_res <- gost(
  query = gene_list,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG")
)

# Save enrichment table
# Convert list columns to character
result <- gost_res$result

result[] <- lapply(result, function(x) {
  if (is.list(x)) {
    sapply(x, function(y) paste(y, collapse = ";"))
  } else {
    x
  }
})

write.csv(
  result,
  "results/deseq2/GO_KEGG_gProfiler.csv",
  row.names = FALSE
)

# Save Manhattan-style enrichment plot
png(
  "results/deseq2/GO_KEGG_gProfiler.png",
  width = 1800,
  height = 1200,
  res = 180
)

p <- gostplot(gost_res)
print(p)

dev.off()

cat("GO and KEGG enrichment completed successfully!\n")
