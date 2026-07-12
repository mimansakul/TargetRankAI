library(DESeq2)
library(apeglm)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

# ----------------------------
# Create output folders
# ----------------------------

dir.create("results/deseq2", recursive=TRUE, showWarnings=FALSE)
dir.create("figures", recursive=TRUE, showWarnings=FALSE)

# ----------------------------
# Read count matrix
# ----------------------------

counts <- read.delim(
"results/counts/gene_counts.txt",
comment.char="#"
)

# Remove annotation columns
countData <- counts[,7:ncol(counts)]

# Clean sample names
colnames(countData) <- sub(
    "results\\.bam\\.(SRR[0-9]+)\\.sorted\\.bam",
    "\\1",
    colnames(countData)
)

rownames(countData) <- counts$Geneid

# ----------------------------
# Read metadata
# ----------------------------

metadata <- read.csv(
"metadata/sample_metadata.csv",
row.names=1
)

metadata$condition <- factor(
    metadata$condition,
    levels = c("Healthy", "Granuloma")
)

# ----------------------------
# Create DESeq object
# ----------------------------

dds <- DESeqDataSetFromMatrix(
countData=countData,
colData=metadata,
design=~condition
)

# Remove genes with almost no counts

dds <- dds[rowSums(counts(dds)) >= 10,]

# ----------------------------
# Differential expression
# ----------------------------

dds <- DESeq(dds)

res <- results(dds)

res <- lfcShrink(
    dds,
    coef = "condition_Granuloma_vs_Healthy",
    type = "apeglm"
)

write.csv(
as.data.frame(res),
"results/deseq2/deseq2_results.csv"
)

# ----------------------------
# Normalized counts
# ----------------------------

norm_counts <- counts(dds, normalized=TRUE)

write.csv(
norm_counts,
"results/deseq2/normalized_counts.csv"
)

# ----------------------------
# PCA Plot
# ----------------------------

vsd <- vst(dds)

pdf("figures/PCA_plot.pdf")

plotPCA(vsd, intgroup="condition")

dev.off()

# ----------------------------
# MA Plot
# ----------------------------

pdf("figures/MA_plot.pdf")

plotMA(res)

dev.off()

# ----------------------------
# Volcano Plot
# ----------------------------

pdf("figures/Volcano_plot.pdf")

EnhancedVolcano(
res,
lab=rownames(res),
x="log2FoldChange",
y="padj",
title="TB vs Control"
)

dev.off()

# ----------------------------
# Heatmap
# ----------------------------

topGenes <- head(
order(res$padj),
50
)

mat <- assay(vsd)[topGenes,]

mat <- assay(vsd)[head(order(res$padj), 50), ]
mat <- mat - rowMeans(mat)

pdf("figures/Heatmap_top50.pdf")

pheatmap(
mat,
scale="row"
)

dev.off()

cat("\nAnalysis Complete!\n")
