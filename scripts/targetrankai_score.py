import pandas as pd
import numpy as np

# Load annotated DESeq2 results
df = pd.read_csv("results/deseq2/significant_genes_annotated.csv")

# Remove genes without symbols
df = df.dropna(subset=["SYMBOL"])

# Feature 1: absolute fold change
df["FC_score"] = abs(df["log2FoldChange"])

# Feature 2: statistical significance
df["P_score"] = -np.log10(df["padj"] + 1e-300)

# Normalize scores
df["FC_score"] /= df["FC_score"].max()
df["P_score"] /= df["P_score"].max()

# Final TargetRankAI score
df["TargetRankAI_score"] = (
    0.5 * df["FC_score"] +
    0.5 * df["P_score"]
)

# Rank
df = df.sort_values(
    "TargetRankAI_score",
    ascending=False
)

df["Rank"] = range(1, len(df)+1)

# Save
df.to_csv(
    "results/TargetRankAI_ranked_genes.csv",
    index=False
)

print(df[[
    "Rank",
    "SYMBOL",
    "TargetRankAI_score",
    "log2FoldChange",
    "padj"
]].head(20))

print("\nSaved:")
print("results/TargetRankAI_ranked_genes.csv")
