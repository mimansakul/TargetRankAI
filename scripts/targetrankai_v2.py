import pandas as pd
import numpy as np

# Ranked genes
genes = pd.read_csv("results/TargetRankAI_ranked_genes.csv")

# GO/KEGG enrichment
enrich = pd.read_csv("results/deseq2/GO_KEGG_gProfiler.csv")

# Existing scores
genes["FC_score"] /= genes["FC_score"].max()
genes["P_score"] /= genes["P_score"].max()

# Count how many enriched pathways each gene participates in
pathway_score = {}

for _, row in enrich.iterrows():
    if pd.isna(row["intersection"]):
        continue

    gene_list = str(row["intersection"]).split(",")

    for g in gene_list:
        g = g.strip()
        pathway_score[g] = pathway_score.get(g, 0) + 1

genes["Pathway_score"] = genes["SYMBOL"].map(pathway_score).fillna(0)

if genes["Pathway_score"].max() > 0:
    genes["Pathway_score"] /= genes["Pathway_score"].max()

genes["TargetRankAI_v2"] = (
      0.40*genes["FC_score"]
    + 0.30*genes["P_score"]
    + 0.30*genes["Pathway_score"]
)

genes = genes.sort_values(
    "TargetRankAI_v2",
    ascending=False
)

genes["Rank"] = range(1, len(genes)+1)

genes.to_csv(
    "results/TargetRankAI_v2.csv",
    index=False
)

print(genes[
    ["Rank","SYMBOL","TargetRankAI_v2"]
].head(20))
