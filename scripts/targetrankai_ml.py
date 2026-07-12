import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    roc_auc_score
)
import matplotlib.pyplot as plt

# ----------------------------
# Load ranked genes
# ----------------------------

df = pd.read_csv("results/TargetRankAI_ranked_genes.csv")

# Keep only complete rows
df = df.dropna(subset=[
    "baseMean",
    "log2FoldChange",
    "padj",
    "TargetRankAI_score"
])

# ----------------------------
# Create labels
# Top 20% genes = High priority
# ----------------------------

cutoff = df["TargetRankAI_score"].quantile(0.80)

df["HighPriority"] = (
    df["TargetRankAI_score"] >= cutoff
).astype(int)

# ----------------------------
# Features
# ----------------------------

X = df[[
    "baseMean",
    "log2FoldChange",
    "padj"
]]

X["neglog10padj"] = -np.log10(X["padj"] + 1e-300)
X = X.drop(columns="padj")

y = df["HighPriority"]

# ----------------------------
# Train/Test
# ----------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.25,
    random_state=42,
    stratify=y
)

# ----------------------------
# Random Forest
# ----------------------------

rf = RandomForestClassifier(
    n_estimators=300,
    random_state=42
)

rf.fit(X_train, y_train)

pred = rf.predict(X_test)

prob = rf.predict_proba(X_test)[:,1]

print("\nClassification Report\n")
print(classification_report(y_test,pred))

print("\nConfusion Matrix\n")
print(confusion_matrix(y_test,pred))

print("\nROC AUC:", roc_auc_score(y_test,prob))

# ----------------------------
# Feature Importance
# ----------------------------

importance = pd.DataFrame({
    "Feature":X.columns,
    "Importance":rf.feature_importances_
})

importance = importance.sort_values(
    "Importance",
    ascending=False
)

print("\nFeature Importance\n")
print(importance)

importance.to_csv(
    "results/ml_feature_importance.csv",
    index=False
)

plt.figure(figsize=(7,4))

plt.bar(
    importance["Feature"],
    importance["Importance"]
)

plt.ylabel("Importance")

plt.tight_layout()

plt.savefig(
    "results/ml_feature_importance.png",
    dpi=300
)

print("\nSaved:")
print("results/ml_feature_importance.csv")
print("results/ml_feature_importance.png")
