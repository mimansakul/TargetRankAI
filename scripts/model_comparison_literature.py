import pandas as pd
import numpy as np
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import StratifiedKFold, cross_validate
import time

from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score
)

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC

# Optional XGBoost
try:
    from xgboost import XGBClassifier
    xgb_available = True
except ImportError:
    xgb_available = False

# ----------------------------
# Load ranked genes
# ----------------------------

df = pd.read_csv("results/TargetRankAI_ranked_genes.csv")

# ----------------------------
# Load literature-supported TB genes
# ----------------------------

with open("data/literature/tb_known_genes.txt") as f:
    tb_genes = set(line.strip() for line in f if line.strip())

# ----------------------------
# Create labels
# ----------------------------

df["HighPriority"] = df["SYMBOL"].isin(tb_genes).astype(int)

print("\nKnown TB genes in dataset:", df["HighPriority"].sum())

# ----------------------------
# Features
# ----------------------------

X = df[[
    "baseMean",
    "log2FoldChange",
    "padj"
]].copy()

X["neglog10padj"] = -np.log10(X["padj"] + 1e-300)
X.drop(columns="padj", inplace=True)

y = df["HighPriority"]

# ----------------------------
# Train/Test Split
# ----------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.25,
    random_state=42,
    stratify=y
)

smote = SMOTE(
    random_state=42,
    k_neighbors=2
)

X_train, y_train = smote.fit_resample(
    X_train,
    y_train
)

# ----------------------------
# Models
# ----------------------------

models = {

    "Logistic Regression":
        LogisticRegression(
            max_iter=1000,
            class_weight="balanced"
        ),

    "Random Forest":
        RandomForestClassifier(
            n_estimators=500,
            max_depth=8,
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=42
        ),

    "Gradient Boosting":
        GradientBoostingClassifier(
            random_state=42
        ),

    "SVM":
        SVC(
            probability=True,
            class_weight="balanced",
            random_state=42
        )
}

if xgb_available:

    scale = (len(y_train) - sum(y_train)) / sum(y_train)

    models["XGBoost"] = XGBClassifier(
        n_estimators=300,
        learning_rate=0.05,
        max_depth=4,
        subsample=0.8,
        colsample_bytree=0.8,
        eval_metric="logloss",
        scale_pos_weight=scale,
        random_state=42
    )

# ----------------------------
# Train and Evaluate
# ----------------------------

results = []

best_model = None
best_auc = 0
best_f1 = -1

for name, model in models.items():

    start = time.perf_counter()

    model.fit(X_train, y_train)

    pred = model.predict(X_test)

    prob = model.predict_proba(X_test)[:, 1]

    runtime = round(time.perf_counter() - start, 4)

    auc = roc_auc_score(y_test, prob)

    current_f1 = f1_score(y_test, pred, zero_division=0)

    results.append({
        "Model": name,
        "Accuracy": accuracy_score(y_test, pred),
        "Precision": precision_score(y_test, pred, zero_division=0),
        "Recall": recall_score(y_test, pred, zero_division=0),
        "F1": current_f1,
        "ROC_AUC": auc,
        "Runtime_sec": runtime
    })

    if (
        best_model is None
        or current_f1 > best_f1
        or (current_f1 == best_f1 and auc > best_auc)
    ):
        best_model = model
        best_auc = auc
        best_f1 = current_f1

results = pd.DataFrame(results)

results = results.sort_values(
    ["F1", "ROC_AUC"],
    ascending=[False, False]
)

print("\n============================")
print("MODEL COMPARISON")
print("============================\n")

print(results)

results.to_csv(
    "results/model_comparison.csv",
    index=False
)

# ----------------------------
# Predict all genes
# ----------------------------

df["Prediction_Probability"] = best_model.predict_proba(X)[:,1]

df = df.sort_values(
    "Prediction_Probability",
    ascending=False
)

df.to_csv(
    "results/final_gene_predictions.csv",
    index=False
)

print("\nBest model:", results.iloc[0]["Model"])

print("\nTop predicted genes:")

print(df[
    [
        "SYMBOL",
        "Prediction_Probability",
        "log2FoldChange",
        "padj"
    ]
].head(20))

print("\nSaved files:")
print("results/model_comparison.csv")
print("results/final_gene_predictions.csv")

# ----------------------------
# Feature Importance (Tree models)
# ----------------------------

if hasattr(best_model, "feature_importances_"):

    imp = pd.DataFrame({
        "Feature": X.columns,
        "Importance": best_model.feature_importances_
    })

    imp = imp.sort_values(
        "Importance",
        ascending=False
    )

    imp.to_csv(
        "results/ml_feature_importance_literature.csv",
        index=False
    )

    print("\nFeature Importance")
    print(imp)
