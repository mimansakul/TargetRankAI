#  TargetRankAI

Machine Learning-Assisted Therapeutic Target Prioritization from Human Tuberculosis RNA-seq Data

---

## Overview

TargetRankAI is an end-to-end RNA-seq bioinformatics workflow that integrates differential gene expression analysis, functional enrichment, and supervised machine learning to prioritize candidate therapeutic targets associated with human tuberculosis.

The pipeline begins with raw RNA-seq data from GEO/SRA and performs quality control, preprocessing, genome alignment, read quantification, differential expression analysis, pathway enrichment, and machine-learning-based target ranking.

---

## Features

- RNA-seq quality assessment
- Read trimming using fastp
- Genome alignment using HISAT2
- BAM processing using SAMtools
- Gene quantification using featureCounts
- Differential expression analysis using DESeq2
- GO and KEGG enrichment analysis
- Gene annotation
- Machine learning–based target prioritization
- Automated ranking of candidate therapeutic targets

---

## Dataset

| Item | Value |
|------|-------|
| GEO Dataset | GSE184537 |
| Disease | Human Tuberculosis |
| Platform | Illumina RNA-seq |
| Samples Used | 4 |

---

## Workflow

```
GEO
 ↓
SRA Toolkit
 ↓
FASTQ
 ↓
FastQC
 ↓
MultiQC
 ↓
fastp
 ↓
HISAT2
 ↓
SAMtools
 ↓
featureCounts
 ↓
DESeq2
 ↓
GO / KEGG
 ↓
Machine Learning
 ↓
Target Ranking
```

---

## Tools

- Ubuntu Linux
- Bash
- Python
- R
- HISAT2
- SAMtools
- FastQC
- MultiQC
- fastp
- featureCounts
- DESeq2
- clusterProfiler
- gProfiler2

---

## Results

### Differential Expression

| Metric | Value |
|---------|------:|
| Genes quantified | 40,829 |
| Annotated genes | 17,394 |
| Significant genes | 486 |
| Upregulated | 283 |
| Downregulated | 259 |

---

### Top Candidate Genes

1. SPP1
2. LTF
3. MMP9
4. ADAMDEC1
5. MMP12
6. PLA2G2D
7. CHIT1
8. COL3A1

---

### Machine Learning Models

- Logistic Regression
- Random Forest
- XGBoost
- Gradient Boosting
- Support Vector Machine

Random Forest was selected as the final prioritization model because it achieved the best balance between precision and recall for this imbalanced classification task.

---

## Repository Structure

```
TargetRankAI/
├── data/
├── docs/
├── figures/
├── metadata/
├── notebooks/
├── results/
├── scripts/
├── README.md
├── LICENSE
└── environment.yml
```

---

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/TargetRankAI.git

cd TargetRankAI

conda env create -f environment.yml

conda activate targetrankai
```

---

## Skills Demonstrated

- Linux
- Bash scripting
- RNA-seq analysis
- Differential expression
- Functional enrichment
- Machine learning
- Data visualization
- Reproducible bioinformatics workflows

---

## Citation

If you use this repository, please cite it using the provided `CITATION.cff`.

---

## License

Released under the MIT License.
