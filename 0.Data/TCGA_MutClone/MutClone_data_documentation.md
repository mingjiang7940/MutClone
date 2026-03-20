# 📦 MutClone Data Documentation

This document describes the directory organization and internal data structures used in the MutClone project. It is intended for reproducibility and ease of use across both R and Python pipelines.

---

# 📁 Directory Structure

```bash
0.MutClone/
├── 0.Data/
│   ├── Results/                                   # Stores model outputs (predictions, etc.)
│   │
│   └── TCGA_MutClone/                             # TCGA-derived mutation-labeled expression datasets
│       ├── CRC_TCGA/                              # Colorectal cancer (CRC)
│       │   ├── mut.exprs.data.indel.RDS           # log2(TPM+1); used for single-cancer SVM models
│       │   ├── mut.pancancer.exprs.data.indel.pkl # Upper-quantile normalized expression; used for dnnTL (Python)
│       │   └── mut.pancancer.exprs.data.indel.RDS # Upper-quantile normalized expression; used for svmTL (R)
│       │
│       ├── BRCA_TCGA/                             # Breast cancer (BRCA)
│       │   ├── mut.exprs.data.indel.RDS
│       │   ├── mut.pancancer.exprs.data.indel.pkl
│       │   └── mut.pancancer.exprs.data.indel.RDS
│       │
│       └── LUAD_TCGA/                             # Lung adenocarcinoma (LUAD)
│           ├── mut.exprs.data.indel.RDS
│           ├── mut.pancancer.exprs.data.indel.pkl
│           └── mut.pancancer.exprs.data.indel.RDS
```

---

## 🔬 Data Processing

- **Single-cancer datasets** (`mut.exprs.data.indel.RDS`) are normalized using **log2(TPM + 1)** and used for **within-cancer SVM models**.
- **Pan-cancer datasets** (`mut.pancancer.exprs.data.indel.*`) are processed using **upper-quantile normalization** (typically the 75th percentile per sample), enabling cross-sample and cross-cancer comparability.
- `.RDS` files are used in **R-based pipelines (svmTL)**.
- `.pkl` files are used in **Python-based pipelines (dnnTL)**.

---

# 🧬 File Data Structure

## `mut.pancancer.exprs.data.indel.RDS`

Each file contains an R list object (`r.object`) with two components:

```r
r.object
├── $all.data : matrix [15042 × 311]
│   ├── rownames: genes (e.g., GTPBP6, A1BG, ...)
│   └── colnames: samples (TCGA IDs)
│
└── $all.label : list (355 genes)
    ├── $APC
    │   ├── $label : integer vector [311]
    │   │   ├── TCGA-3L-AA1B-01 → 1
    │   │   ├── TCGA-4N-A93T-01 → 0
    │   │   └── ...
    │   │
    │   └── $background : character vector [3000]
    │       ├── "CCL25"
    │       ├── "HOXC10"
    │       └── ...
    │
    ├── $TP53
    ├── $KRAS
    └── ...
```

---

## 📊 Component Description

### 🔹 `all.data`

- Type: numeric matrix  
- Dimension: **15,042 genes × 311 samples**  
- Description:
  - Rows represent genes  
  - Columns represent samples (TCGA IDs)  
  - Values represent normalized gene expression  

---

### 🔹 `all.label`

- Type: named list (length = 355 genes)  
- Description: mutation annotations defined per gene  

Each gene contains:

#### ▪ `label`
- Type: named integer vector (length = 311)  
- Description:
  - Mutation status per sample  
  - `1` = mutant  
  - `0` = wild-type  
  - Names must match sample IDs in `all.data`  

#### ▪ `background`
- Type: character vector (~3000 genes)  
- Description:
  - Predefined background gene set   

---

## 📌 Key Properties

- Sample IDs are **consistent across all components**  
- Each gene defines an **independent supervised learning task**  
- Supports both:
  - Single-cancer modeling  
  - Multi-cancer transfer learning (svmTL / dnnTL)  

---

## 📎 Example (R output)

```r
List of 2
 $ all.data : num [1:15042, 1:311] ...
 $ all.label:List of 355
  ..$ ABCA13  :List of 2
  .. ..$ label     : Named int [1:311] ...
  .. ..$ background: chr [1:3000] ...
```

