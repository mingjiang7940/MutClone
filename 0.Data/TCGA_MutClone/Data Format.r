##Directory Structure
0.MutClone/
├── 0.Data/
│   ├── Results/                                   # Stores model outputs (predictions, evaluation metrics, etc.)
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


##File data structure
#mut.pancancer.exprs.data.indel.RDS
list(r.object)
├── $all.data : matrix [15042 × 311]
│   ├── rownames: genes (e.g. GTPBP6, A1BG, ...)
│   └── colnames: samples (TCGA IDs)
│
└── $all.label : list (355 genes)
    ├── $APC
    │   ├── $label : vector int [311] 
    │   │   ├── TCGA-3L-AA1B-01 → 1
    │   │   ├── TCGA-4N-A93T-01 → 0
    │   │   └── ...
    │   │
    │   └── $background : vector chr [3000]
    │       ├── "CCL25"
    │       ├── "HOXC10"
    │       └── ...
    │
    ├── $TP53
    ├── $KRAS
    └── ...
	
	
Example:
List of 2
 $ all.data : num [1:15042, 1:311] 10.02 0 4.53 7.79 7.48 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:15042] "GTPBP6" "EFCAB12" "A1BG" "A1CF" ...
  .. ..$ : chr [1:311] "TCGA-3L-AA1B-01" "TCGA-4N-A93T-01" "TCGA-4T-AA8H-01" "TCGA-5M-AAT4-01" ...
 $ all.label:List of 355
  ..$ ABCA13  :List of 2
  .. ..$ label     : Named int [1:311] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. ..- attr(*, "names")= chr [1:311] "TCGA-3L-AA1B-01" "TCGA-4N-A93T-01" "TCGA-4T-AA8H-01" "TCGA-5M-AAT4-01" ...
  .. ..$ background: chr [1:3000] "CCL25" "HOXC10" "CRYBA2" "ALPP" ...
  ..$ ABCA2   :List of 2
  .. ..$ label     : Named int [1:311] 0 0 0 0 0 1 0 1 0 0 ...
  .. .. ..- attr(*, "names")= chr [1:311] "TCGA-3L-AA1B-01" "TCGA-4N-A93T-01" "TCGA-4T-AA8H-01" "TCGA-5M-AAT4-01" ...
  .. ..$ background: chr [1:3000] "CCL25" "HOXC10" "CRYBA2" "ALPP" ...
  
  