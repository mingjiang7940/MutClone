# Figure and Table Legends

## Main Figures

### Fig. 1 | Direct mutation calling is structurally insufficient for clone labeling in mutation-sparse scRNA-seq
This figure defines the problem setting of the study. Panels summarize the limited overlap between exonic mutations identified by genomic sequencing and those recoverable from pseudo-bulk scRNA-seq in paired tumor datasets, the sparse per-cell variant-supporting read and UMI coverage at detected sites, the broad exon-level capture bias across mutated genes, and the frequent monoallelic expression that further reduces recovery of mutant reads even when the gene is expressed. Together, these analyses show that routine scRNA-seq is structurally underpowered for general cell-level mutant-clone labeling.

### Fig. 2 | Transfer learning recovers mutation-associated transcriptional programs across cancers
This figure introduces the first layer of MutClone. Panels outline the transfer-learning workflow, compare baseline SVM models with the transfer-learning models across colorectal cancer driver-gene tasks, illustrate source-domain selection in a representative CRC_FBXW7 example, and show that performance gains increase when poorly matched source domains are excluded. The figure supports the claim that mutation-associated transcriptional programs learned from multi-cancer bulk datasets are transferable and improve prediction performance in the CRC target domain.

### Fig. 3 | Atlas-calibrated confidence boundaries convert transferred mutation programs into single-cell clone inference
This figure introduces the second layer of MutClone. Panels summarize construction of the CRC single-cell atlas, identification of high-confidence mutant and wild-type cell populations, learning of stable dual-threshold decision regions, and a representative calibrated boundary for a CRC_KRAS model. The figure establishes that transferred mutation programs are not used as bulk-style cutoffs, but are constrained by atlas-derived confidence boundaries that preserve an explicit Intermediate zone.

### Fig. 4 | Atlas-calibrated boundaries support robust mutant-clone inference in independent validation settings
This figure presents the independent-validation layer of the manuscript. One panel summarizes the validation framework spanning multiple external scenarios with genomic or same-cell mutation truth. Additional panels show model-level accuracy, precision, negative predictive value, and Intermediate-cell fractions across validation tasks, together with a direct comparison between atlas-calibrated boundary classification and thresholds derived only from bulk cross-validation. The figure supports the conclusion that boundary learning adds practical cell-level value beyond bulk transfer alone.

### Fig. 5 | MutClone improves accuracy-coverage balance and recovers transcriptionally silent mutant cells missed by competing methods
This figure provides the main benchmarking evidence. Panels compare MutClone with Scissor, scAB, and PIPET across datasets and genes, highlight APC-focused evaluations with detailed classification metrics, and quantify mutant-cell recovery across methods, including cells uniquely recovered by MutClone despite lacking target-gene transcript evidence. The figure positions MutClone as an adoption-level advance by combining strong accuracy with usable coverage and mutation-silent cell recovery.

### Fig. 6 | TP53-mutant clone labels recapitulate canonical TP53-loss transcriptional programs
This figure serves as a biological fidelity check for inferred clone labels. Panels show pathway enrichment in TP53-mutant versus TP53-wild-type cells, recurrently upregulated genes across samples, FOXM1-centered enrichment of the conserved TP53-mutant program, comparison with broader malignant epithelial states, and ribosomal protein gene-score behavior relative to normal and tumor controls. Together, the analyses show that TP53-mutant clone labels recover expected TP53-loss biology rather than arbitrary expression clusters.

### Fig. 7 | APC-mutant clones accumulate in stem-like compartments and align with differentiation blockade
This figure presents the strongest biology-facing payoff of the study. Panels show enrichment of APC-mutant cells in stem-like compartments across cohorts and in the integrated analysis, increased stem-cell proportion in samples containing APC-mutant stem cells, lower differentiation scores and earlier pseudotime for APC-mutant stem cells, and lineage-trajectory analyses showing depletion of APC-mutant cells from downstream differentiated branches. The figure supports the interpretation that APC-mutant clones align with differentiation blockade rather than a simple proliferation gain.

### Fig. 8 | KRAS-mutant clones define a shared program branch and a qualified metastasis-associated extension
This shared figure carries both the KRAS program branch and the metastasis-related extension. Program-focused panels show the conserved KRAS-mutant transcriptional program across cohorts, including the shared upregulated-gene set, cohort-level functional themes, and representative CEACAM/MAPK-linked features. Metastasis-focused panels compare metastasis-related transcriptional scores between KRAS-mutant and KRAS-wild-type cells across cohorts, enrichment of KRAS-mutant cells in the top tail of the score distribution, and the clinical association between KRAS mutation and metastatic-site burden. The figure should be interpreted conservatively: the KRAS-associated program is consistent across cohorts, whereas the metastasis-related branch remains heterogeneous and clinically suggestive rather than fully closed.

## Main Tables

### Table 1 | Selected transfer-model performance for key colorectal cancer driver-gene tasks
This table summarizes the performance of the selected MutClone transfer models for the key CRC driver genes. It reports the chosen modeling configuration for each task together with cross-validation and external-validation performance, providing the compact quantitative backbone for the method claim introduced in Fig. 2.
