# Supplementary Figure and Table Legends

## Extended Data Figures

### Extended Data Fig. 1 | Atlas scale, composition, and boundary-learning support for the CRC single-cell reference
This figure provides supporting detail for the atlas-calibration layer of MutClone. Panels are intended to summarize the size and cohort composition of the CRC atlas, malignant-cell extraction and integration context, and additional boundary-learning support analyses that are informative for the method framework but not essential to keep in the main text.

### Extended Data Fig. 2 | Dataset-specific spread of independent validation performance
This figure expands the independent-validation summary shown in Fig. 4. It presents dataset-level validation performance across the individual external settings so that between-dataset heterogeneity can be inspected without interrupting the main validation narrative.

### Extended Data Fig. 3 | Coverage trade-offs, unique recovered cells, and transcript-silent mutant-cell detail in benchmarking
This figure extends the benchmark analyses from Fig. 5. Panels are intended to document coverage-accuracy trade-offs across methods, sets of mutant cells uniquely recovered by MutClone, and the lack of target-gene transcript expression in many recovered mutant cells. The figure supports the claim that MutClone adds value specifically in mutation-silent settings.

### Extended Data Fig. 4 | Extended TP53 pathway and control analyses
This figure provides additional TP53-focused analyses beyond the main fidelity panels in Fig. 6. It is intended to include extended pathway-enrichment contrasts, control comparisons, and supporting views of the TP53-associated transcriptional program across samples.

### Extended Data Fig. 5 | Cohort-level consistency of APC stem-compartment enrichment
This figure expands the APC analyses from Fig. 7 by showing cohort-level or sample-level consistency of APC-mutant cell enrichment in stem-like compartments. It retains supporting replication detail while leaving the main text focused on the central lineage interpretation.

### Extended Data Fig. 6 | APC pseudotime, trajectory, proliferation, and differentiation-support detail
This figure provides the secondary APC analyses that underpin the differentiation-blockade interpretation. Panels are intended to include proliferation-score comparisons, additional pseudotime views, trajectory-supporting analyses, and supporting expression signatures or regulators linked to impaired differentiation.

### Extended Data Fig. 7 | KRAS module scores, enrichment detail, and CEACAM-associated support
This figure extends the KRAS program branch in Fig. 8. Panels are intended to show cohort-level module-score comparisons, pathway-enrichment detail, and additional support for the CEACAM-linked inflammatory and EMT-like program observed in KRAS-mutant clones.

### Extended Data Fig. 8 | Heterogeneity and controlled analyses of the KRAS metastasis-related signal
This figure provides the supporting detail for the metastasis-related extension summarized in Fig. 8. It is intended to document the cross-cohort heterogeneity of the score-based signal, including controlled analyses such as the KRAS-only background comparison in the Wang dataset and related sensitivity analyses. The figure should be interpreted as boundary-setting support rather than as closure of a causal metastasis claim.

## Extended Data Tables

### Extended Data Table 1 | Summary of datasets included in the CRC single-cell atlas
This table summarizes the cohorts, sample counts, and atlas-related metadata for the datasets included in construction of the CRC single-cell atlas used for boundary calibration.

### Extended Data Table 2 | Summary of independent validation datasets and label definitions
This table summarizes the datasets used for external validation, the source of mutation truth or proxy labels, the level at which labels were assigned, and the key inclusion rules used in the validation framework.

## Supplementary Figures

### Supplementary Fig. 1 | Additional statistics supporting the direct-calling limitation
This figure provides supplementary analyses for the structural limitations of direct mutation calling, such as per-cell mutation-count sparsity, sample-level mutation capture, and related variant-frequency views. It supports the problem-setting claim in Fig. 1 without expanding the main text beyond the minimum persuasive chain.

### Supplementary Fig. 2 | Baseline multi-cancer mutation predictability before transfer learning
This figure provides background for the transfer-learning framework by showing baseline predictability of driver-gene mutation status across cancer types before transfer-based model construction. It motivates why mutation-associated transcriptional programs can be treated as transferable signals.

## Supplementary Tables

### Supplementary Table 1 | Baseline multi-cancer SVM screen for mutation prediction
This table reports the baseline SVM screening results across cancer-gene tasks and provides the background landscape from which transferable mutation-prediction tasks were selected.

### Supplementary Table 2 | Full comparison of CRC model AUROC and AUPRC across candidate methods
This table provides the full performance comparison for CRC mutation-prediction models, including AUROC and AUPRC metrics across the candidate configurations evaluated during method development.

### Supplementary Table 3 | Cross-validation detail for key CRC driver-gene transfer models
This table reports detailed cross-validation results for the selected key CRC driver-gene models, complementing the compact summary presented in Table 1.

### Supplementary Table 4 | Sample-level summary of TP53-focused analyses
This table summarizes the samples included in the TP53 analyses, including the number of TP53-mutant and TP53-wild-type cells or other relevant per-sample inclusion information.

### Supplementary Table 5 | Cohort and sample summary for KRAS program and metastasis analyses
This table summarizes the cohorts and samples included in the KRAS analyses shared across the program branch and the metastasis-related extension.

### Supplementary Table 6 | Functional enrichment detail for the conserved KRAS-mutant transcriptional program
This table contains the pathway- and process-enrichment results associated with the shared KRAS-mutant upregulated gene program described in Section 7.

### Supplementary Table 7 | Per-cohort comparison of metastasis-related transcriptional scores
This table reports the per-cohort statistical comparison of metastasis-related transcriptional scores between KRAS-mutant and KRAS-wild-type cells.

### Supplementary Table 8 | Per-sample detail of metastasis-related score comparisons
This table reports the sample-level metastasis-related score comparisons for samples containing both KRAS-mutant and KRAS-wild-type cells.

### Supplementary Table 9 | Enrichment of KRAS-mutant cells in the top 5% metastasis-score fraction by contingency analysis
This table reports enrichment testing for KRAS-mutant cells within the top 5% of the metastasis-related score distribution using contingency-based analyses.

### Supplementary Table 10 | Enrichment of KRAS-mutant cells in the top 5% metastasis-score fraction by mixed-effects modeling
This table reports mixed-effects-model results for enrichment of KRAS-mutant cells in the top 5% metastasis-related score group across cohorts or samples.

### Supplementary Table 11 | Association between driver mutation status and metastatic-site burden
This table reports the regression analysis testing the association between mutation status of key driver genes and the number of metastatic sites, including the clinically suggestive KRAS signal summarized in Fig. 8.
