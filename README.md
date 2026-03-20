# MutClone

**Inferring mutant clones from scRNA-seq data**

---

## 🧠 Overview

Understanding how somatic mutations shape cellular phenotypes requires identifying *mutant clones*—cell populations sharing the same mutation. However, most single-cell RNA sequencing (scRNA-seq) datasets lack matched genotype information.

**MutClone** addresses this limitation by inferring mutant clones directly from gene expression profiles using a transfer learning framework that integrates bulk and single-cell data.

---

## 🚀 Workflow Overview

MutClone consists of four main steps, corresponding to the project structure:

---

### **Step 1. Model Training (`Step1_ModelTraining`)**

- Train mutation prediction models using **multi-cancer bulk transcriptomic datasets**  
- Learn **mutation-associated transcriptional signals**  

---

### **Step 2. Atlas Construction (`Step2_AtlasConstruction`)**

- Construct a **large-scale single-cell atlas**  

---

### **Step 3. Boundary Learning (`Step3_BoundaryLearning`)**

- Project bulk-derived signals into the single-cell expression space 
- Learn **decision boundaries** in the single-cell atlas  
- Identify stable regions corresponding to mutant and wild-type states  

---

### **Step 4. MutClone Inference (`Step4_MutCloneInference`)**

- Apply selected bulk-trained models to scRNA-seq data  
- Predict **mutation probabilities for each cell**  
- Infer mutation status using **confidence-based decision boundaries**  

---

## 🔧 Supporting Modules

- `Function/` — utility functions used across all steps  
- `0.Data/`   — input datasets and processed data  
- `Results/`  — prediction outputs  

---

## 🎯 Summary

**MutClone integrates Multi-cancer transfer learning with atlas-based boundary optimization to enable accurate and scalable inference of mutant clones directly from scRNA-seq data.**