# MutClone

**Inferring mutant clones from scRNA-seq data**

---

## 🧠 Overview

Understanding how somatic mutations shape cellular phenotypes requires identifying *mutant clones*—cell populations sharing the same mutation. However, most single-cell RNA sequencing (scRNA-seq) datasets lack matched genotype information.

**MutClone** addresses this limitation by inferring mutant clones directly from gene expression profiles using a transfer learning framework.

---

## 🚀 Method Summary

MutClone trained in two stages:

1. **Multi-cancer transfer learning (bulk data)**  
   - Learns mutation-associated transcriptional signals from multi-cancer bulk tumor datasets  

2. **Single-cell inference (scRNA-seq)**  
   - Projects learned signals into a large-scale single-cell mutation atlas  
   - Defines stable decision boundaries for robust mutant clone identification  

---