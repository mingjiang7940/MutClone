#'#############################  ---- Section <1> ----  ###############################
#'Core objective:				perform integration and construct the Atlas using the preprocessed data
#' 
#' Paragraph <1> data integration:
#' 
#' 				Point Input:
#' 							seurat.all.cell merged dataset after preprocessing
#' 
#' 				Point Output:
#' 							seurat_obj_int  integrated dataset
#' 
#' 
#' Special note:
#'Created on: November 3, 2025
#'##############################################################################
rm(list=ls())
library(Seurat)
out.dir = "/WorkSpace/chengmingjiang/TmpData/integrateData"
seurat.all.cell = readRDS(file.path(out.dir,"0.数据预处理","seurat.all.cell.RDS"))
meta.data = readRDS( file.path(out.dir,"0.数据预处理","seurat.all.cell_meta.data_wes.RDS") )
seurat.all.cell@meta.data = meta.data




#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective:		data integration
#' 
#'############ [25-10-29] ############

#'-Part [1]-
#' Identify high-quality samples (batches) as references based on sequencing depth and cell-type distribution uniformity
#'#######
#Point (1.1) calculate and summarize sequencing depth and cell-type distribution uniformity scores
library(dplyr)
meta <- seurat.all.cell@meta.data
# Step 1: CellType.Major distribution for each sample
celltype_dist <- meta %>%
		group_by(Dataset, CellType.Major) %>%
		summarise(n = n(), .groups = "drop") %>%
		group_by(Dataset) %>%
		mutate(freq = n / sum(n))
# Step 2: calculate Shannon entropy (to measure cell-type uniformity)
entropy <- celltype_dist %>%
		group_by(Dataset) %>%
		summarise(entropy = -sum(freq * log(freq)))
# Step 3: calculate sequencing depth for each sample (using mean nCount_RNA here)
depth <- meta %>%
		group_by(Dataset) %>%
		summarise(mean_depth = mean(nCount_RNA))
# Step 4: merge the two metrics
sample_stats <- entropy %>%
		left_join(depth, by = "Dataset")
# Step 5: rank by a combined score of “uniformity” and “sequencing depth” (weights can be adjusted)
sample_stats <- sample_stats %>%
		mutate(score = 0.5*scale(entropy) + 0.5*scale(mean_depth)) %>%  # weighted composite metric
		arrange(desc(score))

#Point (1.2) plot the overall scores to inspect characteristics of high-quality samples (select the 6 samples with clearly higher quality than the others)
library(ggplot2)
library(dplyr)
sample_stats <- sample_stats %>%
		arrange(desc(score))  # sort by score (optional)

#'-Part [2]-
#' Content: use high-quality samples as references for Seurat integration
#'#######
#sort by score from high to low
#seurat.obj.list = SplitObject(seurat.all.cell, split.by="SampleID") #sample IDs are uncommon
seurat.obj.list = SplitObject(seurat.all.cell, split.by="Dataset")
seurat.obj.list = seurat.obj.list[sample_stats$Dataset]

Device.path = "D:/Project/0.MutClone/Function"
source(file.path(Device.path,"4.Seurat","SeuratIntegrate.R"))
seurat_obj_int = SeuratIntegrate(seurat.obj.list, mode="reference-based+RPCA", reference=1:6)
saveRDS(seurat_obj_int, file.path(out.dir,"1.整合数据构建Atlas", "seurat_obj_cell_int.RDS"))
