#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: extract cancer cells from integrated data as the initial Atlas
#' 
#' Point Input:
#' 				sc.dat				integrated dataset
#' 				meta.data			corresponding meta.data
#' 				pct					corresponding patient center data
#' 
#' Point Output：
#' 				initAtlas.RDS				initial Atlas data
#' 				initAtlas.pseudo_bulk.RDS   pseudo-bulk data after smoothing single cells in the initial Atlas (for downstream model prediction)		
#'############ [25-01-07] ############
rm(list=ls())
Device.path = "D:/Project/0.MutClone/Function"
out.dir = "/WorkSpace/chengmingjiang/1.Project_sc2mutTL/3.Section/00.完成故事梳理后_整理版/S3_模型偶联Atlas/0.Data/新版本"
cancer = "CRC"
if (!dir.exists(out.dir))	dir.create(out.dir, recursive = TRUE)

#'-Part [1]-
#' Content: prepare data
#'#######
#' (1) load integrated data
library(Seurat)
sc.dat = readRDS("/IData3/DataCenter/整合数据/PreMut_CRC_2024/OMICSData/scRNAseq.rds") 
meta.data = readRDS("/IData3/DataCenter/整合数据/PreMut_CRC_2024/OMICSData/meta.data.rds") 
pct = readRDS("/IData3/DataCenter/整合数据/PreMut_CRC_2024/PatientCenter/PatientCenter.rds")


#' 1.2) mutation labels: directly use sample-level mutation labels
#GenomeEvents and SampleInfo
common_cols = setdiff(intersect(colnames(pct$SampleInfo), colnames(pct$GenomeEvents)),c("PatientID", "SampleID"))
idx = !( colnames(pct$GenomeEvents) %in% common_cols)
merged_data <- merge(pct$SampleInfo, pct$GenomeEvents[,idx], by = c("PatientID", "SampleID"), all = TRUE)
#SupplementInfo
common_cols = setdiff(intersect(colnames(merged_data), colnames(pct$SupplementInfo)),c("PatientID", "SampleID"))
idx = !( colnames(pct$SupplementInfo) %in% common_cols)
if(!is.null(pct$SupplementInfo)) merged_data <- merge(merged_data, pct$SupplementInfo[,idx], by = c("PatientID", "SampleID"), all = TRUE )
idx = !is.na(merged_data[,"SampleID"])
merged_data = merged_data[idx,]
#add mutation labels into meta.data
mutation_columns <- grep("Mutation$", colnames(merged_data), value = TRUE)
mutation_data <- merged_data[, c("PatientID", "SampleID", mutation_columns)]
mutation_data <- unique(mutation_data)
rownames(mutation_data) = mutation_data$SampleID
meta.dat_wes = cbind(meta.data, mutation_data[meta.data$SampleID,-(1:2)])

sc.dat@meta.data = meta.dat_wes[colnames(sc.dat),]
sc.dat.cc = subset(sc.dat, CancerCell)

saveRDS(sc.dat.cc, file.path(out.dir,"initAtlas_newCellNames.RDS") )



#'(2) construct pseudo-bulk Point Special note: this step is time-consuming; here we only load the precomputed results
sc.dat.cc0 = sc.dat.cc
DefaultAssay(sc.dat.cc0) <- "RNA"  #ensure the data used are log2(cpm+1); the integrated expression matrix only contains highly variable genes and is not suitable for model prediction
source( file.path(Device.path, "3.BLADE/0.Code/ConstructPseudoBulk_IntraPatient.R") )
seurat.object_pseudo_bulk = psudoBulkInPat(sc.dat.cc0, k=50, mc.cores=1)   #Point input/output units are log2(cpm+1) or log2(tpm+1) 
saveRDS(seurat.object_pseudo_bulk, file.path(out.dir,"initAtlas.pseudo_bulk_newCellNames.RDS") )
#seurat.object_pseudo_bulk = readRDS( file.path(dirname(out.dir),"initAtlas.pseudo_bulk.RDS") )

