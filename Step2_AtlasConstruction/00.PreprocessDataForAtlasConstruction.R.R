#'#############################  ---- Section <1> ----  ###############################
#'Core objective： 					  preprocessing before integration
#' 
#' Paragraph <1> dataset quality control
#' Paragraph <2> dataset normalization and merging
#' 
#'##############################################################################
rm(list=ls())
library(Seurat)
out.dir = "/WorkSpace/chengmingjiang/TmpData/integrateData"
#Point Input：
#############
DatasetIDs = c(
		"Joanito_NatureGenetics_2022_CRCSG1",
		"Joanito_NatureGenetics_2022_CRCSG2",
		"Lee_NatureGenetics_2020_KUL3",
		"Joanito_NatureGenetics_2022_KUL5",
		"Lee_NatureGenetics_2020_SMC",
		"Lee_NatureGenetics_2020_SMC2",
		"Qian_CellReshearch_2020",
		"Che_CellDiscovery_2021",
		"Guo_JCIInsight_2022",
		"Li_CancerCell_2023_MSI_CRC",
		"Becker_NatureGenetic_2022",
		"Khaliq_GenomeBiology_2022",
		"DmitrievaPosocco_Nature_2022",
		"Giguelay_Theranostics_2022",
		"Guo_Gastroenterology_2023",
		"Lenos_NatCommun_2022",
		"Pelka_Cell_2021",
		"Poonpanichakul_BiosciRep_2021",
		"Sathe_ClinCancerRes_2023",
		"Sullivan_Gut_2023",
		"Uhlitz_EMBOMolMed_2021" )
names(DatasetIDs) = DatasetIDs


#' Point [0] data import
source("/pub5/xiaoyun/BioU/chengmingjiang/function/0导入数据相关函数[统一标准数据中心]/输出ID.xteam的基础信息.R")
df = retrieve.scRNAseq.ID.xteam(DatasetIDs)
seurat.file.list =  df[,"Path.scRNAseq"]
seurat.list = lapply(seurat.file.list, readRDS)
names(seurat.list) = DatasetIDs
#To prevent duplicate cell names, prepend each cell name with the dataset name
for(t.dat.name in names(seurat.list)){
	new.name = paste(colnames(seurat.list[[t.dat.name]]), t.dat.name,sep="|")
	seurat.list[[t.dat.name]] = RenameCells(seurat.list[[t.dat.name]], new.names =new.name )
}



#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective：       dataset quality control
#' 
#' 
#'-Part [1] Content：for each dataset: cell quality assessment and outlier detection	
#'-Part [2] Content：filter low-depth datasets
#'-Part [3] Content：within each dataset, filter outlier cells and genes
#'-Part [4] Content：filter normal samples and samples with too few cancer cells 				 
#'############ [25-10-27] ############

#'-Part [1]-
#' Content：for each dataset: cell quality assessment
#' #' detect outliers using [Q1 – k(Q3 – Q1),Q3 + k(Q3-Q1)], with k = 1.5	
#'#######
#Point (1.1) extract cell sequencing depth and number of detected genes
library(Seurat)
library(ggplot2)
library(dplyr)
#1.1) extract sequencing depth (nCount_RNA) for each dataset
depth_df.list <- lapply(names(seurat.list), function(name) {
			obj <- seurat.list[[name]]
			data.frame(
					Dataset = name,
					Depth = obj$nCount_RNA,  # or obj@meta.data$nCount_RNA
					Num_Gene = obj$nFeature_RNA,
					pct.mt = obj$percent.mt
			)
		})
depth_df = as.data.frame( Reduce(rbind, depth_df.list) )

#1.2) calculate for each dataset: [Q1 – k(Q3 – Q1),Q3 + k(Q3-Q1)], with k = 1.5, to detect outliers
k = 1.5
limit_df <- depth_df %>%
		group_by(Dataset) %>%
		summarise(
				Q1_depth = quantile(Depth, 0.25, na.rm = TRUE),
				Q3_depth = quantile(Depth, 0.75, na.rm = TRUE),
				lower_depth = Q1_depth - k * (Q3_depth - Q1_depth),
				upper_depth = Q3_depth + k * (Q3_depth - Q1_depth),
				Q1_gene = quantile(Num_Gene, 0.25, na.rm = TRUE),
				Q3_gene = quantile(Num_Gene, 0.75, na.rm = TRUE),
				lower_gene = Q1_gene - k * (Q3_gene - Q1_gene),
				upper_gene = Q3_gene + k * (Q3_gene - Q1_gene)
		)


#Point (1.2) visualize the distribution of cell sequencing depth and gene number for each dataset
limit_df$x <- as.numeric(factor(limit_df$Dataset)) # provide x positions for segments
#1.2.1) visualize Depth
p1 <- ggplot(depth_df, aes(x = Dataset, y = Depth, fill = Dataset)) +
		geom_violin(trim = FALSE, scale = "width") +
		geom_boxplot(width = 0.15, outlier.size = 0.3, color = "black", alpha = 0.5) +
		geom_segment(
				data = limit_df,
				aes(x = x - 0.4, xend = x + 0.4, y = lower_depth, yend = lower_depth),
				color = "blue", linetype = "dashed", size = 0.4
		) +
		geom_segment(
				data = limit_df,
				aes(x = x - 0.4, xend = x + 0.4, y = upper_depth, yend = upper_depth),
				color = "red", linetype = "dashed", size = 0.4
		) +
		scale_y_log10() +
		theme_bw(base_size = 12) +
		theme(
				axis.text.x = element_text(angle = 90, hjust = 1),
				legend.position = "none"
		) +
		labs(
				x = "Dataset",
				y = "Sequencing depth (nCount_RNA)",
				title = "Sequencing depth distribution with per-dataset thresholds"
		)

#1.2.2) visualize Num_Gene
p2 <- ggplot(depth_df, aes(x = Dataset, y = Num_Gene, fill = Dataset)) +
		geom_violin(trim = FALSE, scale = "width") +
		geom_boxplot(width = 0.15, outlier.size = 0.3, color = "black", alpha = 0.5) +
		geom_segment(
				data = limit_df,
				aes(x = x - 0.4, xend = x + 0.4, y = lower_gene, yend = lower_gene),
				color = "blue", linetype = "dashed", size = 0.4
		) +
		geom_segment(
				data = limit_df,
				aes(x = x - 0.4, xend = x + 0.4, y = upper_gene, yend = upper_gene),
				color = "red", linetype = "dashed", size = 0.4
		) +
		scale_y_log10() +
		theme_bw(base_size = 12) +
		theme(
				axis.text.x = element_text(angle = 90, hjust = 1),
				legend.position = "none"
		) +
		labs(
				x = "Dataset",
				y = "Number of detected genes (nFeature_RNA)",
				title = "Gene count distribution with per-dataset thresholds"
		)

png(file.path(out.dir,"0.数据预处理", "测序深度_带线段阈值.png"), width = 2000, height = 1500, res = 200)
print(p1)
dev.off()

png(file.path(out.dir,"0.数据预处理", "表达的基因数量_带线段阈值.png"), width = 2000, height = 1500, res = 200)
print(p2)
dev.off()



#'-Part [2]-
#' Content：filter low-depth datasets
#' 
#'#######
#' Point (2.1) Special note：based on observation, Che_CellDiscovery_2021 and Becker_NatureGenetic_2022 have relatively low sequencing depth
#' Statistical testing further confirmed that Che_CellDiscovery_2021 and Becker_NatureGenetic_2022 generally have lower sequencing depth than other datasets
library(dplyr)
library(openxlsx)
# target datasets
targets <- c("Che_CellDiscovery_2021", "Becker_NatureGenetic_2022")
other_datasets <- setdiff(unique(depth_df$Dataset), targets)
# pairwise one-sided Wilcoxon test function
pairwise_wilcox <- function(df, value_col, targets, others) {
	expand.grid(target = targets, other = others, stringsAsFactors = FALSE) %>%
			rowwise() %>%
			mutate(
					x = list(df[[value_col]][df$Dataset == target]),
					y = list(df[[value_col]][df$Dataset == other]),
					p_value = wilcox.test(unlist(x), unlist(y), alternative = "less")$p.value,
					direction = paste(target, "<", other)
			) %>%
			ungroup() %>%
			select(target, other, direction, p_value) %>%
			mutate(p_adj = p.adjust(p_value, method = "BH"))
}

depth_results <- pairwise_wilcox(depth_df, "Depth", targets, other_datasets)
gene_results <- pairwise_wilcox(depth_df, "Num_Gene", targets, other_datasets)
results_list <- list(
		Depth = depth_results,
		Num_Gene = gene_results
)
write.xlsx(results_list, file.path(out.dir,"0.数据预处理", "Pairwise_depth_gene.xlsx"))


#' Point (2.2) filter datasets
seurat.list = seurat.list[other_datasets]





#'-Part [3]-
#' Content：within each dataset, filter outlier cells and genes
#'#######
#'Point (3.1) filter cells
# 3.1.1) filter using limit_df based on detected outliers
seurat.list.filtered <- lapply(names(seurat.list), function(name) {
			obj <- seurat.list[[name]]
			lim <- limit_df %>% filter(Dataset == name)
			keep.cells <- which(
					obj$nCount_RNA >= lim$lower_depth & obj$nCount_RNA <= lim$upper_depth &
							obj$nFeature_RNA >= lim$lower_gene & obj$nFeature_RNA <= lim$upper_gene
			)
			obj.filtered <- subset(obj, cells = colnames(obj)[keep.cells])
			cat("Filtered", ncol(obj) - ncol(obj.filtered), "cells\n")
			return(obj.filtered)
		})
names(seurat.list.filtered) = names(seurat.list)
#3.1.2) Special note：doublets and mitochondrial gene proportion were uniformly filtered during data collection, so they are not considered here

#'Point (3.2) gene quality control,
seurat.list.filtered = lapply(seurat.list.filtered,function(SeuratObject){
			min.cells = 10
			t.index = (Matrix::rowSums(SeuratObject@assays$RNA@counts > 0) >= min.cells)
			keepGenes = rownames(SeuratObject)[t.index];
			cat("Filtered", nrow(SeuratObject) - length(keepGenes), "genes\n")
			SeuratObject = subset(SeuratObject, features = keepGenes)
		})


#'-Part [4]-
#' Content：filter normal samples and samples with too few cancer cells
#'#######
min.cc = 150  #minimum number of cancer cells
for(datset in names(seurat.list.filtered)){
	
	sc.dat = seurat.list.filtered[[datset]]
	pct.path = file.path("/IData/DataCenter/ColorectalCancer", datset ,"PatientCenter/PatientCenter.rds")
	pct = readRDS( pct.path )
	#1) normal samples
	rm.sample = na.omit( pct$SampleInfo$SampleID[  pct$SampleInfo$SampleType=="Normal"  ] )
	
	#2) samples with fewer than 150 cancer cells
	cell_stats <- sc.dat@meta.data %>%
			dplyr::group_by(SampleID) %>%
			dplyr::summarise(cancer_cells = sum(CancerCell == TRUE))
	rm.sample <- unique( c(rm.sample, cell_stats$SampleID[cell_stats$cancer_cells < min.cc]) )
	
	#3) filter cells
	keep.cell.idx = !(sc.dat@meta.data$SampleID %in% rm.sample)
	
	cat("Filtered", length(rm.sample), "samples\n")
	cat("Involving", ncol(sc.dat) - sum(keep.cell.idx), "cells\n")
	
	if(sum(keep.cell.idx)!=0){
		sc.dat = subset(sc.dat, cells=rownames(sc.dat@meta.data)[keep.cell.idx] )
		seurat.list.filtered[[datset]] = sc.dat
	}else{
		seurat.list.filtered[[datset]] = NULL
	}
}




#'##########  ---- Paragraph <2> ----  ###########  
#' Core objective：   dataset normalization and merging
#' 
#' Point[1] normalization
#' Point[2] merge cells from all datasets
#' Point[3] merge meta.data and add WES mutation status information
#'############ [25-11-03] ############

#' Point[1] normalization
#convert data to log2(cpm+1) [UMI data] or log2(tpm+1) [non-UMI data]
source("/pub5/xiaoyun/BioU/chengmingjiang/function/单细胞数据处理函数/单细胞数据标化.R")
is.UMI = df[,"is.UMI"]
seurat.list.filtered = seurat2normed(seurat.list.filtered, is.UMI)  #Special note: only UMI data are included here

#' Point[2] merge cells from all datasets
seurat.all.cell = merge(seurat.list.filtered[[1]], seurat.list.filtered[2:length(seurat.list.filtered)])
#seurat.all.cell = readRDS(file.path(out.dir,"0.数据预处理","seurat.all.cell.RDS"))

#' Point[3]：merge meta.data and add WES mutation status information
# (1) import meta.data information
meta.data.file.list = df[,"Path.meta"]
meta.data.list = lapply(meta.data.file.list, readRDS)
names(meta.data.list) = DatasetIDs
#use CellID|Dataset as row names of meta.data to prevent duplicate names
for(t.dat.name in names(meta.data.list)){
	rownames(meta.data.list[[t.dat.name]]) = paste(meta.data.list[[t.dat.name]][,"CellID"],
			meta.data.list[[t.dat.name]][,"Dataset"],
			sep="|" )
}
com.col = Reduce(intersect, lapply(meta.data.list, colnames))
meta.data.list = lapply(meta.data.list, function(x) x[,com.col] )
sc.meta.data = Reduce(rbind, meta.data.list)

# (2) import mutation information from datasets with WES
mut.dataset = c("Lee_NatureGenetics_2020_SMC","Joanito_NatureGenetics_2022_CRCSG1","Joanito_NatureGenetics_2022_CRCSG2",
		"Lee_NatureGenetics_2020_KUL3","Qian_CellReshearch_2020","Khaliq_GenomeBiology_2022","Joanito_NatureGenetics_2022_KUL5", 
		"Uhlitz_EMBOMolMed_2021") #"Zheng_SignalTransductTargetTher_2022"
#extract 0/1 mutation profiles from patientCenter
mut.mtx.list = lapply(mut.dataset, function(dat.name){
			tfile = file.path("/IData/DataCenter/ColorectalCancer",dat.name,"PatientCenter/PatientCenter.rds") 
			X = readRDS(tfile)
			merged_data <- merge(X$SampleInfo, X$GenomeEvents, by = c("PatientID", "SampleID"))
			if(!is.null(X$SupplementInfo)) merged_data <- merge(merged_data, X$SupplementInfo,by = c("PatientID", "SampleID"))
			mutation_columns <- grep("Mutation$", colnames(merged_data), value = TRUE)
			mutation_data <- merged_data[, c("PatientID", "SampleID", mutation_columns)]
		} )
names(mut.mtx.list) = mut.dataset

mut.ls = Reduce(union ,lapply(mut.mtx.list, colnames))  #get all mutations with WES
mut.ls = setdiff(mut.ls,c("PatientID","SampleID"))
for(i in 1:length(mut.mtx.list) ){
	mut.mtx.list[[i]][,"Dataset"] = names(mut.mtx.list)[i]
}

library(dplyr)
samp.mut.mtx = dplyr::bind_rows(mut.mtx.list)    #merge mutation matrices at sample level, combining identical columns

# (3) fill the mutation matrix into meta.data
#construct a data.frame to store mutation information, with rows corresponding to sc.meta.data rows and columns corresponding to each mutation
WES.mut = as.data.frame( matrix(NA, nrow(sc.meta.data), length(mut.ls)), col.names = mut.ls)
##fill mutation information
sampID_dataset = paste(sc.meta.data$SampleID ,sc.meta.data$Dataset)
t.indx = match(sampID_dataset, paste(samp.mut.mtx[,"SampleID"], samp.mut.mtx[,"Dataset"]) )
WES.mut.samp = samp.mut.mtx[t.indx,mut.ls]
colnames(WES.mut.samp) = mut.ls

sc.meta.data_wes = cbind(sc.meta.data, WES.mut.samp)
saveRDS(sc.meta.data_wes, file=file.path(out.dir,"0.数据预处理","sc.meta.data_wes0.RDS") ) #temporary storage

##ensure consistent ordering between the integrated data and meta.data
seurat.all.cell@meta.data = sc.meta.data_wes[colnames(seurat.all.cell),]
saveRDS(seurat.all.cell, file=file.path(out.dir,"0.数据预处理","seurat.all.cell.RDS"))
saveRDS(sc.meta.data_wes[colnames(seurat.all.cell),], file=file.path(out.dir,"0.数据预处理","seurat.all.cell_meta.data_wes.RDS") )

