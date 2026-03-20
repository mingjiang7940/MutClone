#'##########  ---- Section <1> ----  ###########  
#' Core objective: considering that cells from wild-type samples are mostly wild-type, while mutant samples may contain mixed wild-type cells.
#'                We use the score density of wild-type cells as a reference, and filter out mutant samples whose cell score distributions are similar to them 
#'                (no difference by rank-sum test or JS similarity > 0.75),
#'                retaining samples potentially enriched with mutant cells so that WES sample labels better approximate cell-level labels.
#' 
#' Paragraph <1>  filter mutant samples whose score density is similar to wild-type samples
#' 
#' Paragraph <2>  visualization: plot density distributions after sample filtering
#'############ [25-11-11] ############
rm(list=ls())
DevicePath = 'D:/Project/0.MutClone/Function'
out.dir = "/WorkSpace/chengmingjiang/1.Project_sc2mutTL/3.Section/00.完成故事梳理后_整理版/S3_模型偶联Atlas/0.Data/新版本/02.过滤突变样本使样本标签近似细胞标签"
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
cancer = "CRC" 
mut.ls = c("APC","KRAS","TP53","BRAF","PIK3CA","SMAD4","FBXW7") #input mutation models to be used


#' load data (with WES labels) and prediction scores
iniAtlas = readRDS( file.path(dirname(out.dir),"initAtlas.RDS") )
pred.res.list = readRDS(file.path(dirname(out.dir),"01.能直接用样本标签代替细胞标签吗","pred.res.list_initAtlas.RDS"))


#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: filter mutant samples whose score density is similar to wild-type samples
#' 
#' Criteria: a mutant sample will be filtered if it meets any of the following
#' 		1) JS similarity > 0.75	
#' 		2) no significant difference by rank-sum test		
#'############ [25-11-17] ############

#' Point [0] load single-cell data
df = iniAtlas@meta.data
#balance the number of cells per sample to avoid thresholds being dominated by a few samples
library(dplyr)
sampling.cells = lapply(unique(df$SampleID), function(samp){
			idx = df$SampleID==samp
			set.seed(8766)  # set random seed
			sample(rownames(df)[idx], 150)
		} )
sampling.cells = unlist(sampling.cells)
df_sub <- df[sampling.cells, ]


keep.mut_samp.list = list()
threshold_js = 0.75              #set similarity threshold
for(mid in names(pred.res.list) ){
	
	#separate cells with and without WES labels
	pred.tab = pred.res.list[[mid]]
	mut <- gsub("CRC_|\\|.*", "", mid)
	t.df = df_sub
	t.df[,colnames(pred.tab)] = pred.tab[rownames(t.df),]
	t.df[,"True.lab"] = t.df[,paste0(mut,"Mutation")]
	library(ggplot2)
	df_na  <- t.df[is.na(t.df$True.lab), ]
	df_val <- t.df[!is.na(t.df$True.lab), ]
	
	
	#' Point [1] filter mutant samples similar to wild-type scores (JS similarity >= 0.75)
	library(dplyr)
	df_wt  <- df_val %>% filter(True.lab == 0)         # all wild-type cells
	df_mut <- df_val %>% filter(True.lab == 1)         # all mutant cells
	keep_samples1 <- c()  # initialize
	for(sid in unique(df_mut$SampleID)){
		# get pred.Prob for this mutant sample and wild-type reference
		mut_scores <- df_mut$pred.Prob[df_mut$SampleID == sid]
		wt_scores <- df_wt$pred.Prob
		source( file.path(Device.path, "3.BLADE/0.Code/DensitySimilarity_BetweenScoreVectors.R") )
		sim <- Get_JensenShannon_similarity(mut_scores, wt_scores)
		# keep if similarity < threshold
		if(sim <= threshold_js){
			keep_samples1 <- c(keep_samples1, sid)
		}
	}
	
	#' Point [2] filter mutant samples with no significant difference from wild-type scores (rank-sum test)
	pvals <- df_mut %>%
			group_by(SampleID) %>%
			summarise(
					p = wilcox.test(pred.Prob, df_wt$pred.Prob, alternative = "greater")$p.value
			) %>%
			mutate(
					FDR = p.adjust(p, method = "BH")
			)
	
	# keep samples with significant difference
	keep_samples2 <- pvals %>% filter(FDR <= 0.05) %>% pull(SampleID)
	
	keep_samples =  intersect(keep_samples1, keep_samples2) #keep_samples1
	keep.mut_samp.list[[mid]] = keep_samples
}

saveRDS(keep.mut_samp.list, file.path(out.dir, "keep.mut_samp.list.RDS") )