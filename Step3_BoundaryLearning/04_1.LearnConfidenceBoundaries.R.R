#'#############################  ---- Section <1> ----  ###############################
#' Core objective: determine model confidence boundaries based on high-confidence cells
#' 
#' Paragraph <1> obtain high-confidence cells retained in the Atlas
#' 
#' Paragraph <2> based on high-confidence cells, determine lower and upper bounds of thresholds, visualize score distributions and boundaries, and compute accuracy-related metrics
#' 
#' 
#' Special note:
#' creation date: February 25, 2025
#'##############################################################################
rm( list= ls()  )
Device.path = "D:/Project/0.MutClone/Function"
out.dir = "/WorkSpace/chengmingjiang/1.Project_sc2mutTL/3.Section/00.完成故事梳理后_整理版/S3_模型偶联Atlas/0.Data/新版本/04.模型置信边界确定"
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
cancer = "CRC"  		# which cancer to select models from #, "LUAD", "BRCA", "CRC"

#load Atlas paths after iterative clustering
atlas.dir = file.path(dirname(out.dir),"03.迭代聚类寻找高置信的细胞簇")
t.file = list.files(atlas.dir, pattern ="^CRC_.*RDS$")              #Atlas file paths


#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: obtain high-confidence cells retained in the Atlas
#' 
#'############ [26-01-17] ############
all.df = list()           #store all cells
keep.df = list()		  #store qualified cells
for(i in 1:ceiling( (length(t.file)/3) ) ){
	#loading all Atlas objects at once is too large; load 3 at a time
	tmp.f = t.file[(3*(i-1)+1) : min((3*i), length(t.file))]
	atlas.list = lapply(file.path(atlas.dir,tmp.f), readRDS)
	names(atlas.list) = gsub("\\.RDS", "", tmp.f)
	
	#' retain only cells from qualified clusters
	for(m.id in names(atlas.list) ){
		#0) merge dimensional reduction and meta.data information, including iterative clustering results
		umap_df <- as.data.frame(atlas.list[[m.id]][["seurat.object"]]@reductions$umap@cell.embeddings)
		atlas.info = atlas.list[[m.id]][["seurat.object"]]@meta.data
		atlas.info[,c("atlas.X","atlas.Y")] = umap_df
		#1) set parameters for qualified clusters: must match those used in iterative clustering
		major.clust.prop = 0.95
		min.clus.count = 20
		#2) retain qualified clusters based on cluster information
		clus.info = atlas.list[[m.id]][["clus.info"]]
		keep.cluster.ls = c()
		for(clust in unique(clus.info$CellCluster) ){
			sub.clus.info = subset(clus.info, CellCluster==clust)
			total.cells = sum(sub.clus.info[,"n"])
			max.prop = max(sub.clus.info[,"prop"])
			
			if(total.cells>=min.clus.count & max.prop>=major.clust.prop){ #sufficient cluster size and high phenotype consistency
				keep.cluster.ls = c(keep.cluster.ls, clust)
			}
		}
		#keep.cluster.ls
		
		#3) cells from retained qualified clusters
		df = atlas.info
		all.df[[m.id]] = df
		idx = df[,"atlas.IterCluster"] %in% keep.cluster.ls
		df = df[idx,]
		if( nrow(df)!=0 ){ keep.df[[m.id]] = df }
	}
	
	cat(i,"\n")
}





#'##########  ---- Paragraph <2> ----  ###########  
#' Core objective: based on high-confidence cells, determine lower and upper threshold bounds, visualize score distributions and boundaries, and compute accuracy-related metrics
#' 
#'############ [26-01-26] ############
#'0) load prediction results on the initial Atlas
t.dir = file.path(dirname(out.dir), "01.能直接用样本标签代替细胞标签吗")
pred.res.list = readRDS(file.path(t.dir,"pred.res.list_initAtlas.RDS"))

#'1) grid settings for sensitivity and specificity
library(tidyr)
grid <- expand_grid(
		x = seq(0.90, 0.99, by = 0.01),
		y = seq(0.90, 0.99, by = 0.01)
)
grid = as.data.frame(grid)

#'2) compute threshold boundaries and accuracy-related metrics for each grid point
for(i in 1:nrow(grid) ){
	x = grid[i,"x"]
	y = grid[i,"y"]
	p.list = list()
	threshold_down_and_up.list = list()
	for(mid in names(keep.df) ){
		
		#2.0) extract data for high-confidence cells
		pred.res = pred.res.list[[gsub("Atlas","",mid)]]
		df = keep.df[[mid]][,c("atlas.CellCluster","True.lab","atlas.IterCluster","atlas.X","atlas.Y")]
		df[,"pred.Prob"] = pred.res[rownames(df),"pred.Prob"]
		df[,"True.lab"] = as.character(df[,"True.lab"])
		df = na.omit(df)
		
		#models requiring log transformation for plotting
		need.log = c("CRC_BRAF|svmAtlas", "CRC_FBXW7|svmAtlas", "CRC_FBXW7|svmTLAtlas",
				"CRC_KRAS|dnnTLAtlas", "CRC_PIK3CA|svmAtlas", "CRC_APC|svmAtlas",
				"CRC_TP53|dnnTLAtlas", "CRC_TP53|svmAtlas", "CRC_TP53|svmTLAtlas")#determined based on density plots
		
		#2.1) compute confidence boundaries and obtain ggplot objects
		source( file.path(Device.path,"3.BLADE/0.Code/RocGetTwoCutoffs.R") )
		res = RocGetTwoCutoffs(score=df$pred.Prob, label=df$True.lab,
				sensitivity=x, specificity=y,
				main=mid,need.log=(mid %in% need.log) )
		p.list[[mid]] = res$p.plot 
		threshold_down_and_up.list[[mid]] = res$threshold_down_and_up
	}
	saveRDS(threshold_down_and_up.list, file.path(out.dir,paste0("threshold_down_and_up.list[",x,"_",y,"].rds") ))
	
	
	#' 3) visualize confidence boundaries and density distributions
	to_move <- c("CRC_SMAD4|svmAtlas", "CRC_SMAD4|svmTLAtlas")  #this gene only has 2 models; move to the end for better display
	new.order = c(setdiff(names(p.list), to_move),to_move)
	p.list = p.list[new.order]
	
	library(ggplot2)
	library(patchwork)
	combined_plot <- wrap_plots(p.list, nrow = 3, ncol = ceiling(length(p.list) / 3), byrow=FALSE)
	path = file.path(out.dir,"Results",paste0("modelscore_density_plots[",x,"_",y,"].pdf"))
	dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
	ggsave(path,combined_plot, width = ceiling(length(p.list) / 3)*6, height = 3*6,limitsize = FALSE)	
	
}

