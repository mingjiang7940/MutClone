#'#############################  ---- Section <1> ----  ###############################
#' Core objective: select models and hyperparameters based on grid-searched confidence boundary results
#' 
#'-Part [1]- Content: filter out hyperparameter settings that violate basic criteria
# 				#1) intermediate region proportion > 25%
#				#2) reversed upper and lower bounds (these settings are too loose and can be replaced by stricter ones)
#' 
#'-Part [2]- Content: rank remaining settings by gene, and compute average ranking of metrics (PPV, NPV, AUROC.CV) for each setting within each gene
#' 
#' 
#'-Part [3]- Content: visualize hyperparameter selection results
#' 
#' 
#' Special note:
#' creation date: February 25, 2025
#'##############################################################################
rm( list= ls()  )
Device.path = "D:/Project/0.MutClone/Function"
out.dir = "/WorkSpace/chengmingjiang/1.Project_sc2mutTL/3.Section/00.完成故事梳理后_整理版/S3_模型偶联Atlas/0.Data/新版本/05.模型预测外部单细胞独立数据集"
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
cancer = "CRC"  		# which cancer to select models from #, "LUAD", "BRCA", "CRC"


#'-Part [0]-
#' Content: load confidence boundaries computed under grid parameters
#'#######
library(tidyr)
grid <- expand_grid(
		x = seq(0.90, 0.99, by = 0.01),
		y = seq(0.90, 0.99, by = 0.01)
)
grid = as.data.frame(grid)

all.df = data.frame()
for(i in 1:nrow(grid) ){
	x = grid[i,"x"]
	y = grid[i,"y"]
	#load confidence boundaries
	threshold_down_and_up.list = readRDS(file.path(dirname(out.dir),"04.模型置信边界确定", paste0("threshold_down_and_up.list[",x,"_",y,"].rds")  ))
	df = as.data.frame( Reduce(rbind ,threshold_down_and_up.list) , row.names=paste0(names(threshold_down_and_up.list) ,"|",x,"_",y) )
	df$Scene.id = names(threshold_down_and_up.list)
	df$Algorithm <- sub("^.*\\|", "", df$Scene.id)
	df$Algorithm <- sub("Atlas$", "", df$Algorithm)
	df$Prediction <- sub("\\|.*$", "", df$Scene.id)
	df$param.x = x
	df$param.y = y
	
	all.df = rbind(all.df, df)
}


#Point load cross-validation AUROC of bulk models and map to corresponding models
cv.proformance = readRDS("/WorkSpace/chengmingjiang/1.Project_sc2mutTL/3.Section/00.完成故事梳理后_整理版/S2_迁移预测CRC突变/0.Data/2.迁移模型/CV_proformance_df_AUROC.RDS")
for(i in 1:nrow(all.df)){
	mod = all.df[i,"Prediction"]
	alg = all.df[i,"Algorithm"]
	all.df[i,"AUROC.CV"] = cv.proformance[mod,alg]
}


#'-Part [1]-
#' Content: filter hyperparameter settings that violate basic criteria
#'#######
library(dplyr)
#1) filter settings with intermediate region proportion > 25%
all.df =  all.df %>% filter( Intermediate.porp<0.25 )

#2) filter reversed upper and lower bounds; these settings are too loose and stricter ones are preferred
all.df =  all.df %>% filter( threshold.down<=threshold.up )


#'-Part [2]-
#' Content: rank remaining settings by gene and compute average ranking of metrics (PPV, NPV, ACC, 1-Intermediate.porp)
#'#######
#(2.1) rank settings within each gene for each metric
all.df$classified.porp = 1-all.df$Intermediate.porp
metrics <- c("PPV", "NPV", "ACC", "classified.porp","AUROC.CV")  #"AUROC.full"
rank.df <- all.df %>%
		group_by(Prediction) %>%   # within each model
		mutate(
				across(
						all_of(metrics),
						~ rank(-., ties.method = "average"),  # standard ranking: larger values get better ranks
						.names = "rank.{.col}"
				)
		) %>%
		as.data.frame()


#(2.2) compute composite score using the inverse of mean rank (better rank → larger value)
rank.df2 <- rank.df %>%
		mutate(
				mean.rank = rowMeans(
						cbind(rank.PPV, rank.NPV, rank.AUROC.CV),
						na.rm = TRUE
				)
		)


#(2.3) select the best model and hyperparameters for each gene
best.setting <- rank.df2 %>%
		group_by(Prediction) %>%
		slice_min(mean.rank, n = 1, with_ties = FALSE) %>%  # smallest rank
		ungroup() %>%
		as.data.frame()

#save selected model information
saveRDS(best.setting, file.path(dirname(out.dir), "select.model.RDS") )



#'-Part [3]-
#' Content: visualize hyperparameter selection results
#'#######
#plot overall performance
figAera <- function(df.sub, main=NULL, metric.x="param.x", metric.y="param.y"){
	
	# find the point with the largest size (smallest mean.rank)
	df.star <- df.sub %>%
			filter(mean.rank == min(mean.rank, na.rm = TRUE))
	
	p <- ggplot(
					df.sub,
					aes(
							x = .data[[metric.x]],
							y = .data[[metric.y]],
							color = Algorithm,
							size = 1/mean.rank   # 👈 smaller rank → larger size
					)
			) +
			geom_point( shape = 21,        # hollow circle (can control fill and color)
					fill  = NA,        # no fill → hollow
					stroke = 1.1,        # border width
			) +
			# ⭐ mark the optimal point with a star (separate layer)
			geom_point(
					data = df.star,
					aes(
							x = .data[[metric.x]],
							y = .data[[metric.y]],
							color = Algorithm
					),
					shape = 8,     # star
					size  = 5,
					stroke = 1.2
			) +
			scale_size_continuous( name = "1/(mean rank)" ) +
			labs(
					title = main,
					x = metric.x,
					y = metric.y
			) +
			theme_classic(base_size = 20) +
			theme(
					axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
			)
	
	p
}


#plot
p.list = list()
for(mod.id in unique(rank.df2$Prediction)) {
	library(ggplot2)
	df.sub <- rank.df2 %>%
			filter(Prediction == mod.id)
	p = figAera(df.sub, main=mod.id)
	p.list[[mod.id]] <- p
}


library(patchwork)
pdf(file.path(out.dir, "Results", "1.总体效能_模型选择", "综合排秩_及模型选择情况.pdf"), width = 8 * 7, height = 6)
wrap_plots(p.list, nrow = 1)
dev.off()

