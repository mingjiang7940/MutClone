#'##########  ---- Paragraph <1> ----  ###########  
#' 核心目标： 特征选择仓库, 只要设定method.select.source, 函数自动完成选源
#' 
#' Point 输入：
#' 				source.data.list 源特征矩阵[基因X样本]
#' 				source.lab.list  源01标签
#' 				target.data      靶特征矩阵[基因X样本]
#' Point 输出：
#' 				source.dte.score.sorting 按可迁移性强弱排序的源
#'############ [23-07-04] ############
SourceSelect <- function(source.data.list=NULL,
						source.lab.list=NULL,
						target.data=NULL,
						method.select.source = c("DTE")[1] ){
	
	#Point [3] 使用融合特征进行源排序
	if(method.select.source=="DTE"){
		source( file.path(Device.path,"2.TransferLearning/svmTL/0.Code/DTE_transferability_between_source_and_target_domains.R") )
		source.dte.score.sorting = names( SortSourceByDTE(source.data.list,source.lab.list,target.data) )
	}
	
	return(source.dte.score.sorting)
	
}


