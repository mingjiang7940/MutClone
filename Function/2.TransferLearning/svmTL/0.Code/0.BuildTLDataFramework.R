#'##########  ---- Paragraph <1> ----  ###########  
#' 核心目标：输入源模型list创建TransferLearningClassifyModel对象
#' 
#' Point 输入：
#' 
#' 		必须参数:
#' 				s.model.list   源模型list, 用于训练迁移模型
#' 				data  		   靶数据, 行为特征（必须有行名）、列为样		
#' 				label 		   靶标签, 样本类别标签，名字顺序与mat列名顺序一致,1为阳性0为阴性
#' 
#' 		其余均为可选参数:
#' 				method.select.source    选源的方法[目前支持"DTE"]	
#' 				method.select.feature	特征选择方法
#' 
#' 				method.select.feature   特征选择方法
#' 				optimal.num.features	最佳特征数目[最佳超参数], 如果不指定默认500
#' 
#' 				data.name				数据集的名称[标识data.path.rds中指定标签的names, 例如ABCA13等]
#' 				data.path.rds			数据集的rds地址
#						List of 2      【总特征矩阵和标签, data和label从它产生, 它被存储在模型外以降低模型大小】
#						$ all.data :'data.frame':	15165 obs. of  311 variables: 
#								..$ TCGA-3L-AA1B-01: num [1:15165] 10.02 0 4.53 7.79 7.48 ...
#								..$ TCGA-4N-A93T-01: num [1:15165] 10.83 1.29 7.43 6.67 8.54 ...
#								..$ TCGA-4T-AA8H-01: num [1:15165] 9.493 0.792 4.459 7.451 8.213 ...
#						$ all.label:List of 134
#								..$ ABCA13         :List of 2
#									.. ..$ label     : Named int [1:306] 0 0 0 0 0 0 0 0 0 0 ...
#									.. .. ..- attr(*, "names")= chr [1:306] "TCGA-3L-AA1B-01" "TCGA-4N-A93T-01" "TCGA-4T-AA8H-01" "TCGA-5M-AAT4-01" ...
#									.. ..$ background: chr [1:3000] "CCL25" "CLDN18" "AMER3" "CALB1" ...
#								..$ ADAMTS16       :List of 2
#									.. ..$ label     : Named int [1:304] 0 0 0 0 0 0 0 0 0 0 ...
#									.. .. ..- attr(*, "names")= chr [1:304] "TCGA-3L-AA1B-01" "TCGA-4N-A93T-01" "TCGA-4T-AA8H-01" "TCGA-5M-AATE-01" ...
#									.. ..$ background: chr [1:3000] "CCL25" "CLDN18" "AMER3" "CALB1" ...
#' Point 输出:
#' 
#' 				TransferLearningClassifyModel数据对象
#特别注意：该函数完全继承 [统一标准]分类模型 构建函数
#source("/pub5/xiaoyun/BioX/Bioc/0.BioData/[[统一数据类型]]/分类模型/0.初始化")
#'############ [24-04-23] ############
CreateTLClassifyModel <- function(
		s.model.list=NULL, method.select.source=NULL,
		data=NULL, label=NULL,
		data.name=NULL, data.path.rds=NULL, classifier=NULL, 
		method.select.feature=NULL, optimal.num.features=NULL)
{
	
	#'-Part [2]-
	#' 内容：创建初始数据对象
	#'#######
	ClassifyModel = list(
			data = data, #matrix, 行为特征（必须有行名）、列为样本，用于训练模型的矩阵谱数据集，通常也用于产生交叉证实的结果
			sample.label = label, #向量1，0
			data.name=data.name, 
			data.path.rds=data.path.rds,			
			classifier=classifier,	#整个模型使用的是随机森林
			method.select.feature=method.select.feature,
			optimal.num.features=optimal.num.features,
			SourceModels = s.model.list,
			method.select.source = method.select.source
	);
	
	return(ClassifyModel)
}





