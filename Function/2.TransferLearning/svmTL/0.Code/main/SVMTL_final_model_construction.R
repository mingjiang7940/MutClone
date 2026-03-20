#'#############################  ---- Section <1> ----  ###############################
#'核心目标：ClassifyModel的使用模型[最终模型]训练，以及使用模型用到新数据集的预测
#'创建日期：2021年9月28日
#'##############################################################################
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/0.Build_ClassifyModel_Object.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/1.Extract_Data_and_Labels_from_Model.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/2.Extract_Model_Functions.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/Feature_Selection_Framework.R") )
source( file.path(Device.path,"2.TransferLearning/svmTL/0.Code/source_domain_selection_system.R") )

#'-Part [1]-
#' 内容：构建最终模型 
#'#######
ConstructFinalModel.SVMTL <- function(cm=NULL, way="mut-exprs", background=FALSE){
	ClassifierXXX = ExtractClassifyFunction(cm$classifier)
	
	tmp = ExtractDataLabel.ClassifyModel(cm, way=way, background=background)
	data = tmp$data
	label = tmp$label
	
	#完成源模型的数据装载
	cm$SourceModels = lapply(cm$SourceModels, function(x){
				tmp = ExtractDataLabel.ClassifyModel(x, way=way, background=background)
				x$data = tmp$data
				x$sample.label = tmp$label
				return(x)
			})  
	#提取源和靶数据
	S.data = lapply(cm$SourceModels, function(x) x$data )   #提取源数据
	S.label = lapply(cm$SourceModels, function(x) x$sample.label )
	T.data = list(data) #提取靶数据

	#取出源和靶特征矩阵同时存在的特征（先取出来大家都有的基因）
	data.list = c(T.data,S.data)
	com.genes = Reduce(intersect, lapply(data.list, rownames)) 
	
	
	#' 构建最终模型数据
	train.X = data
	train.Y = label
		
	#' Point [1] 训练模型,更新model
	#（1.1）特征选择
	method.select.feature = cm$method.select.feature
	#source("特征选择综合系统")
	T.features = FeatureSelect(train.X,train.Y,method.select.feature, topn=cm$optimal.num.features) #选靶特征
	S.features = lapply(cm$SourceModels,  function(x) x$final.model$features ) #提源特征, 预先选好的
	#融合特征
	use.features = union( unlist(S.features), T.features)
	features = intersect(use.features,com.genes)
	cm$final.model$features = features

	#（1.2）源选择
	method.select.source = cm$method.select.source
	#source("/pub5/xiaoyun/BioX/DataScience/机器学习/[深度学习]/Learning迁移学习/[1.源的选择]/选源综合系统.R")
	optimal.source.index = SourceSelect(S.data, S.label, train.X, method.select.source)[1:cm$optimal.num.sources]
	cm$SourceModels = cm$SourceModels[optimal.source.index]
	
	#（1.3）模型构建
	model <- ClassifierXXX(data=train.X[features,], label=train.Y, model=cm, operate=c("construct"));
	cm$final.model$model = model	#更新model元素
	
	#(1.4) 拆卸data和label; 注意
	cm["data"] = list(NULL)
	cm["sample.label"] = list(NULL)
	for(ts in names(cm$SourceModels) ){
		cm$SourceModels[[ts]]["data"] = list(NULL)
		cm$SourceModels[[ts]]["sample.label"] = list(NULL)
	}
	
	return(cm)
}


#'-Part [2]-
#' 内容：使用模型预测新数据
#'#######
PredictFinalModel <- function(cm, predict.data, normalizeUsingCM=FALSE, #是否使用模型中的训练数据对predict.data进行分位数标准化
		imputeZero=FALSE #是否对于缺失特征进行补0缺失值
){

	#' 提取训练数据
	t.data = ExtractDataLabel.ClassifyModel(cm)$data
	
	# 1) 参考靶数据, 伪bulk数据进行分位数标准化
	if(normalizeUsingCM){
		com.genes = intersect( rownames(t.data), rownames(predict.data) )
		predict.data = NormalizeQuantileCrossPlatform(ref.matx=t.data[com.genes,], query.matx=predict.data[com.genes,])
		#Q：跨平台标准化，为什么要统一基因？统一后会导致基因大范围缺失吗？
		#	 Ans 1) 因为这个内置R包 只能接收的两个维度一样的矩阵（它会对样本单独排序基因，然后把样本每个基因当成一个分位数，比如有1万个基因它就会记录万个分位数 ）
		#	 Ans 2) 统一基因过程中，质量好的数据通常丢掉基因比较少(2万个蛋白编码基因，统一后，也就丢几百个)
	}
	

	# 2) 缺失的特征基因用0填充
	if(imputeZero)
	{
		test.X = predict.data
		#features = cm$final.model$features
		features <- cm$final.model$features
		
		test.X2 = matrix(0,nrow = length(features), ncol = ncol(predict.data))
		test.X2 = as.data.frame(test.X2, row.names=features )
		colnames(test.X2) = colnames(test.X)
		com.feature = intersect(features, rownames(test.X))
		test.X2[com.feature,] = test.X[com.feature,]
		predict.data = test.X2
	}

	
	ClassifierXXX = ExtractClassifyFunction(cm$classifier)
	pred.value <- ClassifierXXX(data=predict.data, label=NULL, model=cm$final.model$model, operate=c("predict"));
	pred.value = as.numeric(pred.value)
	names(pred.value) = colnames(predict.data)
	
	#标签预测
	pred.label =  rep(0,length(pred.value))
	pred.label[pred.value>=cm$optimal.threshold] = 1
	names(pred.label) = colnames(predict.data)
	
	return(list(score=pred.value, label=pred.label))
}



