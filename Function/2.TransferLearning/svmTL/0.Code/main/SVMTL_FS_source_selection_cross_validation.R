#'#############################  ---- Section <1> ----  ###############################
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/0.Build_ClassifyModel_Object.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/1.Extract_Data_and_Labels_from_Model.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/2.Extract_Model_Functions.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/ROC_PR_with_AUC_CI.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/Feature_Selection_Framework.R") )
source( file.path(Device.path,"2.TransferLearning/svmTL/0.Code/source_domain_selection_system.R") )

#' 封装的交叉证实，目前支持两类问题
#' @param featureData，特征矩阵，其中必须有一列以表示类标签，并以参数label命名，标签列1表示pos，0表示阴性
#' @param label ，代表特征矩阵中表示标签的那一列名字
#' @param classifier ，分类器，目前支持SVM，RVM, NaiveBayesian(特征矩阵运行NA）， logistic, RandomForest
#' @param cross.repeat : 默认10，the value must less than min(num(case), num(control))
#' @param time: 重复多少次的10被交叉验证
#' @returnType 
#' @return 
#' @author wintarcy
#' @export
ClassifierCV.FeatureSelection.SVMTL <- function(cm=NULL, 
		way="mut-exprs", background=FALSE, #存储的数据中包括了background基因这一信息，是否需要使用
		classifier=c("SVMTL")[1], 
		cross.repeat=5, times=3,
		method.select.feature="SVM-RFE", 
		topn.feature=c(30, 50, 100, 200, 300, 400, 500),
		method.select.source="DTE", 
		topn.source=NULL,
		nCore=1
) #用于到特征选择函数
{
	library(future.apply)
	
	plan(multisession, workers=nCore)
	
	suppressMessages(library(zoo))
	tmp = ExtractDataLabel.ClassifyModel(cm, way=way, background=background)
	data = tmp$data
	label = tmp$label
	
	#完成源模型的数据装载
	cm$SourceModels = lapply(cm$SourceModels, function(x){
				tmp = ExtractDataLabel.ClassifyModel(x, way=way, background=background)
				x$data = tmp$data
				x$sample.label = tmp$label
				return(x)
			})   #提取源数据
	
	
	T.data = list(data) #提取靶数据
	T.label = list(label)
	S.data = lapply(cm$SourceModels, function(x) x$data )   #提取源数据
	S.label = lapply(cm$SourceModels, function(x) x$sample.label )
	
	#仅保留源和靶特征矩阵同时存在的 融合特征（先取出来大家都有的基因）
	data.list = c(T.data,S.data)
	label.list = c(T.label,S.label)
	com.genes = Reduce(intersect, lapply(data.list, rownames)) 
	
	#'-Part [1]-
	#' 内容：建立超参数组合
	#'#######
	if(is.null(topn.source)){
		topn.source = 1:length(cm$SourceModels)
	}
	superparam_combinations <- expand.grid(param.feature = topn.feature, param.source = topn.source)
	cat("超参数组合数目: ", nrow(superparam_combinations))
	
	ClassifierXXX = ExtractClassifyFunction(classifier)
	suppressMessages(library(parallel))

	t.PP = lapply(1:times, function(ti){
				#' Point [1] 交叉验证划分
				suppressMessages(library(caret))
				cvIndex.list = createFolds(factor(label), k=cross.repeat, list = TRUE, returnTrain = TRUE)
				
				#' Point [2] 构建分类对象交叉验证部分
				method.select.feature = cm$method.select.feature 
				cross.validation.models = list()
				for(i in 1:cross.repeat){
					cross.validation.models[[i]] = list(data.index=cvIndex.list[[i]],
							features = rownames(data)
					)
				}
				names(cross.validation.models) = paste0("model.",1:cross.repeat)
				
				
				#'-Part [1]-
				#' 内容：训练交叉验证中的模型
				#'#######
				for( i in 1:cross.repeat){
					cat("第", ti, "次重复计算，第", i, "次交叉验证...")
					indx.train = cross.validation.models[[i]]$data.index
					indx.test = setdiff(1:dim(data)[2],indx.train)
					
					train.X = data[,indx.train]
					train.Y = label[indx.train]
					
					#'##########  ---- Paragraph <1> ----  ###########  
					#' 核心目标：根据topN规则，分别得到不同topN下的AUC
					#' 
					#'############ [24-04-23] ############
					#' Point [1]
					#'内容：特征选择
					#method.select.feature = cm$method.select.feature
					#（1）根据所选的特征选择方法对特征进行排序
					if(!is.null(method.select.feature)){ 
						T.features = FeatureSelect(train.X, train.Y, method.select.feature, topn=max(topn.feature))			  #选靶特征
						S.features = lapply(cm$SourceModels,  function(x) x$final.model$features ) #提源特征, 预先选好的
						cross.validation.models[[i]]$features = T.features  #融合特征
					}else{
						stop("A feature selection method must be provided...")
					}
					
					#' Point [1]
					#'内容：源的选择
					#method.select.source = cm$method.select.source
					#（1）根据所选的方法对源进行排序
					if(!is.null(method.select.source)){ 
						cross.validation.models[[i]]$optimal.source.index = SourceSelect(S.data, S.label, train.X, method.select.source)
					}else{
						stop("A source selection method must be provided...")
					}
					
					#（2）遍历不同前N的特征，计算AUC
					cross.validation.models[[i]]$PP = future_lapply(1:nrow(superparam_combinations), function(t.index){
								#.cat(levelN=9)
								t.cm = cm
								#（1.1）特征选择
								features = cross.validation.models[[i]]$features[1:superparam_combinations$param.feature[t.index]]
								use.features = union( unlist(S.features), features)
								features = intersect(use.features,com.genes)
								#（1.2）源选择
								source.index = cross.validation.models[[i]]$optimal.source.index[1:superparam_combinations$param.source[t.index]]
								t.cm$SourceModels = t.cm$SourceModels[source.index]
								
								#（2）模型构建
								model <- ClassifierXXX(data=train.X[features,], label=train.Y, model=t.cm, operate=c("construct"));
								cross.validation.models[[i]]$model = model	#更新model元素
								
								#（3）模型预测获得效能，更新PP
								test.X = data[features,indx.test]
								test.Y	= label[indx.test]
								pred.value <- ClassifierXXX(data=test.X, label=NULL, model=model, operate=c("predict"));
								
								#（5）构建数据对象PP部分
								tt = list( #
										score = pred.value, #基于模型与特定的数据，得到的预测得分
										label = test.Y #所采用的特定数据的类标签
								)
								return(tt)
							}, future.chunk.size=1)
				}
				
				#（3）不同topN下，将不同交叉证实的结果合并，计算整合的AUC，选取最优AUC对应的topN
				temp = lapply(1:nrow(superparam_combinations), function(x){
							t.score = lapply(cross.validation.models, function(y){
										return(y$PP[[x]]$score)
									})
							t.label = lapply(cross.validation.models, function(y){
										return(y$PP[[x]]$label)
									})
							AUROC=ROC_PR.normal(score=t.score,			#得分(list)
									label=t.label,				#与得分对应的真实标签
									type=c("ROC"), 			#可选ROC或PR
									isCrossValidation=TRUE)		#score和label是否是交叉验证结果
							
							return(mean(AUROC))
						}) 
				
				#得到最大AUC的topN，并得到相应的score与label
				max.index = which.max(temp)
				#optimal.num.features = superparam_combinations$param.feature[max.index]
				#optimal.num.sources = superparam_combinations$param.source[max.index]
				
				#（4）返回该最优topN下的交叉证实结果
				t.score = lapply(cross.validation.models, function(y){
							return(y$PP[[max.index]]$score)
						})
				t.label = lapply(cross.validation.models, function(y){
							return(y$PP[[max.index]]$label)
						})
				
				AUROC=ROC_PR.normal(score=t.score,			#得分(list)
						label=t.label,				#与得分对应的真实标签
						type=c("ROC"), 			#可选ROC或PR
						isCrossValidation=TRUE)		#score和label是否是交叉验证结果
				
				AUPR=ROC_PR.normal(score=t.score,			#得分(list)
						label=t.label,				#与得分对应的真实标签
						type="PR", 			#可选ROC或PR
						isCrossValidation=TRUE)		#score和label是否是交叉验证结果
				
				PP = list( #
						#optimal.num.features=optimal.num.features,
						#optimal.num.sources=optimal.num.sources,
						optimal.superparameter=max.index,
						score=t.score,
						label=t.label,
						AUROC = AUROC,
						AUPR = AUPR)	
				return(PP)
			})
	
	#（5）多次重复计算，特别注意：选取出评选最多的topN
	tmp = sapply(t.PP, function(x) x$optimal.superparameter)
	if(length(tmp) > 1){
		idx = which.max(table(tmp))
		t.max =  as.numeric(names(table(tmp)))[idx] #FindValueWithHighestDensity(t.vector=tmp)
	}else{
		t.max = tmp
	}
	
	t.PP = t.PP[ sapply(t.PP, function(x) x$optimal.superparameter == t.max) ]
	
	#cm$optimal.num.features = t.max
	cm$optimal.num.features = superparam_combinations$param.feature[t.max]
	cm$optimal.num.sources = superparam_combinations$param.source[t.max]
	cm$CV.PP = t.PP #有可能有多次的交叉证实结果
	
	return(cm)
}







