source( file.path(Device.path,"1.TraditionalClassifier/0.Code/0.Build_ClassifyModel_Object.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/1.Extract_Data_and_Labels_from_Model.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/2.Extract_Model_Functions.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/ROC_PR_with_AUC_CI.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/Feature_Selection_Framework.R") )
#' Wrapped cross-validation procedure, currently supporting binary classification problems
#' @param featureData Feature matrix, which must contain one column representing class labels; this column should be named according to the parameter `label`, where 1 indicates positive and 0 indicates negative
#' @param label The name of the column in the feature matrix that represents class labels
#' @param classifier Classifier to use; currently supports SVM, RVM, NaiveBayesian (NA values are not allowed in the feature matrix), logistic, and RandomForest
#' @param cross.repeat Default = 10; the value must be less than min(num(case), num(control))
#' @param time Number of repeated k-fold cross-validations
#' @returnType 
#' @return 
#' @author wintarcy
#' @export
ClassifierCV.FeatureSelection <- function(cm=NULL, 
		data=NULL, #rows are features, columns are samples
		label=NULL, #label column: 1 indicates positive, 0 indicates negative
		data.name=NULL, data.path.rds=NULL,
		
		way="mut-exprs", background=FALSE, #the stored data include information on background genes; whether to use them
		classifier=c("SVM", "RVM", "NaiveBayesian",  "logistic", "RandomForest","linearSVM.svmTL")[1], 
		cross.repeat=5, times=10,
		
		method.select.feature="SVM-RFE", 
		topn=seq(10, 500, 10),
		
		nCore=times
) #used in the feature selection function
{
	if(is.null(cm)){
		cm = CreateClassifyModel(data=data, label=label, data.name=data.name, data.path.rds=data.path.rds,
				classifier=classifier, method.select.feature=method.select.feature)
	}
	
	
	tmp = ExtractDataLabel.ClassifyModel(cm, way=way, background=background)
	data = tmp$data
	label = tmp$label
	
	ClassifierXXX = ExtractClassifyFunction(classifier)
	
	library(parallel)
	t.PP = mclapply(1:times, function(ti){
				
				#' Point [1] Cross-validation partitioning
				library(caret)
				cvIndex.list = createFolds(factor(label), k=cross.repeat, list = TRUE, returnTrain = TRUE)
				
				#' Point [2] Construct the cross-validation component of the classification object
				method.select.feature = cm$method.select.feature 
				cross.validation.models = list()
				for(i in 1:cross.repeat){
					cross.validation.models[[i]] = list(data.index=cvIndex.list[[i]],
							features = rownames(data)
					)
				}
				names(cross.validation.models) = paste0("model.",1:cross.repeat)
				
				
				#'-Part [1]-
				#' Content: Train models within cross-validation
				#'#######
				for( i in 1:cross.repeat){
					
					indx.train = cross.validation.models[[i]]$data.index
					indx.test = setdiff(1:dim(data)[2],indx.train)
					
					train.X = data[,indx.train]
					train.Y = label[indx.train]
					
#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: Perform feature selection and obtain AUC values under different topN settings
#' 
#'############ [24-04-23] ############
					method.select.feature = cm$method.select.feature
					# (1) Rank features according to the selected feature selection method
					if(!is.null(method.select.feature)){ 
						cross.validation.models[[i]]$features = FeatureSelect(train.X, train.Y, method.select.feature, topn=max(topn))
					}else{
						stop("必须提供特征选择方法...")
					}
					
					# (2) Traverse different top-N feature sets and calculate AUC
					cross.validation.models[[i]]$PP = lapply(topn, function(t.topn){
								features = cross.validation.models[[i]]$features[1:t.topn]
								
								# (2) Model construction
								model <- ClassifierXXX(data=train.X[features,], label=train.Y, model=NULL, operate=c("construct"));
								cross.validation.models[[i]]$model = model	#update the model element
								
								# (3) Predict using the model and update performance results in PP
								test.X = data[features,indx.test]
								test.Y	= label[indx.test]
								pred.value <- ClassifierXXX(data=test.X, label=NULL, model=model, operate=c("predict"));
								
								# (5) Construct the PP data object
								tt = list( #
										score = pred.value, #prediction scores obtained from the model on the given data
										label = test.Y #class labels of the corresponding data
								)
								return(tt)
							})
				}
				
				
				# (3) For each topN, merge results across cross-validation folds, calculate the integrated AUC, and select the topN corresponding to the optimal AUC
				temp = lapply(1:length(topn), function(x){
							t.score = lapply(cross.validation.models, function(y){
										return(y$PP[[x]]$score)
									})
							t.label = lapply(cross.validation.models, function(y){
										return(y$PP[[x]]$label)
									})
							AUROC=ROC_PR.normal(score=t.score,			#scores (list)
									label=t.label,				#true labels corresponding to the scores
									type=c("ROC"), 			#ROC or PR
									isCrossValidation=TRUE)		#whether score and label are cross-validation results
							
							return(mean(AUROC))
						}) 
				
				# Obtain the topN with the maximum AUC, and retrieve the corresponding score and label
				max.index = which.max(temp)
				optimal.num.features = topn[max.index]
				
				# (4) Return the cross-validation results under the optimal topN
				t.score = lapply(cross.validation.models, function(y){
							return(y$PP[[max.index]]$score)
						})
				t.label = lapply(cross.validation.models, function(y){
							return(y$PP[[max.index]]$label)
						})
				
				AUROC=ROC_PR.normal(score=t.score,			#scores (list)
						label=t.label,				#true labels corresponding to the scores
						type=c("ROC"), 			#ROC or PR
						isCrossValidation=TRUE)		#whether score and label are cross-validation results
				
				AUPR=ROC_PR.normal(score=t.score,			#scores (list)
						label=t.label,				#true labels corresponding to the scores
						type="PR", 			#ROC or PR
						isCrossValidation=TRUE)		#whether score and label are cross-validation results
				
				PP = list( #
						optimal.num.features=optimal.num.features,
						score=t.score,
						label=t.label,
						AUROC = AUROC,
						AUPR = AUPR)	
				return(PP)
			}, mc.cores=nCore)
	
	# (5) Aggregate results across repeated runs; note that the most frequently selected topN is retained
	tmp = sapply(t.PP, function(x) x$optimal.num.features)
	if(length(tmp) > 1){
		idx = which.max(table(tmp))
		t.max =  as.numeric(names(table(tmp)))[idx] #FindValueWithHighestDensity(t.vector=tmp)
	}else{
		t.max = tmp
	}
	t.PP = t.PP[sapply(t.PP, function(x) x$optimal.num.features == t.max)]
	cm$optimal.num.features = t.max	
	cm$CV.PP = t.PP #there may be multiple cross-validation results
	
	return(cm)
}