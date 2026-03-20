#'#############################  ---- Section <1> ----  ###############################
#' Core objective: Training the final model in ClassifyModel and applying it to predict new datasets
#' Created on: September 28, 2021
#'##############################################################################
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/0.Build_ClassifyModel_Object.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/1.Extract_Data_and_Labels_from_Model.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/2.Extract_Model_Functions.R") )
source( file.path(Device.path,"2.TransferLearning/svmTL/0.Code/main/SVMTL_final_model_construction.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/Feature_Selection_Framework.R") )
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/Quantile_Normalization.R") )

#'-Part [1]-
#' Content: Construct the final model
#'#######
ConstructFinalModel <- function(cm=NULL, way="mut-exprs", background=FALSE){
	
	if(cm$classifier != "SVMTL"){ #non-transfer learning model
		
		ClassifierXXX = ExtractClassifyFunction(cm$classifier)
		
		tmp = ExtractDataLabel.ClassifyModel(cm, way=way, background=background)
		data = tmp$data
		label = tmp$label
		
		#' Construct training data for the final model
		train.X = data
		train.Y = label
		
		#' Point [1] Train the model and update the model object
		# 1) Feature selection
		method.select.feature = cm$method.select.feature
		cm$final.model$features = FeatureSelect(train.X,train.Y,method.select.feature, topn=cm$optimal.num.features)
		
		features = cm$final.model$features
		
		rf.classifier <- ClassifierXXX(data=train.X[features,], label=train.Y, model=NULL, operate=c("construct"));
		cm$final.model$model = rf.classifier	#update the model element
		
	}else{ #transfer learning model
		cm = ConstructFinalModel.SVMTL(cm)
	}
	return(cm)
}


#'-Part [2]-
#' Content: Use the trained model to predict new data
#'#######
PredictFinalModel <- function(cm, predict.data, normalizeUsingCM=FALSE, #whether to perform quantile normalization on predict.data using training data in the model
		imputeZero=FALSE #whether to fill missing features with zero
){
	
	#' Extract training data
	t.data = ExtractDataLabel.ClassifyModel(cm)$data
	
	# 1) Quantile normalization using reference (target) data, treated as pseudo-bulk
	if(normalizeUsingCM){
		com.genes = intersect( rownames(t.data), rownames(predict.data) )
		predict.data = NormalizeQuantileCrossPlatform(ref.matx=t.data[com.genes,], query.matx=predict.data[com.genes,])
		#Q: Why align genes for cross-platform normalization? Will this lead to substantial gene loss?
		#   Ans 1) Because the underlying R package requires matrices with identical dimensions (it ranks genes per sample and assigns quantiles; e.g., 10,000 genes correspond to 10,000 quantiles)
		#   Ans 2) In practice, high-quality datasets lose only a small number of genes during intersection (e.g., among ~20,000 protein-coding genes, typically only a few hundred are removed)
	}
	
	
	# 2) Fill missing feature genes with zeros
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
	
	#Label prediction
	pred.label =  rep(0,length(pred.value))
	pred.label[pred.value>=cm$optimal.threshold] = 1
	names(pred.label) = colnames(predict.data)
	
	return(list(score=pred.value, label=pred.label))
}