#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: Rank feature importance based on the absolute values of feature coefficients from a linear-kernel SVM
#' Feature selection using Support Vector Machine based on Recursive Feature Elimination (SVM-RFE)
#' 
#' Input:
# 			Required arguments:
#' 			train.X  Feature matrix (e.g., gene expression matrix; rows = genes, columns = samples)
#' 			train.Y	 Corresponding sample labels
#' 
# 			Method-specific arguments:
#'			topn     Number of features to select, default = 500
#' 			fs.len   Three available modes:
#' 						1."full"  Remove one feature at a time; not recommended when the number of original features (p) is large
# 					 The following modes are based on the principle that more important features are removed more gradually
#' 						2."half"  Remove features following a geometric sequence (1/2), eliminating half of the features at each step
#' 				 (default) 3."power2"  Retain features following the sequence: (before elimination)p, (after first elimination)2^(log2(p)-1), ..., 2^1, (after final elimination)2^0
#' 
#' Output:
#' 			res$fs.order  Feature importance ranked from highest to lowest
#Example: featureselect.svm.rfe(train.X,train.Y,topn=500)
#'############ [23-05-02] ############
featureselect.svm.rfe <- function(train.X, 
								  train.Y,
								  ...){
								  
	unique.parameters = list(...)
	if(is.null(unique.parameters[["topn"]])) topn=500 else topn = unique.parameters[["topn"]]
	if(is.null(unique.parameters[["fs.len"]])) fs.len="power2" else fs.len = unique.parameters[["fs.len"]]
	library(mt)
	res = fs.rfe(t(train.X),factor(train.Y),fs.len=fs.len)
	
	if( length(res$fs.order)>topn ){
		return(res$fs.order[1:topn])
	}else{
		return(res$fs.order)
	}
}




#'##########  ---- Section <1> ----  ###########  
#' Core objective: Perform feature selection using cross-validation
#' 
#' Input:
#' 				ClassifyModel  Classification object after data partitioning
#' 				topn  Number of top features to select
#' 				
#' 				[Optional] from.features  Subset of features to perform selection on; default uses all genes in the expression matrix
#' Output:
#' 				ClassifyModel  Trained model with selected features
#'############ [22-08-23] ############
svmRFE.ClassifyModel <- function(ClassifyModel, topn=500, from.features=NULL){
	
	#Perform feature selection for each cross-validation fold
	for( t.fold in 1:length(ClassifyModel$cross.validation.models) ){
		cat("feature selection fold ",t.fold,"/",length(ClassifyModel$cross.validation.models)," \n")
		
		indx.test = ClassifyModel$cross.validation.models[[t.fold]]$data.index
		indx.train = setdiff(1:dim(ClassifyModel$data)[2],indx.test)
		if(!is.null(from.features)){
			train.X = ClassifyModel$data[from.features,indx.train]
		}else{
			train.X = ClassifyModel$data[,indx.train]
		}	
		train.Y	= ClassifyModel$sample.label[indx.train]			
		
		features.sort <- featureselect.svm.rfe(train.X, train.Y, topn)
		ClassifyModel$cross.validation.models[[t.fold]]$features = features.sort
	}
	
	return(ClassifyModel)
}

