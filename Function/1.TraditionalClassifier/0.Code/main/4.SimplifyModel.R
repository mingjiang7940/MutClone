#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: Remove training data and related information from the model, while retaining only cross-validation PP and the final model
#' 
#' Input:
#' 				ClassifyModel	A trained model
#' 
#' Output:
#' 				ClassifyModel	A simplified model
#'############ [23-07-05] ############
SimplifiedModel <- function(ClassifyModel){
	# 1) Remove the training matrix
	ClassifyModel$data = NULL
	ClassifyModel$sample.label = NULL
	#ClassifyModel$classifier = NULL
	# 2) Simplify cross-validation results
	for(t.fold in 1:length(ClassifyModel$cross.validation.models) ){
		ClassifyModel$cross.validation.models[[t.fold]]$data.index = NULL
		ClassifyModel$cross.validation.models[[t.fold]]$method.select.feature = NULL
		ClassifyModel$cross.validation.models[[t.fold]]$features = NULL
		ClassifyModel$cross.validation.models[[t.fold]]$model = NULL
		ClassifyModel$cross.validation.models[[t.fold]]$optimal.threshold = NULL
	}
	return(ClassifyModel)
}