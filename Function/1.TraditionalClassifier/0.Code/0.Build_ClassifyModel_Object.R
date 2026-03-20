###Parameter name [data format] description###
#' Create a ClassifyModel object from input feature matrix and label vector
#' mat: rows represent features (must have row names), columns represent samples
#' lab: sample class labels; the order must match the column names of mat; 1 indicates positive and 0 indicates negative
###########################2021年10月8日#######################
CreateClassifyModel <- function(data=NULL, data.name=NULL, data.path.rds=NULL,
		label=NULL, classifier=NULL, method.select.feature=NULL)
{
	
	#'-Part [2]-
	#' Content: Create the initial data object
	#'#######
	ClassifyModel = list(
			data = data, #matrix; rows = features (must have row names), columns = samples; used for model training and typically also for generating cross-validation results
			sample.label = label, #vector: 1 and 0
			data.name=data.name, 
			data.path.rds=data.path.rds,			
			classifier=classifier,	#the classifier used in the model (e.g., random forest)
			method.select.feature=method.select.feature
	);
	
	
	
	
	return(ClassifyModel)
}