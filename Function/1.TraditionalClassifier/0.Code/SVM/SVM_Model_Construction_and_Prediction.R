# TODO: Add comment
# 
# Author: Administrator
###############################################################################

# Traditional SVM model construction and prediction
# The SVM prediction output is a probability value; NA values are not allowed and must be within the range [0, 1]

ClassifierSVM <- function(data=NULL, #feature matrix (rows = features, columns = samples)
		label=NULL, #label column: 1 indicates positive, 0 indicates negative
		model=NULL, operate=c("construct", "predict")[1], ...)
{
	library(e1071)
	if(operate == "construct"){
		featureData = cbind(as.data.frame(t(data)), label=as.factor(label) )
		model <- svm(formula=as.formula(paste("label", " ~ .")), data=featureData, probability = TRUE, ...);
		return(model)
	}
	
	if(operate == "predict"){
		featureData = as.data.frame(t(data))
		temp <- predict(model, as.data.frame(featureData, drop=FALSE), probability = TRUE);
		temp = attr(temp, "probabilities")[,"1"] #extract probability
		return(temp)
	}
}