source( file.path(Device.path,"1.TraditionalClassifier/0.Code","FeatureSelect","SVM-RFE.R") )

FeatureSelect <- function(train.X, train.Y, method.select.feature, ... ){

	if(method.select.feature=="SVM-RFE"){
		features = featureselect.svm.rfe(train.X, train.Y,...)
	}
	
	return(features)
}