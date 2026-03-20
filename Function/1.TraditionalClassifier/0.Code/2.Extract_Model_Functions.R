# TODO: Add comment
# 
# Author: Administrator
###############################################################################
source(file.path(Device.path, "1.TraditionalClassifier/0.Code","SVM", "SVM_Model_Construction_and_Prediction.R"))
source(file.path(Device.path, "1.TraditionalClassifier/0.Code","SVM","linearSVM_Model_Construction_and_Prediction_from_SVMTL.R"))
source(file.path(Device.path, "2.TransferLearning/svmTL/0.Code/SVMTL_model_construction_and_prediction.R"))


ExtractClassifyFunction <- function(
		classifier=c("SVM", "RVM", "NaiveBayesian",  "logistic", "RandomForest", "linearSVM.svmTL", "SVMTL", "MultiSource_TargetLabel_dnnTL")[1]
){
	
	ClassifierXXX = switch(classifier, 
			SVM=ClassifierSVM,
			linearSVM.svmTL=ClassifierSVM.svmTL,
			SVMTL=ClassifierSVMTL,
			RVM=ClassifierRVM,
			NaiveBayesian=ClassifierNaiveBayesian,
			RandomForest=ClassifierRandomForest,
			logistic=ClassifierLogistic,
			MultiSource_TargetLabel_dnnTL=ClassifierDnnTL)

	return(ClassifierXXX)
}
	