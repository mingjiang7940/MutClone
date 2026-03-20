#'#############################  ---- Section <1> ----  ###############################
source( file.path(Device.path,"1.TraditionalClassifier/0.Code/Optimize_Threshold_by_Performance_Grid.R") )

ClassifierCV.OptimizeThreshold <- function(cm=NULL, 
		method.threshold="Fscore" #"ROC", #"Fscore","SensitivityAddSpecificity","SensitivityDiffSpecificity"
) #used in the feature selection function
{
	
	#'##########  ---- Paragraph <1> ----  ###########  
	#' Core objective: Determine the optimal classification threshold
	#' 
	#'############ [24-04-06] ############
	score = lapply(cm$CV.PP, function(x) x$score)
	label = lapply(cm$CV.PP, function(x) x$label)
	
	optimal.threshold.ls = lapply(1:length(score), function(i){
				OptimalThresholdClassifyModel(score=score[[i]], label=label[[i]], certification=method.threshold)
			} )
	
	cm$optimal.threshold = mean(unlist(optimal.threshold.ls), na.rm=TRUE)
	cm$method.threshold = method.threshold
	return(cm)
}