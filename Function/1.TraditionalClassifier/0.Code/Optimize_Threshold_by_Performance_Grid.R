#'#############################  ---- Section <1> ----  ###############################
#' Core objective: Obtain the optimal classification threshold using cross-validation data
#' Grid-search the classification threshold, evaluate performance under each threshold in cross-validation, and finally identify the optimal threshold to store in the data object
#' score is a cross-validation list; the optimal threshold is searched within each fold and then averaged
#' label is a cross-validation list; the optimal threshold is searched within each fold and then averaged
#' 
#' Note: the threshold range here is based on probability outputs between 0 and 1
#' Created on: September 28, 2021
#'##############################################################################

source( file.path(Device.path,"1.TraditionalClassifier/0.Code/PerformanceMetrics.R") )

#'-Part [1]-
#' Content: Main function
#'############ [21-09-29] ############
###Parameter name [data format] description###
#' ClassifyModel Classification framework data object
#' certification 【character】Criterion used to optimize the threshold; options include:
#"Fscore","SensitivityAddSpecificity","SensitivityDiffSpecificity", corresponding to maximum F1 score, maximum sensitivity + specificity, and minimum absolute difference between sensitivity and specificity
###########################2021年9月29日#######################
OptimalThresholdClassifyModel <- function(score, label, certification="Fscore"){
	
	#' Search for the optimal threshold using cross-validation evaluation metrics
	if(certification=="Fscore")  result = OptimalThresholdByFscore(score, label) #search for the optimal threshold using Fscore as the evaluation criterion and store it in the data object
	if(certification=="SensitivityAddSpecificity")  result = OptimalThresholdBySensitivityAddSpecificity(score, label) #search for the optimal threshold using Sensitivity+Specificity as the evaluation criterion and store it in the data object
	if(certification=="SensitivityDiffSpecificity")  result = OptimalThresholdBySensitivityDiffSpecificity(score, label) #search for the optimal threshold using |Sensitivity-Specificity| as the evaluation criterion and store it in the data object
	
	
	return(result)
}




#'-Part [2]-
#' Content: Three subfunctions corresponding to three evaluation criteria
#'############ [21-09-29] ############
OptimalThresholdByFscore <- function(score, label){
	#' Point [1] Function: Search for the optimal threshold when Fscore [maximum] is used as the evaluation criterion
	#' Content:
	#' 1) Grid the probability threshold with a step size of 0.01 to generate a series of candidate thresholds
	#' 2) For each candidate threshold, calculate model cross-validation performance metrics under this threshold (sensitivity, specificity, etc.)
	#' 3) Search for the optimal threshold and retain the model containing the optimal threshold as the final model
	############################################
	best.certification.value=NULL
	
	#Set the threshold grid
	min.score <- min(unlist(score), na.rm = TRUE)
	max.score <- max(unlist(score), na.rm = TRUE)
	thresholds <- seq(min.score, max.score, length.out = 101)[-c(1,101)]  #this ensures that if the optimal threshold is smaller than 0.01, it can still be selected; when the minimum and maximum values are 0 and 1, the threshold range is seq(0.01,0.99,0.01)
	
	final.thre = 0
	for( thre in  thresholds){ #' 1) Grid the threshold
		
		#' 2) Calculate performance under the current threshold
		#' cat("Calculating threshold",thre,"\n")
		tmp = sapply(1:length(score), function(i){
					t.performance = PerformanceMetrics(score=unlist(score[[i]]), label=unlist(label[[i]]), threshold=thre, Beta=1) 
					f.score = t.performance["F.measure"]
					return(f.score)
				})
		
		
		#' Extract the criterion value (evaluation metric)
		tmp.certification.value = mean(tmp,na.rm=TRUE)
		
		#' 3) Search for the optimal threshold and retain the model containing the optimal threshold as the final model
		if(thre == thresholds[1]){#if this is the first threshold, use it as the current optimum
			best.certification.value = tmp.certification.value
			final.thre = thre
			
		}else{#if not the first threshold, determine whether the current threshold can replace the historical optimum; if the current metric is larger, record the current threshold
			if(!is.na(best.certification.value<tmp.certification.value) & (best.certification.value<=tmp.certification.value) ){ #
				best.certification.value = tmp.certification.value
				final.thre = thre
			}	
		}
	}
	
	return(final.thre)
}



OptimalThresholdBySensitivityAddSpecificity <- function(score, label){
	#' Point [1] Function: Search for the optimal threshold when Sensitivity+Specificity [maximum] is used as the evaluation criterion
	#' Content:
	#' 1) Grid the probability threshold with a step size of 0.01 to generate a series of candidate thresholds
	#' 2) For each candidate threshold, calculate model cross-validation performance metrics under this threshold (sensitivity, specificity, etc.)
	#' 3) Search for the optimal threshold and retain the model containing the optimal threshold as the final model
	best.certification.value=NULL
	
	#Set the threshold grid
	min.score <- min(unlist(score), na.rm = TRUE)
	max.score <- max(unlist(score), na.rm = TRUE)
	thresholds <- seq(min.score, max.score, length.out = 101)[-c(1,101)]  #this ensures that if the optimal threshold is smaller than 0.01, it can still be selected; when the minimum and maximum values are 0 and 1, the threshold range is seq(0.01,0.99,0.01)
	
	final.thre = 0
	for( thre in  thresholds){ #' 1) Grid the threshold
		
		#' 2) Calculate performance under the current threshold
		#' cat("Calculating threshold",thre,"\n")
		tmp = sapply(1:length(score), function(i){
					t.performance = PerformanceMetrics(score=unlist(score[[i]]), label=unlist(label[[i]]), threshold=thre, Beta=1) 
					sensitivity = t.performance["sensitivity"]
					specificity = t.performance["specificity"]
					return(sensitivity+specificity)
				})
		
		tmp.certification.value = mean(tmp,na.rm=TRUE)
		
		#' 3) Search for the optimal threshold and retain the model containing the optimal threshold as the final model
		if(thre == thresholds[1]){#if this is the first threshold, use it as the current optimum
			best.certification.value = tmp.certification.value
			final.thre = thre
			
		}else{#if not the first threshold, determine whether the current threshold can replace the historical optimum
			if(!is.na(best.certification.value<tmp.certification.value) & best.certification.value<=tmp.certification.value){
				best.certification.value = tmp.certification.value
				final.thre = thre
			}
		}
		
	}
	
	return(final.thre)
}



#'-Part [1]-
#' Content:
#'#######
OptimalThresholdBySensitivityDiffSpecificity <- function(score, label){
	#' Point [1] Function: Search for the optimal threshold when |Sensitivity-Specificity| [minimum] is used as the evaluation criterion
	#' Content:
	#' 1) Grid the probability threshold with a step size of 0.01 to generate a series of candidate thresholds
	#' 2) For each candidate threshold, calculate model cross-validation performance metrics under this threshold (sensitivity, specificity, etc.)
	#' 3) Search for the optimal threshold and retain the model containing the optimal threshold as the final model
	#Set the threshold grid
	min.score <- min(unlist(score), na.rm = TRUE)
	max.score <- max(unlist(score), na.rm = TRUE)
	thresholds <- seq(min.score, max.score, length.out = 101)[-c(1,101)]  #this ensures that if the optimal threshold is smaller than 0.01, it can still be selected; when the minimum and maximum values are 0 and 1, the threshold range is seq(0.01,0.99,0.01)
	
	final.thre = 0
	for( thre in  thresholds){ #' 1) Grid the threshold
		
		#' 2) Calculate performance under the current threshold
		#' cat("Calculating threshold",thre,"\n")
		tmp = sapply(1:length(score), function(i){
					t.performance = PerformanceMetrics(score=unlist(score[[i]]), label=unlist(label[[i]]), threshold=thre, Beta=1) 
					sensitivity = t.performance["sensitivity"]
					specificity = t.performance["specificity"]
					return(abs(sensitivity-specificity))
				})
		
		tmp.certification.value = mean(tmp,na.rm=TRUE)
		
		#' 3) Search for the optimal threshold and retain the model containing the optimal threshold as the final model
		if(thre == thresholds[1]){#if this is the first threshold, use it as the current optimum
			best.certification.value = tmp.certification.value
			final.thre = thre
			
		}else{#if not the first threshold, determine whether the current threshold can replace the historical optimum
			if(!is.na(best.certification.value>=tmp.certification.value) & best.certification.value>=tmp.certification.value){
				best.certification.value = tmp.certification.value 
				final.thre = thre
			}	
		}
		
	}
	
	return(final.thre)
}