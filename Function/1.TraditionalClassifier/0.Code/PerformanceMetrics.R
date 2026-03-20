# TODO: Calculate evaluation metrics for classification performance 
# 
###############################################################################


#' Compute multiple performance evaluation metrics after classification prediction
#' @param score.label Note: under the unified classifier structure, this should be a list, where each element is a data.frame with two columns: score and label
#'                  score Predicted scores computed by the model, vector
#'                  label True class labels; allowed labels must be numeric 1 or 0
#' @param threshold Samples with probability greater than or equal to this threshold are classified as positive (1); otherwise negative (0)
#' @param Beta Parameter for the “F.measure” statistic. Beta controls the relative emphasis on sensitivity versus precision; when Beta > 1, more emphasis is placed on sensitivity; when Beta < 1, more emphasis is placed on precision; when Beta = 1,
#'             sensitivity and precision are weighted equally. Beta can be set as the ratio of the larger class size to the smaller class size to emphasize sensitivity
#' @returnType Numeric vector
#' @return Values of multiple performance evaluation metrics
#' 
#' @author xiaoyun
#' @export
PerformanceMetrics <- function(score, label, threshold=0.5, Beta=1) 
{
	#score = unlist(sapply(score.label, function(x) x[,1]));
	#label = unlist(sapply(score.label, function(x) x[,2]));
	
	#Determine predicted class labels according to the provided threshold
	predict.label <- as.numeric(score >=  threshold);
	
	TP <- sum((label == predict.label) & (label == 1));
	FP <- sum((label != predict.label) & (label == 0));
	TN <- sum((label == predict.label) & (label == 0));
	FN <- sum((label != predict.label) & (label == 1));		
	
	#sensitivity/Recall/True Positive rate/hit rate
	sensitivity <- TP/(TP+FN);
	recall = sensitivity;
	tpr = sensitivity
	
	#specificity/True negative rate
	specificity <- TN/(TN+FP);
	
	#False Positive rate
	fpr = FP/(TN+FP)
	
	precision <- TP/(TP+FP);
	negative.predictive.value <- TN/(FN+TN);
	accuracy <- (TP+TN)/(TP+FP+TN+FN);
	prevalent = (TP + FN)/(TP + FN + FP + TN) ##add positive class prevalence
	
	#MCC <- (TP*TN-FP*FN)/(sqrt((TP+FN)*(TN+FP)*(TP+FP)*(TN+FN))) #computing (TP+FN)*(TN+FP)*(TP+FP)*(TN+FN) first and then sqrt may cause integer overflow
	MCC <- (TP*TN-FP*FN)/(sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TP+FP)*sqrt(TN+FN))
	
	#For class-imbalance scenarios
	G.mean <- (sensitivity*specificity)^(1/2);
	
	F.measure <- ((1+Beta^2)*sensitivity*precision)/((Beta^2)*precision+sensitivity);
	
	
	
	return(c(TP=TP, FP=FP, TN=TN, FN=FN, 
					sensitivity=sensitivity, 
					specificity=specificity, 
					precision=precision,
					recall=recall,
					tpr=tpr,
					fpr=fpr,
					negative.predictive.value=negative.predictive.value, 
					accuracy=accuracy, 
					prevalent=prevalent,
					MCC=MCC, 
					G.mean=G.mean, 
					F.measure=F.measure))
}

#' Slightly different from PerformanceMetrics: each object in score has a weight, so score must have names, and the weight of each named object can be retrieved from weights
#' Final calculations of TP, FP, and related metrics are obtained by summing weights (PerformanceMetrics is effectively equivalent to assigning weight = 1 to every object)
#' This weighted metric calculation is relatively uncommon; current evidence comes from 【23353650, Nature Methods, 2013】, 【27604469, Genome Biology, 2016】, and 【31744546, Genome Biology, 2019】
#' @param score Predicted scores computed by the model, vector
#' @param label True class labels; allowed labels must be numeric 1 or 0, vector, paired with score
#' @param weights Vector storing the weight of each object in score
#' @param threshold Samples with probability greater than or equal to this threshold are classified as positive
#' @param Beta Parameter for F.measure. Beta controls the relative emphasis on sensitivity versus precision; when Beta > 1, more emphasis is placed on sensitivity; when Beta < 1, more emphasis is placed on precision; when Beta = 1,
#'             sensitivity and precision are weighted equally. Beta can be set as the ratio of the larger class size to the smaller class size to emphasize sensitivity
#' @returnType Numeric vector
#' @return Values of multiple performance evaluation metrics 
#' 
#' @author XTeam 
#' @export
PerformanceMetricsUsingWeight <- function(score, label, weights, threshold=0.5, Beta=1) 
{
	
	#Determine predicted class labels according to the provided threshold
	predict.label <- as.numeric(score >=  threshold);
	
	TP <- sum(weights[names(label[(label == predict.label) & (label == 1)])]); ## TP is the sum of weights for true positive objects
	FP <- sum(weights[names(label[(label != predict.label) & (label == 0)])]); ## FP is the sum of weights for false positive objects
	TN <- sum(weights[names(label[(label == predict.label) & (label == 0)])]); ## TN is the sum of weights for true negative objects
	FN <- sum(weights[names(label[(label != predict.label) & (label == 1)])]); ## FN is the sum of weights for false negative objects
	
	
	#sensitivity/Recall/True Positive rate/hit rate
	sensitivity <- TP/(TP+FN);
	recall = sensitivity;
	tpr = sensitivity
	
	#specificity/True negative rate
	specificity <- TN/(TN+FP);
	
	#False Positive rate
	fpr = FP/(TN+FP)
	
	precision <- TP/(TP+FP);
	negative.predictive.value <- TN/(FN+TN);
	accuracy <- (TP+TN)/(TP+FP+TN+FN);
	
	MCC <- (TP*TN-FP*FN)/(sqrt((TP+FN)*(TN+FP)*(TP+FP)*(TN+FN)))
	
	#For class-imbalance scenarios
	G.mean <- (sensitivity*specificity)^(1/2);
	
	F.measure <- ((1+Beta^2)*sensitivity*precision)/((Beta^2)*precision+sensitivity);
	
	return(c(TP=TP, FP=FP, TN=TN, FN=FN, 
					sensitivity=sensitivity, 
					specificity=specificity, 
					precision=precision,
					recall=recall,
					tpr=tpr,
					fpr=fpr,
					negative.predictive.value=negative.predictive.value, 
					accuracy=accuracy, 
					MCC=MCC, 
					G.mean=G.mean, 
					F.measure=F.measure))
}



#' Given a threshold vector generated by a specified step size, compute the average level of performance metrics across different prediction scores under each threshold
#' Note that whether averaging certain metrics is appropriate should be considered carefully
#' @param scoreList Predicted scores, List
#' @param labelList True class labels, List, paired with scoreList. Allowed class labels must be numeric 1 or 0
#' @param from Starting threshold, e.g., 1
#' @param to Ending threshold, e.g., 0
#' @param step Step size
#' @param Beta Parameter for F.measure. Beta controls the relative emphasis on sensitivity versus precision; when Beta > 1, more emphasis is placed on sensitivity; when Beta < 1, more emphasis is placed on precision; when Beta = 1,
#'             sensitivity and precision are weighted equally. Beta can be set as the ratio of the larger class size to the smaller class size to emphasize sensitivity
#' @returnType matrix 
#' @return Average value of each performance metric under different thresholds. Each row represents one threshold, and each column represents one metric.
#' 
#' @author XTeam
#' @export
MeanPerformanceMetrics <- function(scoreList, labelList, from=1, to=0, step=-0.01, Beta=1)
{
	if(length(scoreList) != length(labelList)){
		stop("scoreList与labelList without the same length")
	}
	
#	if(length(scoreList) == 1){
#		stop("length is 1，not using PerformanceMetricsList")
#	}
	
	#Generate the threshold sliding-window vector. Note that if from > to, step should be negative
	thresholds = seq(from, to, by=step) 
	
	
	num.element = length(scoreList);
	
	#For each threshold, compute all metrics for each element in the list, then average the corresponding metrics across elements
	results = sapply(thresholds, function(x){
				temp = sapply(1:num.element, function(i){
							PerformanceMetrics(scoreList[[i]], labelList[[i]], x, Beta)
						})
				temp = rowMeans(temp, na.rm=TRUE);
				return(temp)
			})
	
	return(t(results))
}