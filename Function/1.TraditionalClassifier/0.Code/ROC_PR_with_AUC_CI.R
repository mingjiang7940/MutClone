# TODO: Plot ROC and PR curves for classification performance evaluation, and return their areas under the curve.
# 
# Author: Administrator
###############################################################################

# Note: ROC sorts scores in ascending order and then calculates sensitivity and specificity threshold by threshold; therefore, higher-ranked scores are more likely to correspond to class 1
#' plot ROC and calculate AUC. If there are multiple ROC curves (not cross-validation curves), they can be plotted one by one and then added to the same figure
#' @param score : a score list. If multiple curves are to be displayed, they can be stored in a list, where each element corresponds to one curve
#' If isCrossValidation is FALSE, names should be provided for each element
#' @param label : a label list
#' @param type=c("ROC", "PR") ROC curve or PR curve
#' @param isCrossValidation=FALSE whether the results are from cross-validation, i.e., different elements in the list correspond to different x-fold results
#' @param all.fold.curvers=TRUE when isCrossValidation is TRUE, whether to plot all x-fold curves
#' 
#' @param legend.info vector of labels assigned to each curve; if NULL, labels are automatically set to AUC values
#' @param ... inherited plotting parameters
#' @returnType numeric
#' @return AUC and its confidence interval.
#'         If isCrossValidation is TRUE, the mean AUC across cross-validation folds is returned.
#'         If isCrossValidation is FALSE, the AUC value and confidence interval corresponding to each element in score are returned.
#' @author Administrator
#' @export
ROC_PR.normal <- function(score, label, type=c("ROC", "PR"), isCrossValidation=FALSE, 
		confidenceInterval=FALSE, #whether to return the confidence interval; if TRUE, return it as a formatted string
		is.plot=FALSE, #whether to draw the ROC curve
		all.fold.curvers=TRUE, 
		col=NULL,  	        	#Point [add] set the vector of line colors
		main=NULL,
		legend.info = NULL, 	#Point [add] set the annotation label for each curve
		lwd = 2,   	        	#Point [add] set line width
		show.baseline = FALSE,  #Point [add] whether to show the baseline of a random classifier in the PR curve
		avg="vertical",...)
{
	
	#score = lapply(score.label, function(x) x[,1]);
	#label = lapply(score.label, function(x) x[,2]);
	
	num.element <- length(score);
	names.element <- names(score);
	
	#library("colortools")
	if(is.null(col)){
		col = rainbow(num.element); #wheel("steelblue", num = 12);
	}
	
	## Plot ROC curve or PR curve
	if(type[1] == "ROC"){
		t.x = "fpr"; t.y = "tpr" 
	}else if(type[1] == "PR"){
		t.x = "rec"; t.y = "prec" 
	}
	
	library(ROCR);
	####Point revision log - April 24, 2025 [added calculation of the maximum y value to determine plotting range; otherwise the plot range would be determined only by the first curve]
	get_max_yval <- function(score, label, type = c("ROC", "PR")) {
		t.x <- ifelse(type[1] == "ROC", "fpr", "rec")
		t.y <- ifelse(type[1] == "ROC", "tpr", "prec")
		
		all.y <- c()
		for (i in 1:length(score)) {
			pred <- ROCR::prediction(score[[i]], label[[i]])
			perf <- ROCR::performance(pred, t.y, t.x)
			all.y <- c(all.y, unlist(perf@y.values))
		}
		return(max(all.y, na.rm = TRUE))
	}
	ylim.max <- get_max_yval(score, label, type = type)  # or "ROC"
	#############
	
	####Point revision log - April 27, 2025 added a function for correctly calculating AUC_PR
	##Previously ROCR was used: auc <- round(unlist(perf.auc@y.values), 3). Under class imbalance, this calculation may be unreliable
	#ROCR uses interpolation to avoid excessive fluctuations in the PR curve, but this may lead to over-interpolation when the data are very sparse (i.e., highly imbalanced positive and negative samples), resulting in inaccurate AUPRC estimates
	getAUC_PR <- function(score, label){
		library(PRROC)
		labels_numeric <- as.numeric(label)
		positive_scores <- score[labels_numeric == 1] #scores of positive samples
		negative_scores <- score[labels_numeric == 0] #scores of negative samples
		# Calculate the PR curve and AUC
		pr <- pr.curve(scores.class0 = positive_scores, #scores of positive samples
				scores.class1 = negative_scores, #scores of negative samples
				curve = TRUE)
		auc = round(pr$auc.integral, 3)
		return(auc)
	}
	########################
	
	##here used for cross-validatioin
	if(isCrossValidation & is.list(score) & is.list(label) & length(score) > 1 & (length(score) == length(label))){
		pred <- prediction(score,label);  ## generate prediction results
		perf <- performance(pred, t.y, t.x); ## calculate x and y coordinates of the curve
		
		## Calculate the area under the curve
		if(type[1] == "ROC"){
			perf.auc <- performance(pred,"auc"); 
			#The mean AUC from cross-validation is the average AUC across folds
			auc <- round(mean(unlist(perf.auc@y.values)),3); 
		}else if(type[1] == "PR"){
			perf.auc <- performance(pred,"aucpr"); 
			auc.list = lapply(1:length(score), function(i){
						getAUC_PR(score[[i]], label[[i]])
					})
			auc = round(mean(unlist(auc.list)),3);
		}
		
		
		
		#Draw curve 
		if(is.plot){
			tryCatch(
					{
						plot(perf,col=col,lwd=lwd, avg=avg, main=main, xaxs = "i", yaxs = "i", ylim = c(0, ylim.max), xlim = c(0, 1), ); ## draw the averaged ROC/PR curve from cross-validation
						if(all.fold.curvers){ plot(perf,col="grey82",lty=3, add=TRUE); } ## draw the curves from individual cross-validation folds
						
						####Revision log - April 23, 2025 [add baseline drawing]
						if (show.baseline) {
							if (type[1] == "ROC") {
								abline(a = 0, b = 1, lty = 2, col = "grey60", lwd = lwd )
							} else if (type[1] == "PR") {
								pos.ratio <- sum(unlist(label) == 1) / length(unlist(label))
								abline(h = pos.ratio, lty = 2, col = "grey60", lwd = lwd )
								# Add a small label on the left y-axis
								axis(side = 2, at = pos.ratio, labels = round(pos.ratio, 3),
										las = 1, tick = TRUE, col.axis = "grey40", cex.axis = 0.8)
							}
						}
						##############
						
						# Add legend 
						if(type[1] == "ROC"){
							info <- paste0("Avegage AUROC: ", auc, collapse="")
							legend("bottomright", legend=info, lwd=2, col=col)
						}else if(type[1] == "PR"){
							info <- paste0("Avegage AUCPR: ", auc, collapse="")
							legend("topright", legend=info, lwd=2, col=col)
						}
						
					},
					warning = function(w) { message('Waring') ; return(NA) },
					error = function(e) { message('Error') ; return(NA) },
					finally = { message('next...') }
			)
		}
		
		
		return(auc);
	}else{#not cross-validation
		
		all.auc <- vector(length=num.element)
		auc.ci.low <- vector(length=num.element)
		auc.ci.up <- vector(length=num.element)
		for(i in 1:num.element){
			pred <- prediction(score[[i]],label[[i]]); ## generate prediction results
			perf <- performance(pred, t.y, t.x);   ## calculate x and y coordinates of the curve
			
			## Calculate the area under the curve
			if(type[1] == "ROC"){
				perf.auc <- performance(pred,"auc"); 
				auc <- round(unlist(perf.auc@y.values), 3); ## extract the area under the curve
			}else if(type[1] == "PR"){
				#perf.auc <- performance(pred,"aucpr"); #this approach may be inaccurate when one class has few samples
				auc = getAUC_PR(score[[i]], label[[i]])
			}
			all.auc[i] = auc;
			
			## Calculate the confidence interval of AUC
			#library(pROC)
			#auc.ci.low[i] <- round(roc(label[[i]], score[[i]], ci = TRUE, auc = TRUE)$ci[1], 3)
			#auc.ci.up[i] <- round(roc(label[[i]], score[[i]], ci = TRUE, auc = TRUE)$ci[3], 3)
			#'direction parameter: You should set this explicitly to “>” or “<” whenever you are resampling or randomizing the data, otherwise the curves will be biased towards higher AUC values.
			#Set a fixed direction such that scores of class 0 are lower than those of class 1 [folds involve randomness, so this parameter must be explicitly set]
			auc.ci.low[i] <- round(pROC::roc(label[[i]], score[[i]], ci = TRUE, auc = TRUE, direction="<")$ci[1], 3)
			auc.ci.up[i] <- round(pROC::roc(label[[i]], score[[i]], ci = TRUE, auc = TRUE, direction="<")$ci[3], 3)
			
			par(mar = c(5, 5, 2, 2) + 0.1)
			if(is.plot){
				if(i == 1){
					plot(perf, col=col[i], lwd=lwd,  main=main, xaxs = "i", yaxs = "i", ylim = c(0, ylim.max), xlim = c(0, 1), ... );
				}else{
					plot(perf, col=col[i], lwd=lwd,  main=main, yaxs = "i", yaxs = "i", ylim = c(0, ylim.max), xlim = c(0, 1), add=TRUE, ... );
				}
				
				####Revision log - April 23, 2025 [draw baseline]
				if (show.baseline) {
					if (type[1] == "ROC") {
						abline(a = 0, b = 1, lty = 2, col = "grey60", lwd = lwd )
					} else if (type[1] == "PR") {
						pos.ratio <- sum(label[[i]] == 1) / length(label[[i]])
						abline(h = pos.ratio, lty = 2, col = "grey60", lwd = lwd )
						# Add a small label on the left y-axis
						axis(side = 2, at = pos.ratio, labels = round(pos.ratio, 3),
								las = 1, tick = TRUE, col.axis = "grey40", cex.axis = 0.8)
					}
				}
				########
			}
		}
		if(confidenceInterval)
		{
			info <- paste(names.element, " AUC: ", all.auc, "(", auc.ci.low, "-", auc.ci.up, ")", sep="")
		}else{
			info = all.auc
		}
		
		####Revision log - March 23, 2024
		#' Whether user-defined curve annotation labels are provided
		if(!is.null(legend.info)){
			for(i.line in 1:length(info)){
				info[i.line] = paste(legend.info[i.line], info[i.line])
			}
		}
		
		if(is.plot){legend("bottomright", legend=info, lwd=lwd, col=col)}
		return(info);
	}
	
}

#' If score and label are unavailable and only precomputed metrics under different thresholds are available, such as tpr and fpr, this function can be used to plot ROC and PR curves
#' @param x.values list, may have names, representing metrics under different thresholds
#' @param y.values list paired with x.values
#' If type is “ROC”, then x.values must be fpr and y.values must be tpr
#' @param type type=c("ROC", "PR") ROC curve or PR curve
#' @param col specify curve colors
#' @param main figure title
#' @returnType numeric vector 
#' @return Area under the curve corresponding to each element in x.values
#' 
#' @author XTeam
#' @export
ROC_PR.metrics <- function(x.values, y.values, type=c("ROC", "PR"), col=NULL, main=NULL)
{
	library(ROCR)
	num.element = length(x.values);
	names.element <- names(x.values);
	
	## Generate curve colors
	if(is.null(col)){
		col = rainbow(num.element);
	}
	
	## Determine whether to plot ROC or PR curves
	if(type[1] == "ROC"){
		x.name = "False positive rate";
		y.name = "True positive rate";
	} else if(type[1] == "PR"){
		x.name = "Recall";
		y.name = "Precision";
	} 
	
	
	## Calculate area under the curve 
	all.auc <- vector(length=num.element)
	for(i in 1:num.element){
		## Generate the same data structure as returned by performance
		perf = new("performance", x.name = x.name, y.name = y.name, 
				alpha.name = "Cutoff", x.values = list(x.values[[i]]), y.values = list(y.values[[i]]), 
				alpha.values = list(NULL)) 
		
		#Calculate AUC using the corresponding metrics, e.g., tpr and fpr
		auc <- round(simple_auc(y.values[[i]], x.values[[i]]), 3);
		all.auc[i] = auc;
		
		# Draw curve
		if(i == 1){
			plot(perf, col=col[i], lwd=2,  main=main, xaxs = "i", yaxs = "i");
		}else{
			plot(perf, col=col[i], lwd=2,  main=main, xaxs = "i", yaxs = "i", add=TRUE);
		}
		
	}
	info <- paste(names.element,": ",  all.auc, sep="")
	legend("bottomright", legend=info, lwd=2, col=col)
	
	return(all.auc);
}



#' Calculate AUC based on matched and ordered "True positive rate" and "False positive rate"
#' In fact, this can also be used to calculate the AUC of a PR curve, except that TPR is set as precision and FPR is set as recall
#' @param TPR vector of ordered "True positive rate" values
#' @param FPR vector of "False positive rate" values paired with TPR
#' @returnType numeric vector
#' @return Area under the curve 
#' 
#' @author XTeam
#' @export
simple_auc <- function(TPR, FPR){
	
	## Remove NA values
	t.index = which(is.nan(TPR));
	if(length(t.index) > 0){
		TPR = TPR[-t.index];
		FPR = FPR[-t.index];
	}
	
	## Remove NA values
	t.index = which(is.nan(FPR));
	if(length(t.index) > 0){
		TPR = TPR[-t.index];
		FPR = FPR[-t.index];
	}
	# inputs already sorted, best scores first 
	dFPR <- c(diff(FPR), 0)
	dTPR <- c(diff(TPR), 0)
	sum(TPR * dFPR) + sum(dTPR * dFPR)/2 ## calculate the area under the curve
}