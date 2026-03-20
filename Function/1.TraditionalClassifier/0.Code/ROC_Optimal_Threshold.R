#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: Obtain the optimal threshold by ROC analysis
#' 
#' Input:
#' 				label    Sample labels
#' 				score    Predicted scores for samples
#' 				pos      Specify the label of the positive class
#' 				best.method  Criterion for obtaining the optimal threshold; default is maximizing the Youden index, i.e., sensitivities + specificities - 1
#' 							 Optional: "closest.topleft" (threshold closest to the upper-left corner of the ROC curve), "Fscore"
#' Output:
#' 
#' 
#To account for class imbalance, an r term can be introduced here, changing the optimization objective to max(sensitivities + r\times specificities)
#'############ [22-11-19] ############
RocGetThreshold <- function(score, label, pos = NULL, best.method = "youden") {
	
	label <- factor(label, levels = c(setdiff(unique(label), pos), pos))
	rocobj <- pROC::roc(label, score, direction = "<")
	
	if (best.method == "Fscore") {
		
		thresholds.tab <- pROC::coords(rocobj, x = "all", ret = c("threshold", "sensitivity","precision"))
		thresholds.tab$f1_scores <- 2 * thresholds[,"sensitivity"] * thresholds[,"precision"] /
				(thresholds[,"sensitivity"] + thresholds[,"precision"])
		thresholds.tab = na.omit(thresholds.tab)
		
		max_f1_index <- which.max(thresholds.tab$f1_scores)[1]
		best_threshold <- thresholds[max_f1_index,"threshold"]
		
	} else {
		roc_result <- pROC::coords(rocobj, "best", best.method = best.method)
		best_threshold <- if (prod(dim(roc_result)) > 3) as.numeric(roc_result[1, 1]) else as.numeric(roc_result[1])
	}
	
	return(best_threshold)
}