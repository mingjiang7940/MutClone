#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: obtain the optimal threshold through ROC analysis
#' 
#' Point Input:
#' 				label    sample labels
#' 				score  predicted scores for samples
#' 				pos    specify the label of the positive class
#' 				best.method  criterion for selecting the optimal threshold; default is the maximum Youden index, i.e. youden index = sensitivities + specificities - 1 is maximized
#' 							 optional: "closest.topleft" (threshold closest to the top-left corner of the ROC curve), "Fscore"
#' Point Output：
#' 
#' 
#To account for sample imbalance, an r term was added here, and the optimization target becomes max(sensitivities + r\times specificities)
#'############ [22-11-19] ############
RocGetThreshold <- function(score, label, pos = NULL, best.method = "youden") {
	
	label <- factor(label, levels = c(setdiff(unique(label), pos), pos))
	rocobj <- pROC::roc(label, score, direction = "<")
	
	if (best.method == "Fscore") {
		
		thresholds.tab <- pROC::coords(rocobj, x = "all", ret = c("threshold", "sensitivity","precision"))
		thresholds.tab$f1_scores <- 2 * thresholds.tab[,"sensitivity"] * thresholds.tab[,"precision"] /
				(thresholds.tab[,"sensitivity"] + thresholds.tab[,"precision"])
		thresholds.tab = na.omit(thresholds.tab)
		
		max_f1_index <- which.max(thresholds.tab$f1_scores)[1]
		best_threshold <- thresholds.tab[max_f1_index,"threshold"]
		
	} else {
		roc_result <- pROC::coords(rocobj, "best", best.method = best.method)
		best_threshold <- if (prod(dim(roc_result)) > 3) as.numeric(roc_result[1, 1]) else as.numeric(roc_result[1])
	}
	
	return(best_threshold)
}
