#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: Extract detailed performance metrics from ClassifyModel
#' Point input:
#' 				ClassifyModel Classification model object or transfer learning model
#' Point output:
#' 				performance   A vector of model performance metrics
#'################################################
ExtractPerformanceMetrics <- function(ClassifyModel){
	
	# 1) Extract model performance
	AUROC.rep = unlist( lapply(ClassifyModel$CV.PP, function(CV.PP) CV.PP$AUROC  ) ) #retrieve performance from each repetition
	AUPR.rep = unlist( lapply(ClassifyModel$CV.PP, function(CV.PP) CV.PP$AUPR  ) )
	
	performance = c(AUROC=mean(AUROC.rep), AUPR=mean(AUPR.rep))
	
	return(performance)
}
