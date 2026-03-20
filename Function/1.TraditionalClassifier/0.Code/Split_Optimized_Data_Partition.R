#'#############################  ---- Section <1> ----  ###############################
#' Core objective: SPlit: An Optimal Method for Data Splitting
#' 
#' Note:
#' Created on: August 21, 2021
#'##############################################################################
#' Parameter dat: rows are genes and columns are samples; samples are split into training and test sets, and the indices of the test set are returned
#' Parameter label: numeric for regression problems and factor for classification problems
#' Point: the data type of label must be correct, because the algorithms for these two problem types are not exactly the same
#' splitRatio: proportion used for splitting (splitRatio = training samples / total samples)
#' #######
#Point output: positions of test-set samples; note that the returned positions correspond to the smaller subset (i.e., the test set)
SPlitPan <- function(dat,label,splitRatio=0.7){
	
	dat=t(dat)
	library(SPlit)
	#The input feature matrix must be normalized by column
	dat.scale = scale(dat)
	temp.dat = data.frame(dat.scale,y=label)
	SPlit_indices = SPlit(temp.dat, splitRatio=splitRatio)
	return(SPlit_indices)
}