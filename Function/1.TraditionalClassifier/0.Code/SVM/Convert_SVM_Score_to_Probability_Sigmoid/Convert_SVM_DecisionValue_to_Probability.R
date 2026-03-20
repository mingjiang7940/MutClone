#'#############################  ---- Section <1> ----  ###############################
#' Core objective: Convert SVM decision values to probabilities using Platt scaling
#' 
#' Part <1> Train Platt's sigmoid using cross-validation test scores and labels to map SVM outputs to probabilities
#' Part <2> Apply Platt's sigmoid to convert model prediction scores into probabilities
#' 
#' Package installation:
#' remotes::install_github("bioinf-jku/platt")
#' 
#'#References:
#' https://github.com/bioinf-jku/platt/blob/master/R/platt.R
#' https://github.com/bioinf-jku/platt/blob/master/R/predictProb.R
#' https://rdrr.io/github/bioinf-jku/platt/f/inst/doc/platt.pdf
#' Probabilistic Outputs for Support Vector Machines and Comparisons to Regularized Likelihood Methods, Platt, 1999
#'##############################################################################

library(platt)
#' ####--- Part <1> ----###
#' Core objective: Train Platt's sigmoid (plattScaling) to convert SVM outputs into probabilities
#' Input:
#' 				predictions Prediction score vector from cross-validation test sets
#' 				labels      Corresponding binary labels (0/1) for the cross-validation test samples
#' Output: 
#' 				plattScalingResult   Trained Platt's sigmoid transformer for converting SVM outputs into probabilities
#'############ [23-12-03] ############
plattScalingResult = plattScaling(predictions, labels)


#'#### --- Part <2> ----###
#' Core objective: Convert new prediction scores into probabilities using Platt's sigmoid
#' Input:
#' 				pS            	  Trained Platt's sigmoid
#' 				new.predictions   Prediction scores to be transformed into probabilities
#' Output:
#' 				Prob          Converted probabilities
#'############ [23-12-03] ############
Prob = predictProb(plattScalingResult, new.predictions)
