#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: Cross-platform quantile normalization to make datasets from different platforms share similar statistical distributions
#' 
#' 					 Perform quantile normalization of matrix columns based on a specified reference distribution	
#' (these function normalizes the columns of a matrix based upon a specified normalization distribution)
#' 
#' Input:
#' 				ref.matx     Reference matrix [rows = genes, columns = samples], used to obtain the reference distribution (i.e., quantile distribution vector)
#' 				query.matx   [rows = genes, columns = samples]; quantile normalization is performed on query.matx using the distribution of ref.matx as reference, so that its distribution becomes similar to ref.matx
#' Output:
#' 				qn.cat 		  Quantile-normalized result matrix [including ref.matx and the normalized query.matx]
#' 
#' Ref: https://github.com/greenelab/RNAseq_titration_results/blob/master/util/normalization_functions.R
#' 【Commun Biol, 2023, Ifs=5.9, Cross-platform normalization enables machine learning model training on microarray and RNA-seq data simultaneously】
#' Note: This quantile normalization does not handle outliers (No special allowances are made for outliers.)
#'#######
NormalizeQuantileCrossPlatform <- function(ref.matx, query.matx){
	
	#Data format transformation
	ref.df <- data.frame( ref.matx , check.names=FALSE)
	query.df <- data.frame( query.matx , check.names=FALSE)
	
	#' Point[1] Obtain quantile information from the reference dataset
	# Sort each column of the reference data separately, combine sorted columns into a matrix, and compute the row-wise mean of the sorted matrix as the reference quantile distribution
	qn.targ <- preprocessCore::normalize.quantiles.determine.target( data.matrix(ref.df) ) 
	#' Point[2] Perform quantile normalization on the query data
	# Using the above quantile information, normalize query.matx so that its distribution matches ref.matx
	qn.seq <- preprocessCore::normalize.quantiles.use.target(
			data.matrix(query.df),   #data to be normalized
			qn.targ,					 #reference distribution from ref.matx
			copy = F)
	#Combine the reference matrix and the quantile-normalized query matrix as the final matrix
	#qn.cat <- cbind(ref.df, qn.seq)
	#return(qn.cat)
	
	return(qn.seq)
}




#'##########  ---- Paragraph <附> ----  ###########  
#' Core objective: Single-platform quantile normalization
#' 
#' 		No reference distribution from other datasets is used
#' 
#' Input:
#' 				matx   Matrix to be quantile normalized within a single platform [rows = genes, columns = samples]
#' Output:
#' 				qn     Quantile-normalized matrix [rows = genes, columns = samples]
#' Ref: https://github.com/greenelab/RNAseq_titration_results/blob/master/util/normalization_functions.R
#' 【Commun Biol, 2023, Ifs=5.9, Cross-platform normalization enables machine learning model training on microarray and RNA-seq data simultaneously】
#' 
#' When using this function, the following paper should be cited: A Comparison of Normalization Methods for High Density Oligonucleotide Array Data Based on Bias and Variance
#' Note: This quantile normalization does not handle outliers (No special allowances are made for outliers.)
#'############ [23-11-17] ############
NormalizeQuantileSinglePlatform <- function(matx){
	
	# Sort each column of matx, combine sorted columns into a matrix, compute the row-wise mean of the sorted matrix as the reference quantile distribution,
	# and normalize matx based on this distribution
	matx = data.frame(matx, check.names=FALSE)
	qn <- preprocessCore::normalize.quantiles(data.matrix(matx), copy = F)
	
	return(qn)
}