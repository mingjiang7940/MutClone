#'#############################  ---- Section <1> ----  ###############################
#'Core objective:					KNN smoothing of the expression matrix
#' 
#' Point Input:
#' 				dat normalized expression matrix (e.g., cpm, tpm, etc.) [rows = genes, columns = cells/samples]
#' 				k	number of nearest neighbors used for smoothing, default 50
#' 				mc.cores number of parallel cores used for computation, default 10
#' 
#' Point Output: 
#' 				knn.mat smoothed matrix
#' Special note:
#'Created on: November 3, 2023
#' Ref: Cancer classification of single-cell gene expression data by neural network, Bioinformatics, 2020
#'##############################################################################
#Calculate the distance matrix, identify the k nearest neighboring cells, and average them to smooth the cell expression profiles
KnnSmoothing <- function(dat, k=50, mc.cores=10){
	
	if( !is.matrix(dat) ) dat=as.matrix(dat)
	#(1) Calculate the distance matrix
	library(distances)
	Dist = distances( t(dat) )
	
	#(2) KNN smoothing
	smoothing.dat = list() #Each element stores the smoothed expression vector for one cell
	library(parallel)
	smoothing.dat = mclapply(as.list(1:ncol(dat)), function(i){
				
				iname = colnames(dat)[i]
				# 1) Find the k nearest neighbors of cell iname
				iManhattan.dist = Dist[which(colnames(dat)==iname),]
				knn.names = colnames(dat)[order(iManhattan.dist)[1:k]]
				
				# 2) Smooth the expression profile of cell iname
				tmp.dat = rowMeans(dat[,knn.names])
				cat(i,"/",ncol(dat)," run...\n")
				
				return(tmp.dat)		
			},
			mc.cores=mc.cores )
	
	#Convert smoothing results to a matrix
	knn.mat = matrix(unlist(smoothing.dat), length(smoothing.dat[[1]]), length(smoothing.dat))
	colnames(knn.mat) = colnames(dat)
	rownames(knn.mat) = rownames(dat)
	
	return(knn.mat)
}

