#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: Construct a parameterized SVM model (e.g., compatible with svmTL format)
#' 
#' Part[operate=="construct"] 
#'		 Input:
#		 				data  Feature matrix (rows = genes, columns = samples)
#		 				label Binary labels (0/1) corresponding to samples
#
#						(Optional parameters)
#						...  
#						Other parameters from alg1_baselineClass can also be provided:
#						For example:
#							 model.names  #If MergeSvmModels will be used later, it is recommended to assign a name to the model
#							 cost         SVM cost parameter; smaller values improve generalization (default = 0.5). It represents the penalty for constraint violations, corresponding to the regularization constant "C" in the Lagrangian formulation
#							 cross        If a positive integer is specified, k-fold cross-validation is performed on training data to evaluate model quality
#
#'		Output: Model format compatible with svmTL
#					List of 2
#							 $ baselineSVM:  		#Model matrix used by svmTL
#							 $ model.list :List of 1  
#							  ..$ LUSC:List of 3    #Model name and structure
#							  .. ..$ W           : num [1:5287] 2.32e-04 6.15e-05 -2.03e-05 -1.68e-05 1.12e-04 ...
#							  .. ..$ b           : num -0.401
#							  .. ..$ use.features: chr [1:5287] "FAM50B" "AMER3" "FAM163A" "DEFA5" ...
#							  .. ..$ probA       : num -4.52    #[internally trained by SVM] sigmoid parameter for converting decision values to probabilities
#							  .. ..$ probB       : num 0.192
#
#' Part[operate=="predict"] 
#' 				Input:
#			 				model Model used for prediction
#			 				data  Dataset for prediction (rows = genes, columns = samples)
#'				Output: 
#							tmp   Predicted probabilities from the SVM model
#'############ [24-04-05] ############
source(file.path(Device.path, "2.TransferLearning", "svmTL","0.Code","BasicAlg", "Alg1_BaselineSVMClassifier.R"))

ClassifierSVM.svmTL <- function(data=NULL, label=NULL, model=NULL, operate=c("construct", "predict")[1], cross=3, cost=0.5, model.names=NULL, ...)
{
	#Data format transformation
	if(operate == "construct"){
		#' Convert data into RTLbase input format
		label[label==0] = -1
		Ds.X = list(t(data))
		Ds.Y = list(label)
		if(!is.null(model.names)){names(Ds.X) = names(Ds.Y) = model.names}
		#' Adjust sample order to ensure consistency of SVM decision value signs across sources (positive = mutant, negative = wild-type); samples are ordered as label 1 -> -1
		Ds.Y = lapply(Ds.Y,function(x) sort(x, decreasing=TRUE) )
		Ds.X = mapply(function(x,y) x[names(y),], Ds.X, Ds.Y ,SIMPLIFY=FALSE,USE.NAMES = TRUE)
		##########Point 1. Construct baseline SVM models###########
		#' alg1_baselineClass: Fit linear SVM hyperplanes in the source domains
		#' The results are (w_1,b_1),(w_2,b_2), ... ,(w_m,b_m)
		#' Stored in alg1_res as input for Point 2
		######################################
		suppressMessages(library('e1071'))
		suppressMessages(library('caret'))
		alg1_res <- alg1_baselineClass(
				TrainXls = Ds.X, #each matrix: rows = samples, columns = features
				TrainYls = Ds.Y,	#labels are 1 and -1
				TestXls = NULL,
				TestYls = NULL,
				X_cols2Keep = rownames(data),  #use all genes in the matrix for modeling
				svmCost = cost,
				K_forCrossV = cross,
				prnt2scr = FALSE,
				transX=F, sampleRed=F, doParalellSVM=F, datatyp="Exprs", ...) #
		return(alg1_res)
	}
	
	
	#'-Part [2]- Prediction
	if(operate == "predict"){
		RTLmodel.predict <- function(model,Data){
			#' Input:
			# 		 model = list(W, 				#linear SVM parameter w
			# 		 			  b,				#linear SVM parameter b
			# 		 			  use.features)	    #features used in the model
			# 		 Data 						matrix with rows = genes and columns = samples
			#' Output:
			# 		 y.value     Predicted decision values
			# 					 y.value = t( as.matrix(Data[model$use.features,]) ) %*% model$W + model$b
			###########################2021年11月14日#######################
			y.value = t( as.matrix(Data[model$use.features,]) ) %*% model$W + model$b
			return(y.value[,1])
		}
		#pred.value.svm = RTLmodel.predict(alg1_res$model.list[[1]], data)
		t.model = model$model.list[[1]]
		y = RTLmodel.predict(t.model, data) #output decision values
		temp = 1/(1+exp(t.model$probA*y+t.model$probB))      #sigmoid transforms decision values into probabilities
		return(temp)
	}
}