#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: perform iterative clustering on a Seurat object
#' 
#' Point Input:
#' 				seurat.object    Seurat data object, where CellCluster stores the current clustering result
#' 				getClustFun      custom function that takes seurat.obj as input and returns a character vector of clusters that need further subdivision; this function is executed at each iteration to update the clusters to be split / if NULL, all current clusters are subdivided at each iteration
#' 				features		 vector of feature genes used for iterative clustering; if NULL, 3000 highly variable genes will be used
#' 
#' 				max.iter         maximum number of iterations
#' 				init.resolution  initial clustering resolution
#'              step.resolution  increment in clustering resolution at each iteration
#' 
#' Point Output:
#' 				seurat.obj after iterative clustering, with the iterative clustering result stored in CellCluster
#'############ [24-09-19] ############
source( file.path(Device.path, "3.BLADE/0.Code/SplitClusters.R") )
IterativeClustering <- function(seurat.object, getClustFun=NULL, features=NULL, max.iter = 10, init.resolution=1, step.resolution=0.2 ){
	
	library(Seurat)
	n.iter = 0
	resolution = init.resolution
	
	#1) Check whether the specified feature genes for iterative clustering are available
	if(!is.null(features)){
		VariableFeatures(seurat.object) = features
	}else if( length(VariableFeatures(seurat.object))==0 ){
		seurat.object = FindVariableFeatures(seurat.object, nfeatures = 3000)
	}
	
	#2) Check whether scaled data are available
	#tmp = GetAssay(sc.dat, assay = DefaultAssay(sc.dat))@scale.data
	tmp = GetAssay(seurat.object, assay = DefaultAssay(seurat.object))@scale.data
	if(is.null(tmp) | nrow(tmp) == 0){
		seurat.object = ScaleData(seurat.object)
	}
	
	#3) Check whether an initial cluster assignment exists
	if( !("CellCluster" %in% colnames(seurat.object@meta.data)) ){ #If there is no initial cluster assignment, assign all cells to one cluster
		seurat.object@meta.data[,"CellCluster"] = "0"  
	}
	
	#Iterative clustering
	repeat{
		n.iter = n.iter+1
		if( n.iter>max.iter )	break  #Stop iterative clustering if the maximum number of iterations is exceeded
		
		#Point Determine which clusters need further clustering: 1) if a rule function is provided for selecting clusters to enter iteration; 2) if no rule is provided, all clusters enter iteration
		if(!is.null(getClustFun)){
			split.cluster.ls = getClustFun(seurat.object)
		}else{
			split.cluster.ls = unique(seurat.object@meta.data[,"CellCluster"])
		}
		
		if( length(split.cluster.ls)==0 ) break
		
		#Point Further subdivide the selected clusters
		for(i in split.cluster.ls){
			
			seurat.object = SubDivideCluster(seurat.object, clusters=i, resolution=resolution) #, default.assay="RNA"
		}
		
		#Increase clustering resolution
		resolution = resolution+step.resolution 
		cat("完成",n.iter,"次迭代,当前精度",resolution,"\n")
	}
	return(seurat.object)	
}