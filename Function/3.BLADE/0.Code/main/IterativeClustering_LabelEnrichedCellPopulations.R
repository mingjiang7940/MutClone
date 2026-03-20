#'#############################  ---- Section <1> ----  ###############################
#'Core objective: 				Perform iterative clustering based on a Seurat object to identify clusters enriched for labels (mutant or wild-type)
#' 
#' 
#' 
#' Point Input:
#' [Required parameters]
#'         		seurat.object    				Seurat data object, a constructed atlas where data slot contains log2(cpm or tpm),
#'												Point Note: meta.data must include a True.lab column; labeled cells are assigned 0/1 (e.g., mutation status), and unlabeled cells are filled with NA
#'
#' Point Output:
#' 				list$
#' 						seurat.object    #Seurat object after iterative clustering of the atlas
#'								@meta.data
#'									.. ..$ atlas.X                    : num [1:84595] -2.585 -1.68 0.272 -1.6 -1.849 ... 
#'									.. ..$ atlas.Y                    : num [1:84595] 6.02 5.18 -1.23 5.17 2.99 ...
#'									.. ..$ atlas.CellCluster          : chr [1:84595] "3" "3" "6" "3" ...               #Atlas clustering result
#'									.. ..$ atlas.IterCluster          : chr [1:84595] "3.2" "3.10.3" "6.6" "3.20.0" ... #Iterative clustering result
#' 
#' 						clus.info	#Cluster-level information after iterative clustering
#' Special note:
#'Created on: September 4, 2024
#'##############################################################################
source( file.path(Device.path, "3.BLADE/0.Code/IterativeClustering.R") )
IterClusteringGetLabEnrichedClusters <- function(
		seurat.object, 
		#Developer parameters
		min.clus.count = 20,         #Do not further split clusters with fewer than 20 labeled cells
		major.clust.prop = 0.95, 	 #Proportion threshold for the dominant class
		max.iter = 10,		  		 #Stop splitting if iteration exceeds 10
		step.resolution=0.2,  		 #Resolution increment per iteration
		init.resolution=0.1   		 #Initial clustering resolution

){
	library(Seurat)
	seurat.object0 = seurat.object	   #Store the initial atlas						
	
	
	cat("1/2.迭代聚类..\n")  
	sc.dat = seurat.object
	
	#1.1) Initial clustering
	#Point Dimensionality reduction and clustering
	run_dimensionality_pipeline <- function(sc.dat, resolution=resolution, dims_use = 1:30) {
		cat("Running Seurat pipeline on seurat.all.cell...\n")
		sc.dat <- RunPCA(sc.dat, features = VariableFeatures(sc.dat))
		sc.dat <- FindNeighbors(sc.dat, dims = dims_use)
		sc.dat <- FindClusters(sc.dat, resolution=resolution)
		sc.dat <- RunUMAP(sc.dat, dims = dims_use, return.model = TRUE) #return.model must be TRUE for mapping
		return(sc.dat)
	}
	seurat.object = run_dimensionality_pipeline(sc.dat, resolution=init.resolution)
	sc.dat@meta.data[,"atlas.CellCluster"] = as.character( sc.dat@meta.data[,"seurat_clusters"] )
	
	
	#1.2) Sub-function getClustFun: define rules for iterative clustering; determine which clusters need further splitting (global variables used by getClustFun)
	t.min.clus.count <<- min.clus.count  				   #Temporary global variable
	t.major.clust.prop <<- major.clust.prop                #Temporary global variable
	
	#1.3) Set parameters and perform iterative clustering
	features = VariableFeatures(sc.dat)   				  #Features are stored in variable genes
	sc.dat@meta.data[,"CellCluster"] = as.character(sc.dat@meta.data[,"atlas.CellCluster"])      #CellCluster used for iterative clustering
	#options(error = recover)
	seurat.object2 = IterativeClustering(sc.dat, getClustFun=getClustFun, features=features, max.iter = max.iter, init.resolution=init.resolution+step.resolution, step.resolution=step.resolution )
	rm(t.min.clus.count, t.major.clust.prop, envir = .GlobalEnv) #Remove temporary global variables
	
	
	cat("2/2.更新seurat.object对象信息..\n")
	#Update seurat.object metadata
	seurat.object0@meta.data[,"atlas.IterCluster"] = seurat.object2@meta.data[colnames(seurat.object0),"CellCluster"] #Record iterative clustering results
	clus.info = getClusterInfo(seurat.object2)   #Get cluster-level information after iterative clustering
	
	res = list(seurat.object=seurat.object0,  #Final atlas after iterative clustering, can be used for mapping new datasets
			clus.info = clus.info )        #Cluster-level information
	
	return(res)
}




#'-Part Sub-function 1: determine which clusters need further splitting
getClustFun <- function(sc.dat){
	
	clus.info <- getClusterInfo(sc.dat)  #'-Part Sub-function 2: summarize cluster composition
	split.cluster.ls = c()
	for(clust in unique(clus.info$CellCluster) ){
		sub.clus.info = subset(clus.info, CellCluster==clust)
		total.cells = sum(sub.clus.info[,"n"])
		max.prop = max(sub.clus.info[,"prop"])
		
		if(total.cells>=t.min.clus.count & max.prop<t.major.clust.prop){ #Clusters with ≥20 cells and low label purity will be further split
			split.cluster.ls = c(split.cluster.ls, clust)
		}
	}
	
	return(split.cluster.ls)
}


#'-Part Sub-function 2: summarize cluster composition based on labels
#' 
#' Point Input:
#' 				seurat.object   Seurat object containing True.lab column
#' Point Output:
#> clus.info
#   CellCluster True.lab    n       prop
#1            0        0   76 0.29457364
#2            0        1  182 0.70542636
#3            1        0   92 0.06056616
#4            1        1 1427 0.93943384
#5           10        0  160 0.16048144
#6           10        1  837 0.83951856
#7           11        0  168 0.09778813
#8           11        1 1550 0.90221187
#9           12        0  178 0.13464448
#'#######
getClusterInfo <- function(seurat.object){
	
	#'-Point Identify cluster composition statistics
	#0) Summarize cluster information
	meta.data = seurat.object@meta.data
	# Group by CellCluster and compute statistics
	library(dplyr) 
	clus.info <- meta.data %>%
			filter(!is.na(True.lab)) %>%
			group_by(CellCluster, True.lab) %>%
			summarise(n = n(), .groups = "drop_last") %>%
			mutate(prop = n / sum(n)) %>%
			ungroup() %>%
			as.data.frame()
	
	return(clus.info)
}


