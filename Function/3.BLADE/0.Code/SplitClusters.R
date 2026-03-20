#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: subdivide the specified cluster(s) in seurat.object
#' 
#' Point Input:
#' 				seurat.object seurat data object
#' 				clusters	  character or character vector; a specific cluster to be refined, or a set of clusters (cells in the set will be reclustered as a whole)
#' 				resolution	  clustering resolution
#'############ [24-09-09] ############
SubDivideCluster <- function(seurat.object, clusters, #cluster vector in CellCluster
		resolution=1){ #, default.assay="RNA"
	
	library(Seurat);library(SeuratObject);
	if(!all(colnames(seurat.object) == rownames(seurat.object@meta.data))){
		rownames(seurat.object@meta.data) = colnames(seurat.object)
	}
	
	#DefaultAssay(seurat.object) = default.assay
	t.seurat.object = subset(seurat.object, CellCluster %in% clusters)
	
	
	#Dimensionality reduction
	#Special note: iterative clustering may encounter a small number of cells, causing dimensionality reduction and clustering with high npcs to fail; the program will automatically choose an appropriate dimension and retry
	npcs_values <- c(30, 20, 10)   
	t.seurat.object <- tryCatch({
				# Try npcs = 30
				message("Running with npcs = 30...")
				t.seurat.object <- RunPCA(t.seurat.object, npcs = npcs_values[1])
				t.seurat.object <- RunUMAP(t.seurat.object, reduction = 'pca', dims = 1:npcs_values[1], n.neighbors = npcs_values[1])
				t.seurat.object <- FindNeighbors(t.seurat.object, dims = 1:npcs_values[1])
				t.seurat.object
			}, error = function(e) {
				message("Warning: RunPCA/RunUMAP/FindNeighbors failure with npcs = 30, retrying with npcs = 20.")
				tryCatch({
							# Try npcs = 20
							message("Retrying with npcs = 20...")
							t.seurat.object <- RunPCA(t.seurat.object, npcs = npcs_values[2])
							t.seurat.object <- RunUMAP(t.seurat.object, reduction = 'pca', dims = 1:npcs_values[2], n.neighbors = npcs_values[2])
							t.seurat.object <- FindNeighbors(t.seurat.object, dims = 1:npcs_values[2])
							t.seurat.object
						}, error = function(e) {
							message("Warning: RunPCA/RunUMAP/FindNeighbors failure with npcs = 20, retrying with npcs = 10.")
							message("Retrying with npcs = 10...")
							# Try npcs = 10
							t.seurat.object <- RunPCA(t.seurat.object, npcs = npcs_values[3])
							t.seurat.object <- RunUMAP(t.seurat.object, reduction = 'pca', dims = 1:npcs_values[3], n.neighbors = npcs_values[3])
							t.seurat.object <- FindNeighbors(t.seurat.object, dims = 1:npcs_values[3])
							t.seurat.object
						})
			})
	#https://github.com/satijalab/seurat/issues/1914
	#The default number of PCs to compute is 50, so you are trying to find more components than there are cells. You need to set npcs to less than the smallest dimension of the dataset
	
	
	
	#Clustering
	#Special note: iterative clustering may encounter a small number of cells, causing clustering with high npcs to fail; the program will automatically choose an appropriate setting and retry
	t.seurat.object <- tryCatch({
				FindClusters(t.seurat.object, resolution = resolution)
			}, error = function(e) {
				message("Warning: FindClusters failure, Retrying with group.singletons = FALSE.")
				message("Retrying with group.singletons = FALSE...")
				FindClusters(t.seurat.object, resolution = resolution, group.singletons = FALSE) #prevent errors caused by very small clusters 
				#https://github.com/satijalab/seurat/issues/2431
				#I have the same problem with a small data set. If I set group.singletons = FALSE then it runs. 
			})
	
	
	t.cluster.name = "seurat_clusters" #paste0(DefaultAssay(t.seurat.object), "_snn_res.", resolution)
	# Use the clustering result at the specified resolution as the final result
	t.seurat.object$CellCluster <- paste(clusters[1], t.seurat.object[[t.cluster.name]][,1], sep=".")
	
	seurat.object@meta.data[colnames(t.seurat.object),"CellCluster"] = t.seurat.object@meta.data[colnames(t.seurat.object),"CellCluster"]
	return(seurat.object)
}