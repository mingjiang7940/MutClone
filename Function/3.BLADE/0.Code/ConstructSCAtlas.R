# TODO: Add comment
#
#Point Special note: the constructed Atlas must be a standard Seurat object. Different methods can be used for construction, but the final result must be stored in a Seurat object
#
# Author: Administrator
###############################################################################
# Point Input: seurat.object is a Seurat data object, with data stored as log2(tpm+1) or log2(cpm+1)
# Point Output: output the mapped query.scDat
source( file.path(Device.path, "3.BLADE/0.Code/Atlas构建与新数据映射和拓展.gficf.R") )
source( file.path(Device.path, "3.BLADE/0.Code/seuratDataLog.R") )

ConstructSCAtlas <- function(
		seurat.object=NULL,   #integrated single-cell data
		seurat.obj.list=NULL, #unintegrated single-cell data (only one of the first two parameters can be specified)
		features=NULL,
		method="gficf",
		
		mode.path=NULL, #path for storing the special gficf model
		... )
{
	
	#(1) Point Determine which features are used to construct the Atlas
	if( is.null(features) ){
		hvgs <- tryCatch({
					VariableFeatures(seurat.object)
				}, error = function(e) {
					stop("No features or highly variable genes were provided")
				})
	}else{
		hvgs = features
	}
	
	
	#(2) Point Construct the Atlas using the specified method
	if(method=="gficf"){
		#Special note: when method="gficf", the data slot in seurat.object is automatically converted back from log scale to cpm or tpm
		seurat.object = scRemoveLog2(seurat.object)
		
		library(gficf)
		atlas = ConstructCellPoolAtlas.gficf(seurat.object, hvgs=hvgs,  method=method, ... )
		
		# store atlas information in the Seurat object
		meta.data = atlas$embedded
		colnames(meta.data)[colnames(meta.data) == "cluster"] <- "CellCluster"  #rename
		colnames(meta.data) = paste("atlas", colnames(meta.data), sep=".")
		
		seurat.object@meta.data[,colnames(meta.data)] = meta.data[colnames(seurat.object),]
		seurat.object@misc[["atlas.ClusteringResolution"]] = atlas$param$clustering.resolution
		
		#store the Atlas model
		saveGFICF(atlas, mode.path)
		
		seurat.object = scAddLog2(seurat.object)  #add log2 back after construction is completed
		
		#record model path
		seurat.object@misc[["atlas.ModePath"]] = mode.path
	}
	
	
	#(3) Point Seurat-based construction
	if(method=="Seurat RPCA"){
		
		#Point Integration
		if( is.null(seurat.object) & !is.null(seurat.obj.list) ){
			
			if(!is.null(features)){
				for(i in 1:length(seurat.obj.list)){ VariableFeatures(seurat.obj.list[[i]])=features }
			}
			source("/pub5/xiaoyun/BioU/chengmingjiang/function/单细胞数据处理函数/1.数据整合/Seurat v3/SeuratIntegrate.R")
			seurat.object = SeuratIntegrate(seurat.obj.list, mode="reference-based+RPCA", ... )
			
		}else if(!is.null(seurat.object) & is.null(seurat.obj.list)){
			
			cat("Input data have already been integrated; skipping the integration step\n")
			
		}else{
			stop("Error: seurat.object and seurat.obj.list cannot be both specified or both missing!")
		}
		
		#Point Dimensionality reduction and clustering
		run_dimensionality_pipeline <- function(sc.dat, dims_use = 1:30) {
			cat("Running Seurat pipeline on seurat.all.cell...\n")
			# [1] Identify highly variable genes
			sc.dat <- RunPCA(sc.dat, features = VariableFeatures(sc.dat))
			sc.dat <- FindNeighbors(sc.dat, dims = dims_use)
			sc.dat <- FindClusters(sc.dat)
			sc.dat <- RunUMAP(sc.dat, dims = dims_use, return.model = TRUE) #return.model must be included, otherwise mapping cannot be performed
			return(sc.dat)
		}
		
		seurat.object = run_dimensionality_pipeline(seurat.object)
		seurat.object@meta.data[,"atlas.CellCluster"] = as.character( seurat.object@meta.data[,"seurat_clusters"] )
		
		#record features
		if(is.null(features)) features = rownames(seurat.object)
	}
	
	seurat.object@misc[["atlas.method"]] = method #record method name
	seurat.object@misc[["atlas.features"]] = features
	
	return(seurat.object)
}





