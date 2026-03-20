#' Point [1]
#'Content: sub-function: add log2 transformation to Seurat object
scAddLog2 <- function( seurat.object ){
	
	tmp.dat = as.matrix( GetAssayData(seurat.object, assay="RNA", layer = "data") )
	tmp.dat = log2(tmp.dat+1)
	tmp.dat = as(tmp.dat, "dgCMatrix")
	seurat.object = SetAssayData(seurat.object, layer="data", new.data=tmp.dat)
	
	return( seurat.object )
}

#'Content: sub-function: remove log2 transformation from Seurat object
scRemoveLog2 <- function( seurat.object ){
	
	tmp.dat = as.matrix( GetAssayData(seurat.object, assay="RNA", layer = "data") )
	tmp.dat = 2^tmp.dat-1
	tmp.dat[tmp.dat<0] = 0
	tmp.dat = as(tmp.dat, "dgCMatrix")
	seurat.object = SetAssayData(seurat.object, layer="data", new.data=tmp.dat)
	
	return( seurat.object )
}