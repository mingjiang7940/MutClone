# TODO: Add comment
# 
###############################################################################
# Seurat v3 integration mainly consists of two steps:
# 1. FindIntegrationAnchors: identifies pairs of cells across datasets that share similar biological states,
#    referred to as anchors. These anchors are used for subsequent integration.
# 2. IntegrateData: integrates datasets using the previously identified anchor set.
###
#' @param seurat.obj.list: A list of Seurat objects.
#' @param ndim: The number of dimensions used when projecting datasets into a low-dimensional space via CCA.
#'              Default is 30.
#' @param mode: Integration mode. Options include "standard", "reference-based", "RPCA",
#'              and "reference-based+RPCA". Default is "standard".
#'              1. standard: The standard workflow. Anchors are identified between every pair of datasets.
#'                 This can be slow for large datasets. For large-scale data, options 2–3 are recommended,
#'                 as results are comparable to the standard approach.
#'              2. reference-based: A subset of datasets is used to construct a reference, and the remaining
#'                 datasets are treated as queries. Queries are only compared to the reference, greatly
#'                 reducing the number of comparisons and improving speed.
#'              3. RPCA (Reciprocal PCA): Uses Reciprocal PCA instead of CCA for dimensionality reduction,
#'                 improving speed and efficiency for large datasets.
#'              4. reference-based+RPCA: Combines reference-based integration with RPCA.
#' @param reference: If using a mode that includes "reference-based", specify the indices of reference datasets.
#'                  Default is the first two datasets, i.e., c(1,2).
#'                  Note: It is recommended to choose datasets with a large number of cells as references;
#'                  otherwise, errors may occur during reference construction due to insufficient cell numbers.
#'
#' @return A Seurat object. The integrated results are stored in the "integrated" assay,
#'         and DefaultAssay(t.seurat.integrated) returns "integrated".
#'         If needed, the default assay can be switched back to normalized data ("data")
#'         or raw counts ("counts") using DefaultAssay.
#'
#' @example
#' load("/pub6/temp/t.scRNAseq.1.31.RData")
#' seurat.obj.list = t.scRNAseq[1:3]
#' t.seurat.integrated = SeuratIntegrate(seurat.obj.list, ndim = 30, mode = "standard")
#' https://satijalab.org/seurat/archive/v3.2/integration
################## November 1, 2020 #############

SeuratIntegrate=function(seurat.obj.list,ndims=30,mode="standard",reference=c(1,2),anchor.features.n=3000,...)
{
  library(Seurat)
  
  
  #############################standard模式
  if(mode=="standard" | is.null(mode)){
    t.anchors = FindIntegrationAnchors(object.list = seurat.obj.list, dims = 1:ndims)
    seurat.obj.integrated = IntegrateData(anchorset = t.anchors, dims = 1:ndims)
    
  }
  
  
  
  ###########################加速模式1--rpca
  
  if(mode=="RPCA"){
    features = SelectIntegrationFeatures(object.list = seurat.obj.list)
    seurat.obj.list = lapply(X = seurat.obj.list, function(x) {
      x = ScaleData(x, features = features, verbose = FALSE)
      x = RunPCA(x, features = features, verbose = FALSE)
    })
    
    anchors = FindIntegrationAnchors(object.list = seurat.obj.list,reduction = "rpca", 
                                      dims = 1:ndims)
    seurat.obj.integrated = IntegrateData(anchorset = anchors, dims = 1:ndims,...)
    
  }
  
  
  #############################加速模式二---reference-based
  if(mode=="reference-based"){
    for (i in 1:length(seurat.obj.list)) {
      seurat.obj.list[[i]] = SCTransform(seurat.obj.list[[i]], verbose = FALSE)
    }
    seurat.features = SelectIntegrationFeatures(object.list = seurat.obj.list, nfeatures = anchor.features.n)
    seurat.obj.list = PrepSCTIntegration(object.list = seurat.obj.list, anchor.features = seurat.features)
    reference_dataset = reference
    anchors = FindIntegrationAnchors(object.list = seurat.obj.list, normalization.method = "SCT", 
                                      anchor.features = seurat.features, reference = reference_dataset)
    seurat.obj.integrated = IntegrateData(anchorset = anchors, normalization.method = "SCT",...)
  }
  
  
  
  #############################加速模式三--rpca+reference-based
  if(mode=="reference-based+RPCA"){
    seurat.obj.list = lapply(seurat.obj.list,SCTransform);
    seurat.features = SelectIntegrationFeatures(object.list = seurat.obj.list, nfeatures = anchor.features.n)
    seurat.obj.list = PrepSCTIntegration(object.list = seurat.obj.list, anchor.features = seurat.features)
    seurat.obj.list = lapply(X = seurat.obj.list,function(x) {
      x = ScaleData(x, features = seurat.features, verbose = FALSE)
	  x = RunPCA(x, features = seurat.features, verbose = FALSE)
      #x = RunPCA(x, features = seurat.features, verbose = FALSE,approx=FALSE) #You're computing too large a percentage of total singular values, use a standard svd instead.
    })
    reference_dataset = reference
    anchors = FindIntegrationAnchors(object.list = seurat.obj.list, normalization.method = "SCT", reduction = "rpca", 
                                      anchor.features = seurat.features, reference = reference_dataset)
    seurat.obj.integrated = IntegrateData(anchorset = anchors, normalization.method = "SCT",...)
  }
  return(seurat.obj.integrated)
}


DirectCombineSeurat <- function(seurat.obj.list)
{
	seurat.obj = merge(seurat.obj.list[[1]], seurat.obj.list[2:length(seurat.obj.list)])
	return(seurat.obj)
}
