# TODO: Add comment
# 
# Author: Administrator
###############################################################################

#source("/pub5/xiaoyun/BioX/Bioc/BioSingleCell/伪bulk构建/2.对normalized后count值取平均值/F.KnnSmoothing.R")
#'-Part 子函数1：seurat在患者内构建伪bulk
#' 
#' Point 输入：
#' 				seurat.object     seurat对象, data的单位为log2(cpm+1) 或者 log2(tpm+1)
#' 				k          选最近邻多少个细胞构建伪bulk
#' Point 输出：
#' 				dat.matrix.knn  伪bulk矩阵, 注意输出矩阵细胞可能比seurat.object少, 因为当患者细胞数少于k时, 不构建其细胞伪bulk
#'#######
source( file.path(Device.path, "3.BLADE/0.Code/F.KnnSmoothing.R") )
psudoBulkInPat <- function(seurat.object, k=50, mc.cores=10){
	
	# 细胞在患者内中用最近邻50个细胞构建伪bulk, 癌细胞数目少于50的患者被排除
	#k = 50
	meta.dat = seurat.object@meta.data
	dat.matrix = 2^GetAssayData(seurat.object, assay="RNA", layer = "data")-1
	dat.matrix[dat.matrix<0] = 0
	patient_PseudoBulk.list = list()
	for(patient in unique( meta.dat$PatientID )){ 
		
		cells = rownames(meta.dat)[ meta.dat$PatientID==patient ]
		tmp.dat = dat.matrix[,cells, drop=F]
		
		if(ncol(tmp.dat)>=k){
			tmp.cc.knn = KnnSmoothing(tmp.dat, k=k, mc.cores=mc.cores)
			patient_PseudoBulk.list[[patient]] = tmp.cc.knn
		}	
	}
	dat.matrix.knn = do.call(cbind, patient_PseudoBulk.list)
	dat.matrix.knn = log2(dat.matrix.knn+1)
	
	return(dat.matrix.knn)
}

