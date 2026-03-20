# TODO: Add comment
# Author: Administrator
#
# ##Original version without considering clonality
# Preprocessing 3.0_predict mutation genes based on expression
#
# ##All versions considering clonality
# Point Standard log2(tpm+1): Preprocessing 3.1_predict clonal mutation genes based on expression
# Point TCGA pan-cancer data are already batch-corrected: Preprocessing 3.1.2_predict clonal mutation genes based on expression
# Point Our own quantile normalization for pan-cancer data: Preprocessing 3.2_predict clonal mutation genes based on expression_quantile normalization version
###############################################################################
rm(list = ls() )
Device.path = "D:/Project/0.MutClone/Function"
ID.xteam = "CRC_TCGA"
out.dir = file.path(dirname(Device.path),"0.Data/Results", ID.xteam, "2.svmTL")
if(!dir.exists(out.dir)){ dir.create(out.dir, recursive = TRUE) }

way="mut-exprs"
background=FALSE
method.select.feature="SVM-RFE"
method.select.source="DTE"
nCore = 1
options(future.globals.maxSize = 5 * 1024^3)  # 5GB

pancancer.driver = readRDS( file.path(dirname(Device.path), "0.Data", "CancerDriverGene.rds") )
pancancer.driver = unique(unlist(pancancer.driver))
file.mut.exprs.data = file.path(dirname(Device.path),"0.Data",'TCGA_MutClone', ID.xteam, "mut.pancancer.exprs.data.indel.RDS")
mut.exprs.data = readRDS(file.mut.exprs.data)
mutgenes = names(mut.exprs.data$all.label)
mutgenes = mutgenes[mutgenes %in% pancancer.driver]

mutgenes = "APC"

if(length(mutgenes) == 0){
	cat("No mutation genes in this dataset satisfy the specified criteria.")
}else{
	
	source( file.path(Device.path,"2.TransferLearning/svmTL/0.Code/0.BuildTLDataFramework.R") )
	source( file.path(Device.path,"2.TransferLearning/svmTL/0.Code/main/SVMTL_FS_source_selection_cross_validation.R") )
	source( file.path(Device.path,"1.TraditionalClassifier/0.Code/main/2.OptimizeThresholdCV.R") )
	source( file.path(Device.path,"1.TraditionalClassifier/0.Code/main/3.ConstructFinalModel.R") )
	out.file.svmTL = file.path(out.dir, "svmTL.rds")
	
	if( file.exists(out.file.svmTL) ){
		cat("svmTL already exists")
	}else{

		ID.xteams = c("ACC_TCGA", "BLCA_TCGA", "BRCA_TCGA", "CESC_TCGA", "CHOL_TCGA", "CRC_TCGA",
					"DLBC_TCGA", "ESCA_TCGA", "GBM_TCGA", "HNSC_TCGA", "KICH_TCGA", "KIRC_TCGA",
					"KIRP_TCGA", "LGG_TCGA", "LIHC_TCGA", "LUAD_TCGA", "LUSC_TCGA", "MESO_TCGA",
					"OV_TCGA", "PAAD_TCGA", "PCPG_TCGA", "PRAD_TCGA", "SARC_TCGA", "SKCM_TCGA",
					"STAD_TCGA", "TGCT_TCGA", "THCA_TCGA", "THYM_TCGA", "UCEC_TCGA", "UCS_TCGA",
					"UVM_TCGA")
		ID.xteams = setdiff(ID.xteams, ID.xteam)#将自己去除
		

		library(parallel)
		result.cm = lapply(mutgenes, function(t.name){
					
			
					#(0) 识别可用的源
					valid_ID_xteams = c()  # 保留包含 'mutgene' 的 ID_xteams
					for(tmp.id in ID.xteams){
						s.file.mut.exprs.data = file.path(dirname(Device.path),"0.Data",'TCGA_MutClone', tmp.id, "mut.pancancer.exprs.data.indel.RDS")
						s.mut.exprs.data = readRDS(s.file.mut.exprs.data)
						if(t.name %in% names(s.mut.exprs.data$all.label)){
							valid_ID_xteams = c(valid_ID_xteams, tmp.id)
						}
					}
					
					if( length(valid_ID_xteams)==0 ) return(NULL)
					
					#(1.1) 建立源模型
					cat("源模型建立中...\n")
					s.model.list = list()
					for(tmp.id in valid_ID_xteams){

						s.file.mut.exprs.data = file.path(dirname(Device.path),"0.Data",'TCGA_MutClone', tmp.id, "mut.pancancer.exprs.data.indel.RDS")
						source( file.path(Device.path,"1.TraditionalClassifier/0.Code/main/1.FeatureSelectionCV.R") )
						cm = ClassifierCV.FeatureSelection(data=NULL,
								label=NULL,
								data.name=t.name,
								data.path.rds=s.file.mut.exprs.data,
								classifier="linearSVM.svmTL",
								topn=seq(10, 500, 10),
								times=1)
						
						source( file.path(Device.path,"/1.TraditionalClassifier/0.Code/main/2.OptimizeThresholdCV.R") )
						cm = ClassifierCV.OptimizeThreshold(cm)

						source( file.path(Device.path,"1.TraditionalClassifier/0.Code/main/3.ConstructFinalModel.R") )
						cm = ConstructFinalModel(cm)
						
						s.model.list[[tmp.id]] = cm
					}
						
					
			
					cat("迁移模型建立中...\n")
					#(2) Assemble into a "standard transfer learning model"
					#Point Assemble the complete transfer learning data object
					tl.cm = CreateTLClassifyModel(
							data=NULL, #rows are features, columns are samples
							label=NULL, #label vector: 1 indicates positive class, 0 indicates negative class
							data.name=t.name, data.path.rds=file.mut.exprs.data,
							classifier="SVMTL", 
							method.select.feature=method.select.feature,
							s.model.list=s.model.list, 
							method.select.source="DTE")
					
					
					cat("超参数选择...\n")
					#'##########  ---- Paragraph <1> ----  ###########  
					#' Core objective: hyperparameter selection based on cross-validation
					#' 
					#'############ [24-04-30] ############
					cm = ClassifierCV.FeatureSelection.SVMTL(cm=tl.cm, 
							way=way, background=background,
							classifier=c("SVMTL")[1], 
							cross.repeat=5, times=1,
							method.select.feature=method.select.feature, 
							topn.feature=c(30, 50, 100, 200, 300, 400, 500),
							method.select.source=method.select.source, 
							topn.source=1:length(valid_ID_xteams),
							nCore=nCore
						)
					
		
					cm = ClassifierCV.OptimizeThreshold(cm=cm, method.threshold="Fscore")
					
					
					cm = ConstructFinalModel.SVMTL(cm=cm, way="mut-exprs", background=FALSE)
					
					t.file = file.path(out.dir, "单个基因", paste0(t.name, ".rds"))
					if(!dir.exists( dirname(t.file) )){ dir.create(dirname(t.file), recursive = TRUE) }
					
					saveRDS(cm, file=t.file)
					return(cm)
				})
		names(result.cm) = mutgenes
		
		result.cm = result.cm[!sapply(result.cm, is.null)]
		saveRDS(result.cm, file=out.file.svmTL)
	}
}



