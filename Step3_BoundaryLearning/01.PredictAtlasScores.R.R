#'##########  ---- Section <1> ----  ###########  
#' Core objective: examine the distribution of predicted scores at the cell level when directly using sample-level labels
#'                      based on model prediction score distributions
#' 
#' Paragraph <1>    predict cells in the iniAtlas
#'              
#' 
#'############ [25-11-11] ############
rm(list=ls())
Device.path = "D:/Project/0.MutClone/Function"
out.dir = "/WorkSpace/chengmingjiang/1.Project_sc2mutTL/3.Section/00.完成故事梳理后_整理版/S3_模型偶联Atlas/0.Data/新版本/01.能直接用样本标签代替细胞标签吗"
if (!dir.exists(out.dir))	dir.create(out.dir, recursive = TRUE)
cancer = "CRC"
mut.ls = c("APC","KRAS","TP53","BRAF","PIK3CA","SMAD4","FBXW7")




#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective:    predict cells in the iniAtlas
#' 
#' -Part [0]- load data and models
#' 
#' -Part [1]- predict the initial Atlas
#'############ [26-01-16] ############
#'-Part [0]-
#' Content: load data and models
#'#######
#iniAtlas = readRDS( file.path(dirname(out.dir),"initAtlas.RDS") )
#seurat.object_pseudo_bulk = readRDS( file.path(dirname(out.dir),"initAtlas.pseudo_bulk.RDS") )
iniAtlas = readRDS( file.path(dirname(out.dir),"initAtlas_newCellNames.RDS") )
seurat.object_pseudo_bulk = readRDS( file.path(dirname(out.dir),"initAtlas.pseudo_bulk_newCellNames.RDS") )
model.list = list()
for(i in 1:length(mut.ls) ){ 
	
	mutation = mut.ls[i]
	mod.id = paste("CRC",mutation,sep="_")
	mod.typs = c("svm", "svmTL", "dnnTL")
	for(mod.typ in mod.typs){
		
		#load models
		if(mod.typ=="svmTL"){
			t.file = file.path("/IData/DataCenter/TCGA",paste0(cancer,"_TCGA"),"Results/BioGenomics/01.突变/93.预测突变/11.基于表达预测突变/120.基因突变谱与表达谱[带分类标签]/2.SVMTL/svmTL.rds")
			svmTL.list <- tryCatch({ readRDS(t.file) }, error = function(e){NULL} )
			cm.bulk = svmTL.list[[mutation]]
		}else if(mod.typ=="svm"){
			t.file = file.path("/IData/DataCenter/TCGA",paste0(cancer,"_TCGA"),"Results/BioGenomics/01.突变/93.预测突变/11.基于表达预测突变/120.基因突变谱与表达谱[带分类标签]/1.SVM/SVM.indel.rds")
			cm.bulk.list = readRDS(t.file)
			cm.bulk = cm.bulk.list[[mutation]]
		}else if(mod.typ=="dnnTL"){
			library(reticulate)
			use_condaenv("/Software/PythonENV/DeepLearning.py3.9", required = TRUE) 
			py_config()
			DevicePath = 'D:/Project/0.MutClone/Function'
			source_python( file.path(DevicePath,"2.TransferLearning/dnnTL/0.Code/main", "ClassifyModel.py") ) #	
			#py$directory = file.path("/WorkSpace/chengmingjiang/TmpData/TL_model",paste0(cancer,"_TCGA"),"2.迁移学习dnnTL/单个基因3/")
			py$directory = file.path("/IData/DataCenter/TCGA",paste0(cancer,"_TCGA"),"Results/BioGenomics/01.突变/93.预测突变/11.基于表达预测突变/120.基因突变谱与表达谱[带分类标签]/3.dnnTL/单个基因/")
			py$cancer = cancer
			py$target_genes = mutation
			
			# Use reticulate's py_run_string to execute Python code
			# Note: The code inside py_run_string must be left-aligned.
			# Do not change the indentation, otherwise Python will raise an error.
			py_run_string("
import torch
import os
import pickle
import pandas as pd
import sys
import glob
DevicePath = 'D:/Project/0.MutClone/Function'
sys.path.append(DevicePath+'/2.TransferLearning/dnnTL/0.Code/main')
from ClassifyModel import *
from load_ClassifyModel import load_ClassifyModel
from pathlib import Path


directory = Path(directory)
pkl_files = [p for p in directory.iterdir() if p.is_dir()]


models = {}
for m_file in pkl_files:
	gene_name = os.path.splitext(os.path.basename(m_file))[0]
	if gene_name in target_genes:
		model = load_ClassifyModel(export_dir=m_file)
		model.final_model.GPU_id = None
		# Extract the filename without the extension and prepend 'CRC_' as the key
		model_key = cancer+'_' + os.path.splitext(os.path.basename(m_file))[0]
		models[model_key] = model

print(f'Loaded {len(models)} models in total')
" )
			# Aggregate models
			dnnTL.list = py$models
			cm.bulk = dnnTL.list[[mod.id]]
		}
	
	
		model.list[[paste(mod.id, mod.typ, sep="|")]] = cm.bulk
	}
}


#'-Part [1]-
#' Content: Predict the initial Atlas
#'#######
dat.matrix.knn = seurat.object_pseudo_bulk
pred.res.list = list()
for(mid in names(model.list)){
	
	cm.bulk = model.list[[mid]]
	# 1.2) Cross-platform prediction: normalize using training-data quantiles, then predict
	if( class(cm.bulk)[1]=="list" ){
		source( file.path(Device.path, "1.TraditionalClassifier/0.Code/main/3.ConstructFinalModel.R") )
		pred.res = PredictFinalModel(cm=cm.bulk, predict.data=dat.matrix.knn,
				normalizeUsingCM=TRUE, 
				imputeZero=TRUE)
		pred.tab = data.frame( pred.Prob=pred.res$score, row.names=colnames(dat.matrix.knn))
	}else{
		py$cm_bulk = cm.bulk
		py$test_data = as.data.frame(dat.matrix.knn)
		py_run_string("
import torch
import pandas as pd
import sys
DevicePath = 'D:/Project/0.MutClone/Function'
sys.path.append(DevicePath+'/2.TransferLearning/dnnTL/0.Code/main')
from ClassifyModel import *
pred_value3 = cm_bulk.final_model_Predict(predict_data=test_data, normalizeUsingCM=True, imputeZero=True)
cm_bulk.compute_best_threshold(pos_label=1, best_method='youden')					
" )
		
		pred.res = py$pred_value3
		pred.tab = data.frame( pred.Prob=pred.res$Prediction, row.names=colnames(dat.matrix.knn))
	}
	
	pred.res.list[[mid]] = pred.tab
	cat(mid,"\n")
}

#saveRDS(pred.res.list, file.path(out.dir,"pred.res.list_initAtlas.RDS") )
saveRDS(pred.res.list, file.path(out.dir,"pred.res.list_initAtlas_newCellNames.RDS") )
#pred.res.list = readRDS(file.path(out.dir,"pred.res.list_initAtlas.RDS"))


