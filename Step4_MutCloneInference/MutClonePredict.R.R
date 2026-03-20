#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective:    MutClone prediction on single-cell data
#' 
#'-Part [0]- load: datasets to be predicted and perform smoothing with 50 nearest neighbors
#'-Part [1]- load: models
#'-Part [2]- use models to predict new single-cell datasets
#' 
#' 
#' Point Input：
#' 				seurat.object.list				datasets to be predicted, log2(tpm+1) or log2(cpm+1)
#' 				model.list						models	
#' 
#' 				select.model					selected models
#' 				threshold_down_and_up.list		model-specific confidence intervals
#' 				
#' 
#' Point Output：
#' 				pred.dataset.res.list.RDS  stores prediction scores
#' 				pred.dataset.res.list_predLab.RDS  stores prediction scores + predicted labels
#'############ [25-01-07] ############
rm(list=ls())
Device.path = "D:/Project/0.MutClone/Function"
out.dir = "/WorkSpace/chengmingjiang/1.Project_sc2mutTL/3.Section/00.完成故事梳理后_整理版/S3_模型偶联Atlas/0.Data/新版本/05.模型预测外部单细胞独立数据集"
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
cancer = "CRC"


#'-Part [0]-
#' load: datasets to be predicted and perform smoothing with 50 nearest neighbors
#' 		 expression as log2(cpm+1) or log2(tpm+1)
#'#######
#' from "00.external independent dataset preprocessing"
#' t.dir = "/WorkSpace/chengmingjiang/TmpData/integrateData/2.Atlas偶联迁移模型/模型预测独立测试集/1.多区域WES"
t.dir = file.path(dirname(out.dir),"00.外部独立数据集预处理")
seurat.object.list = readRDS( file.path(t.dir,"seurat.object.list.RDS") )

#' smooth data; patients with fewer than 50 cells will not be smoothed
t.dir = file.path(dirname(out.dir),"00.外部独立数据集预处理")
seurat.object_pseudo_bulk.list = list()
for(dataset in names(seurat.object.list) ){
	library(Seurat)
	seurat.object = seurat.object.list[[dataset]]
	if(dataset %in% c("Li_Gut_2020") ) seurat.object[["PatientID"]] = seurat.object[["SampleID"]]
	
	source( file.path(Device.path, "3.BLADE/0.Code/ConstructPseudoBulk_IntraPatient.R") )
	seurat.object_pseudo_bulk.list[[dataset]] = psudoBulkInPat(seurat.object, k=50, mc.cores=1)   #Point 输入输出数据的单位均为log2(cpm+1)或者log2(tpm+1)
	cat(dataset,"\n")
}
#saveRDS(seurat.object_pseudo_bulk.list, file.path(t.dir,"seurat.object_pseudo_bulk.list.RDS") )
#seurat.object_pseudo_bulk.list = readRDS(file.path(t.dir,"seurat.object_pseudo_bulk.list.RDS"))



#'-Part [1]-
#' Content: perform score prediction
#'#######
#(1) load corresponding models
mut.ls = c("APC","KRAS","TP53","BRAF","PIK3CA","SMAD4","FBXW7") #输入要使用的突变模型
model.list = list()
for(i in 1:length(mut.ls) ){ 
	
	mutation = mut.ls[i]
	mod.id = paste("CRC",mutation,sep="_")
	mod.typs = c("svm", "svmTL", "dnnTL")
	for(mod.typ in mod.typs){
		
		#load bulk models
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
			use_condaenv("/Software/PythonENV/DeepLearning.py3.9", required = TRUE)  # 指定 conda 环境
			py_config()
			DevicePath = '/WorkSpace/chengmingjiang/1.Project_sc2mutTL/BladePipline/0.Code'# 检查 Python 是否已经正确设置
			source_python( file.path(DevicePath,"2.TransferLearning/dnnTL/0.Code/main", "ClassifyModel.py") ) #	
			#py$directory = file.path("/WorkSpace/chengmingjiang/TmpData/TL_model",paste0(cancer,"_TCGA"),"2.迁移学习dnnTL/单个基因3/")
			py$directory = file.path("/IData/DataCenter/TCGA",paste0(cancer,"_TCGA"),"Results/BioGenomics/01.突变/93.预测突变/11.基于表达预测突变/120.基因突变谱与表达谱[带分类标签]/3.dnnTL/单个基因/")
			py$cancer = cancer
			py$target_genes = mutation
			
			# use reticulate py_run_string to execute Python code 
			#Point code inside py_run_string must start at the left margin; do not change indentation, otherwise Python will fail
			py_run_string("
import torch
import os
import pickle
import pandas as pd
import sys
import glob
DevicePath = '/WorkSpace/chengmingjiang/1.Project_sc2mutTL/BladePipline/0.Code'
sys.path.append(DevicePath+'/2.TransferLearning/dnnTL/0.Code/main')
from ClassifyModel import *
from load_ClassifyModel import load_ClassifyModel
from pathlib import Path
							
directory = Path(directory)
pkl_files = [p for p in directory.iterdir() if p.is_dir()]
							
# specify gene models to load						
models = {}
for m_file in pkl_files:
	gene_name = os.path.splitext(os.path.basename(m_file))[0]
	if gene_name in target_genes:
		model = load_ClassifyModel(export_dir=m_file)
		model.final_model.GPU_id = None
		model_key = cancer+'_' + os.path.splitext(os.path.basename(m_file))[0]
		models[model_key] = model

print(f'总共加载了 {len(models)} 个模型')
" )
			#汇总模型
			dnnTL.list = py$models
			cm.bulk = dnnTL.list[[mod.id]]
		}	
		model.list[[paste(mod.id, mod.typ, sep="|")]] = cm.bulk
	}
}


	
	
#'-Part [2]-
#' Content: model prediction on new single-cell datasets
#'#######
#(2) prediction
library(Seurat)
pred.dataset.res.list = list()
for(dataset in names(seurat.object.list)){
	
	dat.matrix.knn = seurat.object_pseudo_bulk.list[[dataset]]
	pred.res.list = list()
	for(mid in names(model.list)){
		
		cm = model.list[[mid]]
		# 1.2) cross-platform prediction: normalize using training data quantiles, then predict pseudo-bulk
		if( class(cm)[1]=="list" ){ #standard models, not dnnTL
			source( file.path(Device.path, "1.TraditionalClassifier/0.Code/main/3.ConstructFinalModel.R") )
			pred.res = PredictFinalModel(cm=cm, predict.data=dat.matrix.knn,
					normalizeUsingCM=TRUE, #
					imputeZero=TRUE)
			pred.tab = data.frame( pred.Prob=pred.res$score, row.names=colnames(dat.matrix.knn))
		}else{
			py$cm = cm
			py$test_data = as.data.frame(dat.matrix.knn)
			py_run_string("
import torch
import pandas as pd
import sys
DevicePath = '/WorkSpace/chengmingjiang/1.Project_sc2mutTL/BladePipline/0.Code'
sys.path.append('/WorkSpace/chengmingjiang/1.Project_sc2mutTL/BladePipline/0.Code/2.TransferLearning/dnnTL/0.Code/main')
from ClassifyModel import *
pred_value3 = cm.final_model_Predict(predict_data=test_data, normalizeUsingCM=True, imputeZero=True)					
" )
			
			pred.res = py$pred_value3
			pred.tab = data.frame( pred.Prob=pred.res$Prediction, row.names=colnames(dat.matrix.knn))
		}
		pred.res.list[[mid]] = pred.tab

		cat(mid,"\n")
	}
	
	pred.dataset.res.list[[dataset]] = pred.res.list
	cat(dataset,"\n")
}

path = file.path(out.dir,"RData",paste0("pred.dataset.res.list.RDS") )
dir.create(dirname(path),recursive = TRUE, showWarnings = FALSE)
saveRDS(pred.dataset.res.list, path)




#'-Part [3]-
#' Content: convert prediction scores to labels based on confidence boundaries
#'#######
pred.dataset.res.list0 = readRDS( file.path(out.dir,"RData",paste0("pred.dataset.res.list.RDS") ) )

#load selected model information
select.model = readRDS(file.path(dirname(out.dir), "select.model.RDS") )
select.model$mod.id = paste(select.model$Prediction, select.model$Algorithm, sep="|")

library(Seurat)
pred.dataset.res.list = list()
for(dataset in names(pred.dataset.res.list0)){
	pred.res.list = pred.dataset.res.list0[[dataset]]
	
	for(mid in select.model$mod.id ){
		
		#(1) model confidence intervals
		pred.tab = pred.res.list[[mid]]
		idx = select.model$mod.id==mid
		x = select.model[idx,"param.x"]
		y = select.model[idx,"param.y"]
		threshold_down_and_up.list = readRDS(file.path(dirname(out.dir),"04.模型置信边界确定", paste0("threshold_down_and_up.list[",x,"_",y,"].rds")  ))
		
		#(2) assign labels based on confidence intervals
		threshold_down_and_up = threshold_down_and_up.list[[paste0(mid,"Atlas")]]
		pred.tab$pred.Lab <- NA_integer_  
		pred.tab$pred.Lab[pred.tab$pred.Prob <= threshold_down_and_up["threshold.down"]] <- 0
		pred.tab$pred.Lab[pred.tab$pred.Prob >= threshold_down_and_up["threshold.up"]] <- 1
		
		pred.dataset.res.list[[dataset]][[mid]] = pred.tab
	}
	
}


path = file.path(out.dir,"RData",paste0("pred.dataset.res.list_predLab.RDS") )
dir.create(dirname(path),recursive = TRUE, showWarnings = FALSE)
saveRDS(pred.dataset.res.list, path)






