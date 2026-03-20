# TODO: traning SVM
###############################################################################
rm(list = ls() )
Device.path = "D:/Project/0.MutClone/Function"
ID.xteams = c("CRC_TCGA", "BRCA_TCGA", "LUAD_TCGA", "ACC_TCGA", "BLCA_TCGA", "CESC_TCGA", "CHOL_TCGA",
			"DLBC_TCGA", "ESCA_TCGA", "GBM_TCGA", "HNSC_TCGA", "KICH_TCGA", "KIRC_TCGA",
			"KIRP_TCGA", "LGG_TCGA", "LIHC_TCGA", "LUSC_TCGA", "MESO_TCGA",
			"OV_TCGA", "PAAD_TCGA", "PCPG_TCGA", "PRAD_TCGA", "SARC_TCGA", "SKCM_TCGA",
			"STAD_TCGA", "TGCT_TCGA", "THCA_TCGA", "THYM_TCGA", "UCEC_TCGA", "UCS_TCGA",
			"UVM_TCGA")


for(ID.xteam in ID.xteams ){
	
	#ID.xteam = "CRC_TCGA"
	out.dir = file.path(dirname(Device.path),"0.Data/Results", ID.xteam, "1.SVM")
	if(!dir.exists(out.dir)){ dir.create(out.dir, recursive = TRUE) }
	
	way="mut-exprs"
	background=FALSE
	method.select.feature="SVM-RFE"
	
	pancancer.driver = readRDS( file.path(dirname(Device.path), "0.Data", "CancerDriverGene.rds") )
	pancancer.driver = unique(unlist(pancancer.driver))
	file.mut.exprs.data = file.path(dirname(Device.path),"0.Data",'TCGA_MutClone', ID.xteam, "mut.exprs.data.indel.RDS")
	mut.exprs.data = readRDS(file.mut.exprs.data)
	mutgenes = names(mut.exprs.data$all.label)
	mutgenes = mutgenes[mutgenes %in% pancancer.driver]
	
	if(length(mutgenes) == 0){
		cat("No mutation genes in this dataset satisfy the specified criteria.")
	}else{
		
		source( file.path(Device.path,"1.TraditionalClassifier/0.Code/main/1.FeatureSelectionCV.R") )
		source( file.path(Device.path,"1.TraditionalClassifier/0.Code/main/2.OptimizeThresholdCV.R") )
		source( file.path(Device.path,"1.TraditionalClassifier/0.Code/main/ExtractPerformance.R") )
		source( file.path(Device.path,"1.TraditionalClassifier/0.Code/main/3.ConstructFinalModel.R") )
		source( file.path(Device.path,"1.TraditionalClassifier/0.Code/main/4.SimplifyModel.R") )
		out.file.svm = file.path(out.dir, "SVM.rds")
		
		if( file.exists(out.file.svm) ){
			cat("[SVM]")
		}else{
	
			library(parallel)
			result.cm = lapply(mutgenes, function(t.name){
						
						#(1) FeatureSelection
						cm = ClassifierCV.FeatureSelection(data=NULL,
								label=NULL,
								data.name= t.name, #"TP53",
								data.path.rds=file.mut.exprs.data,
								classifier="SVM",
								topn=seq(10, 500, 10),
								times=1)
						


						#(2) get OptimizeThreshold
						cm = ClassifierCV.OptimizeThreshold(cm=cm, method.threshold="Fscore")
						
						#(3) ConstructFinalModel
						cm = ConstructFinalModel(cm=cm, way="mut-exprs", background=FALSE)
						
						cm = SimplifiedModel(cm)
						
						t.file = file.path(out.dir, "单个基因", paste0(t.name, ".rds"))
						if(!dir.exists( dirname(t.file) )){ dir.create(dirname(t.file), recursive = TRUE) }
						
						saveRDS(cm, file=t.file)
						return(cm)
					})
			names(result.cm) = mutgenes
			
			result.cm = result.cm[!sapply(result.cm, is.null)]
			saveRDS(result.cm, file=out.file.svm)
		}
	}
	
	cat(ID.xteam,"\n")
}




