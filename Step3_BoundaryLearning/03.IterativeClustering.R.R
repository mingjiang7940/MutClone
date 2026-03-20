#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: identify high-confidence cells via iterative clustering:
#' 					we assume that mutant or wild-type cell populations that form clusters in expression space
#' 					have neighboring cells with consistent phenotypes,
#' 					and are more likely to represent true mutant or wild-type cells (high-confidence cells).
#' 
#' Point Input：
#' 				initAtlas				initial Atlas
#' 				keep.mut_samp.list      samples retained after filtering mutant samples whose scores are similar to wild-type cells
#' Point Output：
#' 				Atlas.RDS				reference Atlas after iterative clustering
#'############ [25-01-07] ############
rm(list=ls())
Device.path = "D:/Project/0.MutClone/Function"
out.dir = "/WorkSpace/chengmingjiang/1.Project_sc2mutTL/3.Section/00.完成故事梳理后_整理版/S3_模型偶联Atlas/0.Data/新版本/03.迭代聚类寻找高置信的细胞簇"
cancer = "CRC"
if (!dir.exists(out.dir))	dir.create(out.dir, recursive = TRUE)
mut.ls = c("APC","KRAS","TP53","BRAF","PIK3CA","SMAD4","FBXW7") #input mutation models to be used


#data loading
iniAtlas = readRDS( file.path(dirname(out.dir),"initAtlas.RDS") )  #select input
keep.mut_samp.list  = readRDS(file.path(dirname(out.dir), "02.过滤突变样本使样本标签近似细胞标签", "keep.mut_samp.list.RDS"))


#'-Part [2]-
#' Content: iterative clustering of Atlas
#'#######
#set error tolerance to 20%
#options(error = recover)
options(error = NULL)
seurat.object0 = iniAtlas
for(i in 1:length(mut.ls) ){ 
	
	#Point [0] extract information
	mutation = mut.ls[i]
	mod.id = paste("CRC",mutation,sep="_")
	mod.typs = c("svm", "svmTL", "dnnTL")
	seurat.object = seurat.object0
	Mutation.info = seurat.object@meta.data[, paste0(mutation,"Mutation")]
	seurat.object@meta.data[,"True.lab"] = Mutation.info
	
	#obtain sample information
	md = seurat.object@meta.data[,c("True.lab","SampleID")]
	is.NA = unique(md$SampleID[is.na(md$True.lab)])
	md = md[!is.na(md$True.lab),]
	samples_by_label <- list(
			MT = unique(md$SampleID[md$True.lab == 1]),
			WT = unique(md$SampleID[md$True.lab == 0]),
			is.NA = is.NA
	)
	
	
	for(mod.typ in mod.typs){
		modid_type = paste(mod.id, mod.typ, sep="|")
		
		#Point [1] filter mixed mutant samples
		t.samples_by_label = samples_by_label
		t.samples_by_label$MT = keep.mut_samp.list[[modid_type]]
		if(length(t.samples_by_label$MT)==0) next
		keep.samps = unlist(t.samples_by_label)
		t.seurat.object = subset(seurat.object, SampleID %in% keep.samps)
		
		#Point [2] train Atlas
		DefaultAssay(t.seurat.object) <- "integrated" #must use integrated data for clustering
		source( file.path(Device.path,"3.BLADE/0.Code/main/IterativeClustering_LabelEnrichedCellPopulations.R") )
		atlas.res = IterClusteringGetLabEnrichedClusters(seurat.object=t.seurat.object, 
				major.clust.prop = 0.95,       #cluster consistency proportion
				init.resolution=0.1)
		
		#save results
		saveRDS(atlas.res, file =  file.path(out.dir,paste0(modid_type,"Atlas.RDS") ))
	}
	
	cat(i,"/",length(mut.ls),"\n")
}

