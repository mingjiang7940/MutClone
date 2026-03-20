#'##########  ---- Paragraph <3> ----  ###########  
#' 核心目标：  	   构建svmTL模型【无选源, 输出为非概率的决策值】
#' 					           ClassifierSVMTL
#' Part[operate=="construct"] 
#'		Point 输入：
#'				靶数据与标签:
#		 				data    #靶特征矩阵[基因x样本] 
#		 				label   #靶样本标签01向量, names与T.data列名一致
#						model	#来自标准的迁移学习模型, 包含源模型和靶模型, 注意模型data和label必须装配
#
#'		Point 输出:     
#' 				List of 2	#训练好的svmTL模型
#'      [迁移模型]		$ TL.model  :List of 3
#							..$ W           : Named num [1:755] -0.1483 0.0899 -0.1121 -0.1449 0.0909 ...
#							.. ..- attr(*, "names")= chr [1:755] "EDA2R" "KIRREL1" "TEX2" "COL9A1" ...
#							..$ b           : num 2.06
#							..$ use.features: chr [1:755] "EDA2R" "KIRREL1" "TEX2" "COL9A1" ...
#'		[基础模型]		$ base.model:List of 2
#							..$ baselineSVM: num [1:4, 1:756] -0.03694 -0.55644 0 0 0.00314 ...
#							.. ..- attr(*, "dimnames")=List of 2
#							.. .. ..$ : chr [1:4] "CRC" "LUAD" "ESCA" "PAAD"
#							.. .. ..$ : chr [1:756] "EDA2R" "KIRREL1" "TEX2" "COL9A1" ...
#							..$ model.list :List of 1
#							.. ..$ CRC:List of 5
#							.. .. ..$ W           : num [1:755] -0.036943 0.003139 0.000129 0.028322 0.020539 ...
#							.. .. ..$ b           : num 0.28
#							.. .. ..$ use.features: chr [1:755] "EDA2R" "KIRREL1" "TEX2" "COL9A1" ...
#							.. .. ..$ probA       : num -5.48
#							.. .. ..$ probB       : num 0.0894
#' Part[operate=="predict"] 
#' 				Point 输入:
#			 				model   用于预测的模型svmTL模型
#			 				data    用于预测的数据集, [基因x样本] 
#							base.model=FALSE    是否用未迁移的基础模型预测  
#'				Point 输出: 
#							tmp     模型的决策值
#
#'Point 特别注意 尽管函数内部有datatyp只能设置"scRNASeqLogNormCount"或"FC"参数, 但datatyp设置"scRNASeqLogNormCount"依然能适用于bulk表达矩阵
# 函数内部datatyp="scRNASeqLogNormCount"参数是指输入的表达矩阵需要是log后的标化数据, 例如：log(tpm+1), log(cpm+1),
# datatyp="scRNASeqLogNormCount"只作用于alg4_BiasUpdate, alg6_NormalVectorUpdate (alg3_shiftComp输入也出现该参数但不对其有任何影响), 只会影响模型微调时选择的超参数
# 经过实践观察datatyp="scRNASeqLogNormCount"设定测超参数在bulk数据也适用
#'############ [24-04-19] ############
source(file.path(Device.path,"1.TraditionalClassifier/0.Code/1.Extract_Data_and_Labels_from_Model.R"))
source( file.path(Device.path, "1.TraditionalClassifier/0.Code","SVM","linearSVM_Model_Construction_and_Prediction_from_SVMTL.R"))
ClassifierSVMTL <- function(data=NULL, #行为特征（特别注意：已经完成特征选择的）列为样本的矩阵形式
		label=NULL, #标签列1表示pos，0表示阴性
		model=NULL, #来自标准的迁移学习模型, 注意模型data和label必须装配
		operate=c("construct", "predict")[1],
		cross=3,
		...
){			
	#'-Part [1]- 建模
	if(operate == "construct"){
		#使用了原来的源模型的：为了让源模型和靶模型能正常融合，需要让它们特征一致，源模型没有用到的特征等价于模型系数w为0，
		#因此用0填充源模型缺失的特征，保证源模型和靶模型特征一致，后进行源模型和靶模型融合的。
		#'Part [1.1] 训练模型临时函数
		tmpTrain <- function(data=data,label=label,model=model, ... )
		{
		
			#Point [1] 提取数据集
			T.data = list(data) #提取靶数据
			T.label = list(label)
			#if(is.null(model$features)){
				features = rownames(data)
#			}else{
#				features = model$features
#			}
			#最后迁移模型源和靶特征必须一致
			S.data = lapply(model$SourceModels, function(x) x$data[features, ] )   #提取源数据
			S.label = lapply(model$SourceModels, function(x) x$sample.label )
			
			#仅保留源和靶特征矩阵同时存在的 融合特征
			data.list = c(T.data,S.data)
			label.list = c(T.label,S.label)
			
			
			##########Point 1.重构靶域基本SVM模型作为迁移模型基础模型###########
			#' alg1_baselineClass：源域中拟合计算基本的线性SVM分类器的超平面
			#' 结果为(w_1,b_1),(w_2,b_2), ... ,(w_m,b_m)
			#' 存储于alg1_res，作为 Point 2的输入
			######################################
			#(1.1) 重构靶模型
			t.alg1_res = ClassifierSVM.svmTL(data=data, label=label, model.names=model$data.name, ...)  #, ...
			
			#(1.2) 融合源和靶模型参数, 形成迁移初始模型
			#s.alg1_res.list = lapply(model$SourceModels, function(x) x$final.model$model )
			s.alg1_res.list = list()
			for(s.name in names(data.list)[-1]){
				s.alg1_res.list[[s.name]] = ClassifierSVM.svmTL(data=data.list[[s.name]], label=S.label[[s.name]], model.names=s.name, ...)
				#cat(s.name,"\n")
			}
			##特别注意 合并靶模型和预训练好的源模型[子函数]
			MergeSvmModels <- function(t.alg1_res, s.alg1_res.list){
				
				# 1.2.1 特征对齐, 确保源和靶特征权重向量的顺序一致
				alg1_res = t.alg1_res
				t.wb = t.alg1_res$baselineSVM; colnames(t.wb) = c(t.alg1_res$model.list[[1]]$use.features, "b.int")
				s.wb.mtx = matrix(0, length(s.alg1_res.list), length(t.alg1_res$baselineSVM) ); rownames(s.wb.mtx) = names(s.alg1_res.list)
				
				# 1.2.2 填充 源模型 训练好的权重
				for(s.idx in 1:length(s.alg1_res.list ) ){
					s.wb.name = c(s.alg1_res.list[[s.idx]]$model.list[[1]]$use.features, "b.int")
					idx = match(colnames(t.wb), s.wb.name)
					s.wb.mtx[s.idx,] = s.alg1_res.list[[s.idx]]$baselineSVM[1,idx]
				}
				#源模型, 无权重特征, 填充0 #[源模型, 仅有部分融合特征权重, 而无权重基因重要性低, 故设初始系数为0]
				s.wb.mtx[ is.na(s.wb.mtx) ] = 0 
				
				# 1.2.2 合并baselineSVM矩阵
				alg1_res[["baselineSVM"]] = rbind(t.wb, s.wb.mtx)
				return(alg1_res)
			}
			##函数End
			alg1_res = MergeSvmModels(t.alg1_res, s.alg1_res.list)
			
			
			##########Point 2."平均"源域的w和b###########
			#' alg2_rob_meanNCov：“平均”源域模型获得在不同数据域之间不变的参数
			#' 结果为(w_0,b_0)存储于alg2_res，作为 Point 3的输入
			######################################
			suppressMessages(library('GSE'))
			suppressMessages(library('caret'))
			
			source(file.path(Device.path, "2.TransferLearning/svmTL/0.Code","BasicAlg", "Alg2_RobustMeanCov.R"))
			alg2_res <- alg2_rob_meanNCov(alg1_res$baselineSVM)
			
			
			############Point 3.对齐数据集#############
			#' alg3_shiftComp：利用最大互相关在方向上预先对齐源和靶
			#' 		  		   将超平面送到靶域样本类别真实分界的所在的局部
			#' 输出第一步调整好的b存于alg3_res，将作为 Point 4.1的输入
			######################################
			#' 0) 数据预处理：数据调整为RTLbase数据输入格式, 进入后续步骤
			#【训练样本顺序与ClassifierSVM.svmTL内部处理一致, 保证模型输出决策值的正负号一致】
			label.list = lapply(label.list, function(label){
						label[label==0] = -1
						return(label)
					} )
			Ds.X = lapply(data.list, t)
			Ds.Y = label.list
			#' 调整样本顺序 #为保证每个源所计算的SVM决策值正负号的一致性(正类为突变，负类为野生型)，将样本按照lab为1 -> -1排列
			Ds.Y = lapply(Ds.Y,function(x) sort(x, decreasing=TRUE) )
			Ds.X = mapply(function(x,y) x[names(y),], Ds.X, Ds.Y ,SIMPLIFY=FALSE,USE.NAMES = TRUE)
			Dt.X = Ds.X[1]
			
			
			#' 1) 所有的源域数据在基线法线向量w0的方向上都被对齐到靶域数据
			suppressMessages(library('caret'))
			
			source(file.path(Device.path, "2.TransferLearning/svmTL/0.Code","BasicAlg", "Alg3_MaxCrossCor.R"))
			alg3_resCor <- alg3_shiftComp(
					task_list   = Dt.X,  #每个矩阵，行是样本列是特征，
					source_list = Ds.X,
					alg2_result = alg2_res,
					print2screen = FALSE,
					ImpFeats = features,
					save2file = F,
					maximumLag = 0,
					CoreClassifier="LinSVM",
					datatyp="scRNASeqLogNormCount",    #Point 特别注意 该参数不影响此步骤
					useAbsCor = T,
					medianMediansBL = F)
			s.tCorThresh = alg3_resCor$s.tCorThresh    #alg3中最初设为源与靶互相关阈值为0.9,然后寻找对齐位置,若无法对齐到大于该阈值，则逐步按0.1步长放宽阈值，继续对齐，重复直到能对齐的阈值为止
			alg3_res = alg3_resCor$b_updated		   #因此，这里的如果s.tCorThresh值为r代表在该对齐位置源与靶相关程度大于r
			
			
			##########Point 4.1.调整适配靶数据############
			#' alg4_BiasUpdate：①调整b使得超平面穿过靶域密度最稀疏处
			#' 
			#' 输出第二步调整的b存于alg4_res, 将作为 Point 4.2的输入
			######################################
			suppressMessages(library('data.table'))
			suppressMessages(library('caret'))

			source(file.path(Device.path, "2.TransferLearning/svmTL/0.Code","BasicAlg", "Alg4_BiasUpdate.R"))
			alg4_res <- alg4_BiasUpdate(
					task_list = Dt.X,  #每个矩阵，行是样本列是特征，
					alg1_result = alg1_res,
					alg2_result = alg2_res,
					alg3_result = alg3_res,
					goodColumns = features, alg4MinFx = "gd",
					Marg = 1,
					save2file =F, ADM=F,
					useMedian = T, ZnormMappingBL=F, datatyp="scRNASeqLogNormCount",
					RCSmodeBL = F,
					CoreClassifier = "LinSVM")
	
			
			##########Point 4.2.调整适配靶数据############
			#' alg6_NormalVectorUpdate：②调整w使得超平面穿过靶域密度最稀疏处
			#' 输出调整后的w
			######################################
			source(file.path(Device.path, "2.TransferLearning/svmTL/0.Code","BasicAlg", "Alg6_NormVecUpdate.R"))
			alg6_res <- alg6_NormalVectorUpdate(
					task_list = Dt.X,
					alg1_result = alg1_res,
					alg2_result = alg2_res,
					alg3_result = alg3_res,
					alg4_result = alg4_res,
					X_feat_cols = features,
					save2file = F,
					Marg = 1, 
					ADM=F, datatyp="scRNASeqLogNormCount",
					RCSmodeBL = F,
					CoreClassifier = "LinSVM")
			
			
			##########Point 5.最终模型成型############
			#' 将调整好的(W.new,b.Updated)组装成模型
			######################################
			W.new =  alg6_res$alg6_w_new[1,]
			b.Updated = -alg4_res$b_alg_norm
			final.model = list(W=W.new, b=b.Updated, use.features = features)
			
			return(final.model)
		} 
		
		
		#'Part [1.2] 内部交叉验证模型获得Sigmoid参数
		##1.2.1) 划分靶域数据
		suppressMessages(library(caret))
		cvIndex.list = createFolds(factor(label), k=cross, list = TRUE, returnTrain = TRUE)
		##1.2.2) 训练和预测模型
		score.list = list()
		label.list = list()
		for(i in 1:length(cvIndex.list) ){		
			###1.2.2.1) 训练
			train.idx = cvIndex.list[[i]]
			my.model = tmpTrain(data[,train.idx],
								label[train.idx],
								model=model,
								... )
			###1.2.2.2) 预测
			test.idx = setdiff(1:ncol(data), train.idx)
			tmp.data = data[,test.idx]
			y = t( as.matrix(tmp.data[my.model$use.features,]) ) %*% my.model$W + my.model$b
			
			score.list[[i]] = y
			label.list[[i]] = label[test.idx]	
		}	
		##1.2.3) 训练决策值转概率的sigmoid
		suppressMessages(library(platt))
		t.lab = unlist(label.list); t.lab[t.lab==-1] = 0;
		plattScalingResult = plattScaling(unlist(score.list),  t.lab)

		
		#'Part [1.3] 训练模型并添加Sigmoid参数
		final.model = tmpTrain(data,label,model=model, ... )
		final.model$plattScalingResult = plattScalingResult
		
		return(final.model)
	}
	
	
	#'-Part [2]- 预测
	if(operate == "predict"){
		#临时小函数
		RTLmodel.predict <- function(model,Data){
			y.value = t( as.matrix(Data[model$use.features,]) ) %*% model$W + model$b
			return(y.value)
		}
		
		#预测
		tmp = RTLmodel.predict(model, data) #输出决策值
		tmp = platt::predictProb(model$plattScalingResult, tmp)
		return(tmp)
	
	}
}	







