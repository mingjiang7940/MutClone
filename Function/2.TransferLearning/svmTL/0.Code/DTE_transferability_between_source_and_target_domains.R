#'##########  ---- Paragraph <1> ----  ###########  
#' 核心目标：  计算类间离散矩阵Sbetween
#'     该指标衡量矩阵数据中特征可以区分样本不同类别的能力
#' 以下Paragraph<2> 和Paragraph <3> 使用该指标完成迁移学习中的源域选择
#' 
#' Point 输入：
#' 				dat.matrix 	数据特征矩阵, 行为特征列为样本
#' 				lab		   	样本的标签向量
#' 
#' Point 输出： 
#' 				Sbetween 	数据的类间离散矩阵
#'############ [22-07-24] ############
CalculateScatterMatrix <-function(dat.matrix, lab){
	
	#dat.matrix = source.data
	#lab=source.lab
	#dat.matrix=tmp.data
	#lab=tmp.lab
	#计算各个类别的中心点
	tmp.dat = cbind(data.frame(lab=lab),t(as.data.frame(dat.matrix)) )
	mean.each.class = aggregate(tmp.dat[,-1],by=list(lab),FUN=mean,na.rm=TRUE)
	mean.each.class.matrix = t(mean.each.class[,-1])
	
	#计算数据总中心
	mean.total = rowMeans(dat.matrix)
	mean.total.matrix = matrix(rep(mean.total,ncol(mean.each.class.matrix) ),  dim(mean.each.class.matrix))
	
	#计算各个类别样本数
	samp.num.each.class= table(lab)[as.character( mean.each.class[,1] ) ]
	
	#计算散度矩阵
	tpm.deviation.matrix = mean.each.class.matrix-mean.total.matrix
	#每一类的偏差乘其转置，计算每类的Sb矩阵
	Sb.matrix.each.class = lapply(as.list(as.data.frame(tpm.deviation.matrix)),function(X) X %*% t(X) )  #
	Sb.matrix.each.class = mapply(function(x,y) x*y,
			as.list(samp.num.each.class),
			Sb.matrix.each.class,
			SIMPLIFY = FALSE) #各类别的矩阵按样本数目加权
	
	#叠加各个类的类别间离散矩阵矩阵
	Sbetween = Reduce('+' , Sb.matrix.each.class)
	
	return(Sbetween)
}


#########2022年7月24日########
##' -----------测试案例-----------
#dat.matrix = tcga.exp_mut[[3]]$Exp
#lab = tcga.exp_mut[[3]]$Mut["TP53",]
##' -----------结束测试-----------
######################






#'##########  ---- Paragraph <2> ----  ###########  
#' 核心目标：计算两个数据域(矩阵)的可迁移性指标DTE(Domain Transferability Estimation)
#' 						的分子和分母, 这里分别称为DIS和DIF
#' 
#' Point 输入：
#' 				source.data 单个源域的矩阵, 行是基因列是样本(行名必须是基因名)
#' 				source.lab  源域的标签向量, 顺序与source.data矩阵的样本顺序一致
#' 				target.data 靶数据的矩阵, 行是基因列是样本(行名必须是基因名)
#' 
#' Point 输出：	DTE的分子和分母两个数值组成的向量
# 					c(DIS,DIF)  
#								DIS为DTE的分子, DIF为DTE的分母
# 								Point 注意，这里不能直接使用DIS/DIF 生成最终DTE的值来排序各源, 因为DIS和DIF需要归一化后才能比较
# 								Point 推荐：把各源域DIS归一化到0-1, DIF归一化到1-0[相当于取倒数], 然后二者相乘获得DTE的值[此时DTE的范围在0-1], 
# 									  此时才能排序源域，然后选择DTE较大的源域进入后续训练
#'									  因此，这里不返回DIS/DIF的值选源, 而是在此函数之外归一化以后才能排序选源
#' 
#' 代码参考来源 https://github.com/chamwen/MEKT/blob/master/lib/DTE.m
#' https://github.com/chamwen/MEKT/blob/master/dte_rsvp.m 【79行-89行】
#' #######################
CalculateDIS_DIF <- function(source.data,source.lab,
		target.data)
{
	#source.data=source.data.list[[1]]
	#source.lab=source.lab.list[[1]]
	#target.data=Tdata
	#修整数据
	commgene = intersect(rownames(source.data) , rownames(target.data))
	source.data = source.data[commgene,]
	target.data = target.data[commgene,]
	
	#计算源域类别间离散矩阵
	#cat("计算类别间离散矩阵\n")
	Sbs = CalculateScatterMatrix(source.data, source.lab)
	
	#计算源与靶间离散矩阵
	#cat("计算源与靶间离散矩阵\n")
	tmp.data = cbind(source.data,target.data)
	tmp.lab = rep(c("s","t"), times=c(ncol(source.data),ncol(target.data)))
	Sbst = CalculateScatterMatrix(tmp.data,tmp.lab)   #Paragraph <1>的函数

	
	#计算Domain Transferability Estimation
	#cat("计算DTE中的DIS和DIF\n")
	DIS = norm(Sbs,type="O")
	DIF = norm(Sbst,type="O")	#离散矩阵一阶范数
	
	return( c(DIS=DIS,DIF=DIF) )	#THCA与SKCM计算的DET=0.05
	
}

#CalculateDET(source.data,source.lab,target.data,target.lab)








#'##########  ---- Paragraph <3> ----  ###########  
#' 核心目标：使用DTE(Domain Transferability Estimation)选源
#' 			   按DTE得分对各源域的可迁移性由强到弱进行排序
#' 
#' Point 输入：
#' 				source.data.list  候选的源域list，注意list的name为各源域的名称
#' 				source.lab.list	  与源域数据对应的标签list, Point 注意，list里的lab每个元素可以没有名字，如果有名字必须与样本名一致
#' 				target.data       靶域数据
#' 
#' Point 输出：
#' 				DTE  			  按值DTE得分由大到小排序的DTE得分向量，向量名称为对应源域名称
#' 
# 
#' 代码参考来源 https://github.com/chamwen/MEKT/blob/master/lib/DTE.m
#' https://github.com/chamwen/MEKT/blob/master/dte_rsvp.m 【79行-89行】
#'############ [22-07-24] ############
SortSourceByDTE <- function(source.data.list,   #源域数据list
		source.lab.list, 	#与源域数据对应的标签list
		target.data){		#靶域数
		
	#source.data.list=Sdata.list
	#source.lab.list=Slab.list
	#target.data=Tdata
	DIS_DIF.list = mapply(function(s.dat,s.lab) CalculateDIS_DIF(s.dat,s.lab,target.data), #Paragraph <2>的函数
			source.data.list,
			source.lab.list,
			SIMPLIFY = FALSE) #各类别的矩阵按样本数目加权
	names(DIS_DIF.list) = names(source.data.list)
	
	DIS.ls = unlist( lapply(DIS_DIF.list,function(X) X[1] ) )
	DIF.ls = unlist( lapply(DIS_DIF.list,function(X) X[2] ) )
	
	#归一化
	library("scales")
	DIS.ls = rescale(DIS.ls,c(0,1))   #此函数相当于源代码rk(2,:)=mapminmax(rk(2,:),0,1); 86行
	DIF.ls = rescale(DIF.ls,c(1,0))	  #Point 此处归一化，有取1/DIF的效果
	
	DTE = DIS.ls*DIF.ls; names(DTE) = names(DIS_DIF.list);
	
	DTE = sort(DTE,decreasing=TRUE)
	return(DTE)
	
}





#########2022年7月25日########
##' -----------测试案例-----------
##' 测试源代码涉及到的数据看结果是否一致
#setwd("/pub6/temp/mingjiang/RTL/DTE")
#library(openxlsx)
#library(tidyverse)
##获得一个excel的多个工作本, 各个源域的特征矩阵
#name_sheet = paste0("Sdata",1:10)
#Sdata.list = map(name_sheet, ~ read.xlsx("Sdata.xlsx",sheet=.,colNames = FALSE))
#names(Sdata.list)=name_sheet
#
##获得各源域的标签
#name_sheet = paste0("Slab",1:10)
#Slab.list = map(name_sheet, ~ read.xlsx("Slab.xlsx",sheet=.,colNames = FALSE))
#Slab.list = lapply(Slab.list, function(X) X[,1] )
#names(Slab.list)=names(Sdata.list)
#
#
##获得靶域数据
#Tdata = read.xlsx("Tdata.xlsx",sheet="Tdata",colNames = FALSE)
#res.dte.score = SortSoueceByDTE(Sdata.list,Slab.list,Tdata)
#
##' -----------结束测试-----------
######################
