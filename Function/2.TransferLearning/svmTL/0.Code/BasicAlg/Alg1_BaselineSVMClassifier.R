#'#############################  ---- Section <1> ----  ###############################
#' 核心目标:训练基本线性svm函数
#' 对应 流程图中(1)部分
#' 它初始训练出各个源域模型,svm(w_1,b_1),(w_2,b_2),... (w_m,b_m)
#' 注意：这里标签为+1,-1
#'创建日期：2021年12月22日
#'##############################################################################

#' Part [1]
#'内容：核心函数，训练线性svm分类器
alg1_baselineClass <- function(TrainXls,
                               TrainYls,
                               TestXls,
                               TestYls,
                               K_forCrossV = 3,
                               svmGamma=.02,
                               svmCost=0.5,
                               prnt2scr=FALSE,
                               X_cols2Keep=NULL,
                               transX=F, sampleRed=F, doParalellSVM=F, datatyp="FC"){


  #initializations
  errorFlag = F


  #Source task data Dm for m = 1, ...., M datasets

  M = length(TrainXls); if(prnt2scr) print(paste("# of datasets = ", M, sep=""))

  if (length(TrainYls) != M) {
    errorFlag = T
    print("training X and Y lengths dont match")
  }



  #Assuming that all df will have the same # of cols;
  if((class(TrainXls[[1]]) == "numeric")[1]){
    nDims = 1
    if(prnt2scr) print("datasets have 1 dim/feat each")

  } else {
    nDims = ncol(as.data.frame((TrainXls[[1]]))[,X_cols2Keep])
    if(prnt2scr) paste("datasets have ", nDims, " dims/feats each", sep="")
  }


  #this will hold the hyperplane parameters for each SVM
  baselineSVM <- matrix(0, nrow=M, ncol=nDims+1)

  #empty vectors for storage of results
  kfoldSVM_ls    <- vector(mode="list")
  results.all    <- vector(mode="list")

  #-------------------end of inputs and initializations


  if(prnt2scr) print(paste("Starting training with a SVM classifier || cost: ",
                           svmCost, ", gamma: ", svmGamma,
                           ", Kernel: linear" , ", cross: ", K_forCrossV, sep=""))

  #' 记录决策值转概率的sigma函数参数
  probA.list = list()
  probB.list = list()
  for (m in 1:M) {
    #m=1

    if(prnt2scr) print(m)

	####修改日志 - 2021年3月26日
	#' 这应该是错误的
	#' if(is.na(X_cols2Keep[1])) Dm.train  <- as.data.frame(TrainXls[[m]])
	#' if(!is.na(X_cols2Keep[1])) Dm.train  <- as.data.frame(TrainXls[[m]]) 修改为
	if(is.null(X_cols2Keep)) Dm.train  <- as.data.frame(TrainXls[[m]])
	if(!is.null(X_cols2Keep)) Dm.train  <- as.data.frame(TrainXls[[m]])[,X_cols2Keep]
	##

    if(transX){
      Dm.train <- as.data.frame(RTL::AllDataManipulations(Dm.train,
                                                          Scale=MaxPerCh,
                                                          Center=ScalePerCh,
                                                          X_cols2Keep = X_cols2Keep,
                                                          globalASINhTrans = ifelse(datatyp=="FC", T, F),
                                                          globalRange1010Scale = F,
                                                          globalScaleByVal = F))[,X_cols2Keep]
    } else{
		####修改日志 - 2021年3月26日
		#' 这应该是错误的
		#' if(!is.na(X_cols2Keep[1])) { 修改为
		if(!is.null(X_cols2Keep)) {
		##	
        Dm.train <- Dm.train[,X_cols2Keep]
      } else {
        Dm.train <- as.data.frame(Dm.train)
      }
    }


    if(!is.factor(TrainYls[[m]])) {
      TrainYls[[m]] <- factor(TrainYls[[m]])
    }

    if("TRUE" %in% levels(TrainYls[[m]])) TrainYls[[m]] <- factor(ifelse(TrainYls[[m]]==T,1,-1))

    if("Pos" %in% levels(TrainYls[[m]])) POS = c("Neg", "Pos")
    if("Neg" %in% levels(TrainYls[[m]])) POS = c("Neg", "Pos")

    if("-1" %in% levels(TrainYls[[m]]))  POS = c(-1, 1)
    if("1" %in% levels(TrainYls[[m]]))  POS = c(-1, 1)


    if(ncol(Dm.train)==1){
      Dm.train <- as.data.frame(cbind(Dm.train, factor(TrainYls[[m]], levels = POS)))
      colnames(Dm.train) <- c("X", "labelSVM")
    } else {
      Dm.train$labelSVM <- factor(TrainYls[[m]], levels = POS)
    }




    if(nDims == 1) colnames(Dm.train) <- c(paste(rep("X", nDims), 1:nDims, sep=""), "labelSVM")

    if(!(sampleRed==F)){
      if(nrow(Dm.train)>sampleRed){
        Dm.train <- Dm.train[sample(1:nrow(Dm.train), sampleRed),]
      }
    }

    #print(nrow(Dm.train))


    if(doParalellSVM==T){
      Dm.train <<- Dm.train
      #print("intiating paralell SVM")
      #for now running in to memory issues. for future dev
      Dm.train_model <- parallelSVM(labelSVM ~. , data=Dm.train,
                                    cost = svmCost , gamma = svmGamma,
                                    type="C-classification",
                                    kernel = 'linear',
                                    cross = K_forCrossV,
                                    scale = F,
									probability=TRUE, numberCores = 4)
    } else {
      #print("intiating non-paralell SVM")

	#' Point [1]
	#'内容：此处为算法核心部分，主要是调用svm函数对模型进行训练
      Dm.train_model <- svm(labelSVM ~. , data=Dm.train,
                            cost=svmCost, gamma=svmGamma,
                            type="C-classification",
                            kernel = 'linear',
                            cross = K_forCrossV,
                            scale = F,
							probability=TRUE
							)	
    }

#	####修改日志 - 2021年3月26日
#	#' 尝试用自己写的预测,目前仍然发现计算的yhat与预测的标签不一致其中m,已经解决【它是把第一次碰到的类当成正类，其余的为负类】
#	Wm <- drop(t(Dm.train_model$coefs) %*% Dm.train_model$SV)
#	bm <- drop(Dm.train_model$rho)
#	#yhat <- sign(as.matrix(t(Ds[[i]]$Exp)[,use.features]) %*% as.matrix(alg1_res$baselineSVM[i,1:length(use.features)]) - alg1_res$baselineSVM[i,length(use.features)+1])
#	yhat <- sign(-as.matrix(Dm.train[,1:(ncol(Dm.train)-1)]) %*% as.matrix(Wm) + bm)
#	table(Dm.train.pred,yhat)
#	table(Dm.train.pred,Dm.train[,ncol(Dm.train)])
#	##

    #print("predicting on train set")

    Dm.train.pred  <- predict(Dm.train_model, Dm.train)

    if("Pos" %in% levels(Dm.train.pred)) POS = "Pos"
    if("1" %in% levels(Dm.train.pred)) POS = "1"


    results.all[[m]] <- list(train=conf.mat.stats(conf.mat.pred=Dm.train.pred, conf.mat.truth=Dm.train$labelSVM, POS),
                             DmTrainPred = Dm.train.pred,
                             train_model=Dm.train_model)



    if(doParalellSVM){
      #Wm is the hyperplane coeffs
      Wm <- drop(t(Dm.train_model[[1]]$coefs) %*% Dm.train_model[[1]]$SV)
      #rho is the negative intercept
      bm <- drop(Dm.train_model[[1]]$rho)
    } else {
      #Wm is the hyperplane coeffs
      Wm <- drop(t(Dm.train_model$coefs) %*% Dm.train_model$SV)
      #rho is the negative intercept
      bm <- drop(Dm.train_model$rho)

    }

    baselineSVM[m,] <- as.vector(cbind(t(Wm), bm))
	probA.list[[m]] = Dm.train_model$probA
	probB.list[[m]] = Dm.train_model$probB
	
    remove(Dm.train)

  }



  colnames(baselineSVM) <- c(paste(rep("Wm", nDims), 1:nDims, sep=""), "b.int")


  rownames(baselineSVM) <- names(TrainXls)
  names(results.all) <- names(TrainXls)

  ####修改日志 - 2022年8月23日
  #' 改造输出, 和model框架一致
  ##
#  return(list(baselineSVM = baselineSVM,
#              results.all = results.all,
#              M = M,
#              nDims = nDims))
  #'改为
  W_b.list = as.list(as.data.frame(t(baselineSVM)))
  model.list = lapply(W_b.list,function(W_b){
			  tmp.W = W_b[-length(W_b)]
			  tmp.b = -W_b[length(W_b)]
			  use.features = X_cols2Keep
			  return(list(W=tmp.W, b=tmp.b, use.features = use.features))
		  })
  for(m in 1:length(model.list) ){
	  model.list[[m]][["probA"]] = probA.list[[m]]
	  model.list[[m]][["probB"]] = probB.list[[m]]
  }
  return(list(baselineSVM = baselineSVM,
			  model.list = model.list))

}

####修改日志 - 2021年3月31日
#' 由于算法1需要，迁移RTL_genericFXs.R中函数到这里
#' 
##
conf.mat.stats <- function(conf.mat.pred, conf.mat.truth, POS="1"){
	
	#conf.mat.pred   =  factor(sign(y_hat[,x_m]), levels=c(-1,1))
	#conf.mat.truth  =  TASK$TrueClass
	
	temp.confMatdat <- confusionMatrix(data=conf.mat.pred, reference=conf.mat.truth, positive=POS, mode="everything")
	return(list(overall= temp.confMatdat$overall, table=temp.confMatdat$table, byClass=temp.confMatdat$byClass))
	
}











