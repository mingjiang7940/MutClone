#'#############################  ---- Section <1> ----  ###############################
#'核心目标：平均源域模型
#' 对应流程图PPT中(2)平均源域模型 部分
#' 这里对训练来的模型𝑆𝑉𝑀(w_1,𝑏_1)，(w_2,𝑏_2)，... (w_m,𝑏_m)参数进行“平均”
#'创建日期：2021年12月22日
#'##############################################################################


#Mahyari, Eisa: Implementation of Lee et al.'s TL-FC into R
#Feb, 2017


##### Robust Mean and Covariance (Algorithm 2) ##########
###
### Input: (W_m, b_m) for m = 1, ..., M
###    Concatenate: u_m <- [W_m, b_m], for all m
###    Initialize: mu <- mean(u_m), C <- cov(u_m)
###    repeat:
###       d_m <- sqrt((u_m - mu)^T C^-1 (u_m-mu))
###       w_m <- phi(d_m)/d_m
###       Update:
###          mu_new <- sum_m(w_m U_m)/sum_m(w_m)
###          C_new <- (sum_m(w_m^2(u_m - mu_new)(u_m - mu_new)^T))/sum_m(w_m^2-1)
###    until: stopping conditions are satisfied
### Output: mu = [w_0, b_0], C_0 = C(1:d, 1:d)
###
############################################################
### Lee, Gyemin, L Stoolman, and C Scott.
### "Transfer Learning for Auto-Gating of Flow Cytometry Data.”
### JMLR(workshop), 2012, 155–66.
############################################################
###
###   Using the GSE package's HuberPairwise function,
###    we estimate the robust mean and covariance.
###    A tuning constant c0 = 1.345 is used in the huber pairwise estimation function
###
############################################################
###  Alqallaf, F.A., Konis, K. P., R. Martin, D., Zamar, R. H.
###  "Scalable Robust Covariance and Correlation Estimates for Data Mining."
###  In Proceedings of the Seventh ACM SIGKDD International Conference
###     on Knowledge Discovery and Data Mining. Edmonton. 2002.
############################################################


#########
######
###
#' Part [1]
#'内容：核心函数，平均模型参数函数
alg2_rob_meanNCov <- function(alg1_result_baselineSVM, print2screen = F){
  #alg1_result_baselineSVM <- alg1_res$baselineSVM; print2screen = T


  if(is.null(alg1_result_baselineSVM)) print("alg1_result_baselineSVM is Null") else {

    if(nrow(alg1_result_baselineSVM)>1) {

      if(print2screen) print("starting Alg 2 to obtain robust mean and covariance of SVM hyperplane parameters")


      #baselineSVM is a matrix of concatination of the normal vector and the y-intercept (bias)
      baselineSVM <- alg1_result_baselineSVM



      #INITIALIZE: C is covariance and U is mean

      C <- cov(baselineSVM); C
      try.U <- try(CovEM(baselineSVM), silent = T)
      if(!(class(try.U)=="try-error")) {
        res <- CovEM(baselineSVM) #Gaussian MLE of mean and covariance
        U <- getLocation(res); #equal to getting mean as U <- sapply(baselineSVM, mean); U
        names(U) <- colnames(C); U
        colnames(C) <- names(U)
        rownames(C) <- names(U)
      } else {
        U <- as.vector(mode="numeric", colMeans(baselineSVM))
        names(U) <- colnames(C); U
      }

      #C <- getScatter(res);C # different to cov()


      # respmh <- GSE::partial.mahalanobis(baselineSVM, mu=U, S=C)
      # plot(respmh, which="index")
      # getDist(respmh)

      #huber pairwise estimation; using prev. built for now
      #c0 is the tuning constant for the huber function.
      #c0=0 would yield QC. Default is c0=1.345

	  #' Point [1]
	  #'内容：算法核心部分
      #' resHub <-try(HuberPairwise(as.matrix(baselineSVM), psi=c("huber"), c0=1.345, computePmd=TRUE), silent = T)
	  ####Point 修改a - 2025年9月2日: 增加了迭代过程, HuberPairwise函数不含迭代
	  tmp.list = as.list( t(baselineSVM) ) #将每个模型的参数分别装进list中
	  resHub = try(robust_model_average(tmp.list, delta = 1.345, tol = 1e-6, max_iter = 100), silent = T)
	
      if(!(class(try.U)=="try-error")) {
        #getDist(resHub)
        #plot(resHub, which="index")

        # resHub@R #correlation matrix
        # resHub@mu #the mean weighted by the Huber loss function, a robust loss fx
        # resHub@pmd #partial mahalanobis
        # resHub@S #the cov weighted by the Huber loss function, a robust loss fx
        # resHub@x


        #OUTPUT:
        #U #the simple mean
        #'U_robust <- resHub@mu; U_robust
		U_robust <- resHub$mu    #Point 修改a：对接robust_model_average
	
        #C #the simple cov
        #'C_robust <- resHub@S; colnames(C_robust) <- colnames(C);
		C_robust <- resHub$C; colnames(C_robust) <- colnames(C);   #Point 修改a：对接robust_model_average
        rownames(C_robust) <- colnames(C); C_robust

        remove(res)

      } else {

        #OUTPUT:
        #U #the simple mean
        U_robust <- U
        #C #the simple cov
        C_robust <- C;

      }



#' 高维数据会很慢
      alg2_CalcMore_res <- alg2_CalcMore(list(U_simple=U, U_robust=U_robust, C_simple=C, C_robust=C_robust))

      names(alg2_CalcMore_res$U_robust_norm) <- names(baselineSVM)

      if(print2screen) print("Alg 2 has completed!")
    }

    if(nrow(alg1_result_baselineSVM) == 1) {





    }






    return(alg2_CalcMore_res)

  }

}


#' Part [2]
#'内容：此处涉及到后续选择超平面部分，给超平面提供旋转方向v_t_norm
alg2_CalcMore <- function(alg2_res) {

  #alg2_res = list(U_simple=U, U_robust=U_robust, C_simple=C, C_robust=C_robust)
  tempC <- alg2_res$C_robust[1:(nrow(alg2_res$C_robust)-1),1:(ncol(alg2_res$C_robust)-1)]

	####修改日志 - 2021年6月23日
	#' alg2_res$v_0 <- eigen(tempC)$values	#'alg2_res$v_0 <- eigen(tempC,only.values=T)$values 这能否加快速度
	#' 修改为，加快速度
	  t1=proc.time()
  alg2_res$v_0 <- eigen(tempC,only.values=T)$values
	  t2=proc.time()
	  t=t2-t1
	  #print(paste0('执行时间：',t[3][[1]]/60,'分钟'))
  
  
  #alg2_res$U_robust # <w_0, b_0>
  alg2_res$w_euc_mag <- norm(alg2_res$U_robust[-length(alg2_res$U_robust)], type="2") #euclidean norm

  #w_t,b_t
  alg2_res$U_robust_norm <- alg2_res$U_robust/alg2_res$w_euc_mag


  #'添加
  t1=proc.time()
  #orthonormalize v_0 with respect to w_0
  alg2_res$v_t <- orthnormal(rbind(as.numeric(alg2_res$v_0[]),
                                   as.numeric(alg2_res$U_robust_norm[-length(alg2_res$U_robust_norm)])))
  #'添加
   t2=proc.time()
   t=t2-t1
   #print(paste0('执行时间：',t[3][[1]]/60,'分钟'))
		

  alg2_res$v_t_norm <- alg2_res$v_t[1]#/alg2_res$w_euc_mag

  return(alg2_res)
}


####修改日志 - 2021年3月31日
#' 从RTL_SecondaryFXs.R 搬到这
#' 
##

# Gram-Schmidt orthonormalization, vectors in rows
# implemented from an online matlab version of the code

orthnormal <- function(X)
{
	#X <- cbind(as.numeric(alg2_res$v_0[-length(alg2_res$v_0)]), as.numeric(alg2_res$U_robust[-length(alg2_res$U_robust)]))
	X<-as.matrix(X)
	n<-nrow(X)
	p<-ncol(X)
	
	W<-NULL
	if(p > 1) {
		W<-cbind(W, X[,1])
		for(k in 2:p) {
			#k = 2
			gw<-rep(0, n)
			
			for(i in 1:(k-1)) {
				#i=1
				gki<-as.vector((t(W[,i]) %*% X[,k])/(t(W[,i]) %*% W[,i]))
				gw<-gw + gki * W[,i]
			}
			W<-cbind(W, X[,k] - gw)
		}
	} else {
		W<-cbind(W, X[,1])
	}
	
	W <- replace(W, is.na(W), 0)
	
	####修改日志 - 2021年4月6日
	#' 修改
	#' W <- - apply(W, 2, function(x) norm(x, type="2"))
	#' 为
	W <- apply(W, 2, function(x) norm(x, type="2"))
	##
	
	W
}





########Point 修改a - 2025年9月2日: 增加了迭代过程, 此前HuberPairwise函数不含迭代########
#' -----------FUN: GetRobustAvg-----------
robust_model_average <- function(mu_list, delta = 1.345, tol = 1e-6, max_iter = 100){
	# mu_list: list of model parameter vectors, each mu_m = c(w_m, b_m)
	# delta: Huber loss threshold
	# tol: 收敛阈值
	# max_iter: 最大迭代次数
	M <- length(mu_list)
	p <- length(mu_list[[1]])
	
	# (1) 初始化均值向量 mu 和协方差矩阵 C
	mu <- Reduce("+", mu_list) / M
	mat <- do.call(rbind, mu_list)
	C <- cov(mat)
	
	for(iter in 1:max_iter){
		# 1) 计算每个模型相对均值的平方马氏距离
		d_sq <- sapply(1:M, function(m){
					diff <- mu_list[[m]] - mu
					as.numeric(t(diff) %*% solve(C) %*% diff)
				})
		d <- sqrt(d_sq)
		
		# 2) 基于 Huber ψ 函数计算权重
		psi <- ifelse(d <= delta, d, delta)
		weights <- psi / d
		weights[d == 0] <- 1  # 避免除以0
		
		# 3) 更新均值向量 mu_new
		mu_new <- Reduce("+", lapply(1:M, function(m) weights[m] * mu_list[[m]])) / sum(weights)
		
		# 4) 更新协方差矩阵 C_new
		C_new <- Reduce("+", lapply(1:M, function(m) weights[m]^2 * 
											(mu_list[[m]] - mu_new) %*% t(mu_list[[m]] - mu_new))) / 
				(sum(weights^2) - 1)
		
		# 5) 检查收敛
		if(max(abs(mu_new - mu)) < tol){
			mu <- mu_new
			C <- C_new
			break
		}
		
		mu <- mu_new
		C <- C_new
	}
	
	return(list(mu = mu, C = C, weights = weights, iterations = iter))
}

#' -----------FUN: GetRobustAvg-----------
#####################
if(FALSE){
	
	# 假设有三个模型，每个模型参数向量是 w,b
	mu_list <- list(
			c(0.5, 1.2),
			c(0.6, 1.1),
			c(10, 20)   # 一个异常值
	)
	
	res <- robust_model_average(mu_list)
	res$mu       # 稳健平均的参数向量 (w0, b0)
	res$C        # 对应协方差矩阵
	res$weights  # 每个模型的权重
	res$iterations
}


