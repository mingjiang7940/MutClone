#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: determine a dual-threshold prediction boundary based on score and label; retain predictions on both sides of the thresholds and reject predictions in the middle region
#'            On the ROC curve:
#'							find the threshold at 95% sensitivity as the lower bound
#'							find the threshold at 95% specificity as the upper bound
#' Point Input:
#' 			  score  model output scores
#'            label  true 0/1 labels
#' 
#' Point Output：
#' 			  list(threshold_down_and_up,  #lower and upper boundary thresholds
#' 					p.plot=p	)          #ggplot data object
#' 
#'          		 where threshold_down_and_up includes:
#'			 				threshold_down_and_up <- c(
#														threshold.down = threshold.down,    #lower threshold
#														threshold.up   = threshold.up,		#upper threshold
#
#														AUROC.full = AUROC.full,            	#AUROC without considering the lower/upper thresholds, representing the ranking ability of the model score
#														Intermediate.porp = Intermediate.prop, 	#proportion of the middle region among all samples
#
#														#Point Special note: these metrics are all calculated after excluding the middle region
#														PPV = PPV,								#positive predictive value
#														NPV = NPV,								#negative predictive value
#														ACC = ACC								#accuracy negative predictive value
#												)
#' 
#' ref：
#' 		https://github.com/noellewarmenhoven/Lumipulse-Nat-Med/blob/main/Calculating%20thresholds.R
#'    	Plasma phospho-tau217 for Alzheimer’s disease diagnosis in primary and secondary care using a fully automated platform, NatureMedicine, 2025
#'############ [26-01-15] ############
RocGetTwoCutoffs <- function(score, label, pos=c("1"),
		sensitivity = 0.95,
		specificity = 0.95,
		n.boot = 1000,
		
		#plotting parameters
		main = NULL,
		need.log = FALSE,
		base_size = 24) {
	
	library(boot)
	library(pROC)
	library(tidyverse)
	library(dplyr)
	label = as.character(label)
	#'Part [0] Calculate the AUROC before considering the boundaries
	roc.full <- roc(
			response  = label,
			predictor = score,
			levels    = c(setdiff(unique(label), pos), pos),
			direction = "<",
			quiet     = TRUE
	)
	AUROC.full <- as.numeric(auc(roc.full))
	
	
	#'-Part [1]-
	#' Content: calculate confidence boundaries and related performance metrics
	#'#######
	if (length(table(label)) == 2 && min(table(label)) >= 10) {
		
		dat <- data.frame(score = score, label = label)
		dat$label <- factor(label, levels = c(setdiff(unique(label), pos), pos))
		#--------------------------------------------------
		# Point [1.1] Bootstrap calculation of lower and upper thresholds
		#            directly written with reference to https://github.com/noellewarmenhoven/Lumipulse-Nat-Med/blob/main/Calculating%20thresholds.R
		#--------------------------------------------------
		boot_res <- replicate(n.boot, {
					
					idx <- sample(seq_len(nrow(dat)), replace = TRUE)
					d <- dat[idx, ]
					roc.obj <- roc(d$label, d$score, quiet = TRUE, direction = "<")
					coords <- coords(
							roc.obj,
							x = "all",
							ret = c("threshold", "sensitivity", "specificity"),
							transpose = FALSE
					)
					
					# #Sensitivity at 95
					epsilon  <- 0.005   # for example ±0.005, making the coverage a small interval around the sensitivity neighborhood
					find_threshold_sens95 <- coords[coords$sensitivity >= (sensitivity-epsilon) & coords$sensitivity <= (sensitivity+epsilon), ]
					threshold_sens95 <- if (nrow(find_threshold_sens95) > 0) {
								find_threshold_sens95[which.max(find_threshold_sens95$specificity), 1]
							} else {
								NA
							}
					threshold.down = threshold_sens95
					
					# #Specificity at 95%
					find_threshold_spec95 <- coords[coords$specificity >= (specificity-epsilon) & coords$specificity <= (specificity+epsilon), ]
					threshold_spec95 <- if (nrow(find_threshold_spec95) > 0) {
								find_threshold_spec95[which.max(find_threshold_spec95$sensitivity), 1]
							} else {
								NA
							}
					threshold.up = threshold_spec95
					
					c(threshold.down = threshold.down,
							threshold.up   = threshold.up)
				})
		# Average Bootstrap results to obtain the final lower and upper thresholds
		boot_res <- t(boot_res)
		threshold.down <- mean(boot_res[, "threshold.down"], na.rm=T)
		threshold.up   <- mean(boot_res[, "threshold.up"], na.rm=T)
		
		#--------------------------------------------------
		# Point [2] Calculate performance using the lower and upper thresholds
		#--------------------------------------------------
		pred <- ifelse(
				score <= threshold.down, 0,
				ifelse(score >= threshold.up, 1, NA)
		)
		
		Intermediate.prop = sum(is.na(pred))/length(pred)
		keep <- !is.na(pred)
		
		TP <- sum(pred[keep] == 1 & label[keep] == "1")
		TN <- sum(pred[keep] == 0 & label[keep] == "0")
		FP <- sum(pred[keep] == 1 & label[keep] == "0")
		FN <- sum(pred[keep] == 0 & label[keep] == "1")
		
		PPV <- TP / (TP + FP)
		NPV <- TN / (TN + FN)
		ACC <- (TP + TN) / sum(keep)
		
	} else {
		threshold.down = threshold.up = PPV = NPV = ACC = AUROC.full = Intermediate.prop = NA
	}
	
	threshold_down_and_up <- c(
			threshold.down = threshold.down,
			threshold.up   = threshold.up,
			PPV = PPV,
			NPV = NPV,
			ACC = ACC,
			AUROC.full = AUROC.full,
			Intermediate.porp = Intermediate.prop
	)
	
	
	
	#'-Part [2]-
	#' Content: visualize density distributions
	#'#######
	df_val_sub = data.frame(pred.Prob=score, True.lab=label)
	density.fig <- function(df_val_sub, model_name, AUROC.full, Intermediate.prop){
		library(ggplot2)
		p <- ggplot(df_val_sub, aes(x = pred.Prob, color = True.lab, fill = True.lab)) +
				#1. Plot mutant and wild-type separately
				geom_density(alpha = 0.05) +
				scale_color_manual(values = c("0" = "blue", "1" = "red")) +
				scale_fill_manual(values = c("0" = "blue", "1" = "red")) +		
				labs(
						title = model_name,
						subtitle = paste0(
								"AUROC = ", sprintf("%.3f", AUROC.full),
								"\n Intermediate = ", sprintf("%.1f%%", Intermediate.prop * 100)
						),
						x = "Score",
						y = "Density"
				) +
				theme_classic(base_size = base_size)		
		return(p)
	}
	
	p = density.fig(df_val_sub=df_val_sub, model_name=main, AUROC.full, Intermediate.prop)
	if(need.log) p = p+scale_y_continuous(trans="log1p")  #
	
	## Lower bound: blue
	if (!is.null(threshold.down) && is.finite(threshold.down)){
		p <- p + geom_vline(xintercept = threshold.down,color = "blue",linewidth = 1.1,linetype = "dashed")
	}
	## Upper bound: red
	if (!is.null(threshold.up) && is.finite(threshold.up)){
		p <- p + geom_vline(xintercept = threshold.up,color = "red",linewidth = 1.1,linetype = "dashed")
	}
	
	
	return(list(
					threshold_down_and_up = threshold_down_and_up,
					p.plot = p
			))
}






