# TODO: Plot ROC and PR curves for classification performance evaluation, and return their areas under the curve.
# 
# Author: Administrator
###############################################################################

# Note: ROC sorts scores in ascending order and then calculates sensitivity and specificity across thresholds; therefore, higher-ranked scores correspond more strongly to class 1
#' plot ROC and calculate AUC. If there are multiple ROC curves (not cross-validation curves), they can be plotted one by one and then added together
#' @param score : a score list. If multiple curves are to be displayed, scores can be stored as a list, where each element corresponds to one curve
#' If isCrossValidation is FALSE, names must be provided for each element
#' @param label : a label list
########FUNCTION########
library(ggplot2)
library(pROC)
library(PRROC)
plot_performance_curves <- function(score.list, label.list, 
		curve.type = "ROC", 
		line.colors = NULL, 
		show.baseline = TRUE,
		legend.info = NULL,
		main = NULL,
		lwd = 1,
		size = 20,
		legend.position = "top" #none
){
	# 0) Set colors and model names
	n_models <- length(score.list)
	model.names = legend.info
	
	# Process model names
	if (is.null(model.names)) {
		model.names <- paste0("Model ", 1:n_models)
	} else if (length(model.names) < n_models) {
		warning("提供的模型名称数量少于模型数量，将使用默认名称补充")
		model.names <- c(model.names, paste0("Model ", (length(model.names)+1):n_models))
	}
	
	# Set colors
	if (is.null(line.colors)) {
		line.colors <- rainbow(n_models)
	} else if (length(line.colors) < n_models) {
		warning("提供的颜色数量少于模型数量，将循环使用颜色")
		line.colors <- rep(line.colors, length.out = n_models)
	}
	
	# 1) Prepare data
	plot_data <- data.frame()
	baseline_data <- data.frame()  # store baseline data
	
	for (i in 1:n_models) {
		scores <- score.list[[i]]
		labels <- label.list[[i]]
		
		if (curve.type == "ROC") {
			roc_obj <- roc(labels, scores)
			temp_data <- data.frame(
					FPR = 1 - roc_obj$specificities,
					TPR = roc_obj$sensitivities,
					Model = model.names[i],
					AUC = as.numeric(roc_obj$auc)
			)
		} else if (curve.type == "PR") {
			pr_obj <- pr.curve(scores.class0 = scores[labels == 1],
					scores.class1 = scores[labels == 0], 
					curve = TRUE)
			temp_data <- data.frame(
					Recall = pr_obj$curve[, 1],
					Precision = pr_obj$curve[, 2],
					Model = model.names[i],
					AUPRC = pr_obj$auc.integral
			)
			
			# Calculate the baseline for each model (positive class proportion)
			pos_prop <- mean(labels == 1)
			baseline_data <- rbind(baseline_data, 
					data.frame(
							Model = model.names[i],
							pos_prop = pos_prop,
							color = line.colors[i]
					))
		} else {
			stop("curve.type必须是'ROC'或'PR'")
		}
		plot_data <- rbind(plot_data, temp_data)
	}
	
	plot_data$Model = factor(plot_data$Model, levels=unique(plot_data$Model) )
	
	
	# 2) Plot
	if (curve.type == "ROC") {
		p <- ggplot(plot_data, aes(x = FPR, y = TPR, color = Model)) +
				geom_line(linewidth = lwd) +
				scale_color_manual(values = line.colors,
						labels = paste0(unique(plot_data$Model), 
								" ", 
								round(unique(plot_data$AUC), 3))) +
				labs(x = "1 - Specificity", 
						y = "Sensitivity",
						title = main) +
				coord_equal() +
				scale_x_continuous(expand = c(0, 0)) +  # remove padding on the X-axis
				scale_y_continuous(expand = c(0, 0)) +  # remove padding on the Y-axis and set range
				theme_classic()+
				theme(legend.position = legend.position,
						axis.line = element_blank(),
						axis.text.y = element_text(size = size-2, angle = 90, hjust = 0.4, vjust = 0.4),
						axis.text.x = element_text(size = size-2, hjust = 0.4, vjust = -0.4),
						axis.title = element_text(size = size+1),
						legend.text = element_text(size = size),
						legend.title = element_blank(), 
						plot.title = element_text(size = size+2),
						panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.8) # gray border
				)
		
		if (show.baseline) {
			p <- p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray")
		}
		
	} else if (curve.type == "PR") {
		# Calculate the y-axis range (reserve 5% margin)
		y_range <- range(plot_data$Precision, na.rm = TRUE)
		y_padding <- 0#diff(y_range) * 0.05
		y_limits <- c(max(0, y_range[1] - y_padding), 
				min(1, y_range[2] + y_padding))
		
		p <- ggplot(plot_data, aes(x = Recall, y = Precision, color = Model)) +
				geom_line(size = lwd) +
				scale_color_manual(values = line.colors,
						labels = paste0(unique(plot_data$Model), 
								" ", 
								round(unique(plot_data$AUPRC), 3) )) +
				labs(x = "Recall", y = "Precision",
						title = main) +
				scale_x_continuous(expand = c(0, 0)) +  # remove padding on the X-axis
				scale_y_continuous(expand = c(0, 0), limits = y_limits) +  # remove padding on the Y-axis and set range
				coord_cartesian(ylim = y_limits) +  # use coord_cartesian to avoid clipping data
				theme_classic()+
				theme(legend.position = legend.position,
						axis.line = element_blank(),
						axis.text.y = element_text(size = size-2, angle = 90, hjust = 0.4, vjust = 0.4),
						axis.text.x = element_text(size = size-2, hjust = 0.4, vjust = -0.4),
						axis.title = element_text(size = size+1),
						legend.text = element_text(size = size),
						legend.title = element_blank(), 
						plot.title = element_text(size = size+2),
						panel.border = element_rect(color = "gray60", fill = NA, linewidth = 0.8) # gray border
				)
		
		if (show.baseline && nrow(baseline_data) > 0) {
			# Add a baseline with the corresponding color for each model
			for (i in 1:nrow(baseline_data)) {
				p <- p + 
						geom_hline(yintercept = baseline_data$pos_prop[i],
								linetype = "dashed", 
								color = baseline_data$color[i],
								alpha = 0.7)+
						annotate("text", 
								x = max(plot_data$Recall),  # or use max(plot_data$Recall)
								y = baseline_data$pos_prop[i], 
								label = paste0("Baseline: ", round(baseline_data$pos_prop[i], 3),"(AUPRCratio:", round(unique(plot_data$AUPRC)[i]/baseline_data$pos_prop[i],3),")" ),
								hjust = 1.1, vjust = -0.3,   # adjust position
								size = 6,
								color = baseline_data$color[i])
			}
		}
	}
	
	return(p)
}
