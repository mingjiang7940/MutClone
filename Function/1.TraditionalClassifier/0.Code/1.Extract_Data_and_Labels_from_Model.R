# TODO: Add comment
# 
# Author: Administrator
###############################################################################

# Extract feature data and label data from the model
ExtractDataLabel.ClassifyModel <- function(cm, way="mut-exprs", background=FALSE)
{
	# Data is directly stored in the model
	if(!is.null(cm$data) & !is.null(cm$sample.label)){
		return(list(data=cm$data, label=cm$sample.label))
	}
	
	
	if(way == "mut-exprs")#use standard data storage format
	{
		t.data = readRDS(cm$data.path.rds)
		all.data = t.data$all.data
		all.label = t.data$all.label
		
		temp = all.label[[cm$data.name]]
		
		if(background){
			data = all.data[temp$background, names(temp$label)]
		}else{
			data = all.data[ ,names(temp$label)]
		}
		
		
		
		return(list(data=data, label=temp$label))
	}
	
	
	stop("Model data not found. Check data path or global variable.")
	
}
