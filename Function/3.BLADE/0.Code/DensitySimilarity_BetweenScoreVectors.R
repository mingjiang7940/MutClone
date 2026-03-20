#'##########  ---- Paragraph <1> ----  ###########  
#' Core objective: Jensen-Shannon divergence
#' 
#'############ [25-11-18] ############
# Calculate the Jensen-Shannon similarity between two distributions based on two probability vectors
Get_JensenShannon_similarity <- function(score.ls1, score.ls2){
	
	# Estimate kernel density and convert to discrete probability distributions
	x_range <- range(c(score.ls1, score.ls2))
	dens_1 <- density(score.ls1, from=x_range[1], to=x_range[2], n=512)
	dens_2  <- density(score.ls2,  from=x_range[1], to=x_range[2], n=512)
	
	sim <- JS_similarity(dens_1$y, dens_2$y)
	return(sim)
}


#Point【Sub-function】: compute JS similarity from density distributions
JS_similarity <- function(p, q){
	library(philentropy)
	# Normalize probabilities
	p <- p / sum(p)
	q <- q / sum(q)
	# Use philentropy::JSD
	Djs <- JSD(rbind(p,q), unit="log2")
	sim <- 1 - Djs   # Convert to similarity in the range 0~1
	return(sim)
}