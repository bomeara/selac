
######################################################################################################################################
######################################################################################################################################
### Adaptive Bootstrap -- Simulating confidence intervals for parameters estimated in selac
######################################################################################################################################
######################################################################################################################################
source("selac.R")

SupportRegion <- function(selac.obj, n.points=10000, scale.int=0.1, desired.delta=2, nuc.model="JC", edge.length="optimize", verbose=TRUE){
	phy <- selac.obj$phy
	data <- selac.obj$codon.data
	data.new<-data[phy$tip.label,]
	if(edge.length == "optimize"){
		par <- c(selac.obj$mle.pars, selac.obj$phy$edge.length[-1])
	}else{
		par <- c(selac.obj$mle.pars)
	}
	lower <- rep(exp(-20), length(par))
	upper <- rep(exp(20), length(par))
	
	interval.results <- AdaptiveConfidenceIntervalSampling(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data.new, nuc.model=nuc.model, scale.int=scale.int)
	interval.results.in <- interval.results[which(interval.results[,1] - min(interval.results[,1])<=2),]
	interval.results.out <- interval.results[which(interval.results[,1] - min(interval.results[,1])>2),]
	
	obj = NULL	
	obj$ci = apply(interval.results.in ,2, quantile)
	obj$points.within.region = interval.results.in
	obj$all.points = interval.results
	class(obj) = "selac.support"	
	return(obj)
}


AdaptiveConfidenceIntervalSampling <- function(par, lower, upper, desired.delta=2, n.points=5000, verbose=TRUE, phy, data, nuc.model, scale.int) {
	
	phy$node.label <- NULL
	codon.data <- SitePattern(data)
	aa.data <- ConvertCodonNumericDataToAAData(codon.data$unique.site.patterns, numcode=1)
	aa.optim <- apply(aa.data[, -1], 2, GetMaxName)
	aa.properties=NULL
	#############################################################
	model.pars.index = 5
	if(nuc.model == "JC"){
		model.pars.index = model.pars.index+1
	}
	#Now assess the likelihood at the MLE:
	starting <- -GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(par[1:model.pars.index]), codon.data, phy, aa.optim_array=aa.optim, root.p_array=NULL, numcode=1, aa.properties=NULL, nuc.model=nuc.model, logspace=TRUE)
	#Generate the multipliers for feeling the boundaries:
	min.multipliers <- rep(1, length(par))
	max.multipliers <- rep(1, length(par))
	results <- data.frame(data.frame(matrix(nrow=n.points+1, ncol=1+length(par))))
	results[1,] <- unname(c(starting, par))
	write.table(results[1,], file="supportPointsSet", quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)
	for (i in sequence(n.points)) {
		sim.points <- NA
		while(is.na(sim.points[1])) {
			sim.points <- GenerateValues(par, lower=lower, upper=upper, scale.int=scale.int, examined.max=max.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, max, na.rm=TRUE), examined.min=min.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, min, na.rm=TRUE))
		}
		phy$edge.length = c(1, sim.points[(model.pars.index+1):length(par)])
		second <- -GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(sim.points[1:model.pars.index]), codon.data, phy, aa.optim_array=aa.optim, root.p_array=NULL, numcode=1, aa.properties=NULL, nuc.model=nuc.model, logspace=TRUE)
		results[i+1,] <- c(second, sim.points)
		write.table(results[i+1,], file="supportPointsSet", quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE, append=TRUE)
		if(i%%20==0) {
			for (j in sequence(length(par))) {
				returned.range <- range(results[which((results[,1]-min(results[,1], na.rm=TRUE))<desired.delta), j+1], na.rm=TRUE)
				total.range <- range(results[,j+1], na.rm=TRUE)
				width.ratio <- diff(returned.range)/diff(total.range)
				if(is.na(width.ratio)) {
					width.ratio=1	
				}
				if(width.ratio > 0.5) { #we are not sampling widely enough
					min.multipliers[j] <- min.multipliers[j] * (1-scale.int)
					max.multipliers[j] <- max.multipliers[j] * (1+scale.int) #expand the range
				} else {
					min.multipliers[j] <- 1
					max.multipliers[j] <- 1
				}
			}
		}
		if (verbose && i%%100==0) {
			print(paste(i, "of", n.points, "points done"))	
		}
	}
	return(results)
}


GenerateValues <- function(par, lower, upper, scale.int, max.tries=100, expand.prob=0, examined.max, examined.min) {
	pass=FALSE
	tries=0
	while(!pass && tries<=max.tries) {
		tries <- tries+1
		pass=TRUE
		new.vals <- rep(NA, length(par))
		for(i in sequence(length(par))) {
			examined.max[i] <- max(0.001, examined.max[i])
			new.vals[i] <- runif(1, max(lower[i], (1-scale.int)*examined.min[i]), min(upper[i], (1+scale.int)*examined.max[i]))
			if(new.vals[i]<lower[i]) {
				pass=FALSE
			}
			if(new.vals[i]>upper[i]) {
				pass=FALSE
			}
		}
	}
	if(tries>max.tries) {
		return(NA)
	}
	return(new.vals)
}


print.selac.support <- function(x,...){

	cat("\nSupport Region\n")
	print(x$ci[,-1])
	cat("\n")
}
	

