
######################################################################################################################################
######################################################################################################################################
### Adaptive Bootstrap -- Simulating confidence intervals for parameters estimated in selac
######################################################################################################################################
######################################################################################################################################
#source("selac.R")

#' @title Adaptive sampling of support region
#'
#' @description
#' Adaptively samples points for each parameter to obtain an of the confidence intervals
#'
#' @param selac.obj An object of class "selac" that contains the MLE from a given model run.
#' @param n.points Indicates the number of points to sample.
#' @param scale.int The scaling multiplier that defines the interval to randomly sample. By default it is set to 0.1, meaning that values are drawn at random along an interval that encompasses 10 percent above and below the MLE.
#' @param desired.delta Defines the number of lnL units away from the MLE to include. By default the value is set to 2.
#' @param parallel.type Designates whether a parallel run should occur by gene ("by.gene") or by site ("by.site").
#' @param n.cores The number of cores to run analyses over.
#' @param verbose A logical indicating whether progress should be printed to the screen. The default is TRUE.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.
#'
#' @details
#' This function provides a means for sampling the likelihood surface quickly to estimate confidence intervals that reflect the uncertainty in the MLE. The function starts with the MLE from the hisse run. It then uses a scaling multiplier to generate an interval by which to randomly alter each parameter. However, the algorithm was designed to \dQuote{feel} the boundaries of the random search. In other words, when the algorithm begins to sample the hinterlands of the surface, it will know to restrict the boundary to allow sampling of more reasonable values based on the currently sampled set. The goal of this sampling process is to find points within some desired distance from the MLE; by default we assume this distance is 2 lnLik. The confidence interval can be estimated directly from these points. The full set of points tried are also provided and can be used to generate contour plots (though, it is not entirely straightforward to do so -- but certainly doable).
SupportRegion <- function(selac.obj, n.points=10000, scale.int=0.1, desired.delta=2, parallel.type="by.gene", n.cores=NULL, verbose=TRUE, fasta.rows.to.keep=NULL) {
	
    phy <- selac.obj$phy
    partitions <- selac.obj$partitions
    n.partitions <- length(partitions)
    pars.mat <- selac.obj$mle.pars
    pars.index <- selac.obj$index.matrix
    pars <- c()
    #Ok, this is all going to be rather clunky. So bear with me...
    base.freqs.names <- c()
    if(!selac.obj$nuc.model == "JC"){
        rate.names <- c()
    }

    Ne.name <- "Ne"
    alpha.name <- "alpha"
    beta.name <- "beta"
    edge.length.names <- rep("edge.length", length(phy$edge.length))
    par.names <- c()
    
    for(partition.index in 1:n.partitions){
        if(selac.obj$nuc.model == "JC"){
            c.phi.q.name <- paste("C.phi.q", partition.index, sep="_")
            base.freqs.names <- paste(c("freqA", "freqC", "freqG"), partition.index, sep="_")
            tmp.names <- c(c.phi.q.name, alpha.name, beta.name, Ne.name, base.freqs.names, edge.length.names)
            par.names <- c(par.names, tmp.names[which(pars.index[partition.index,]>=pars.index[partition.index,1])])
            pars <- c(pars, pars.mat[partition.index,which(pars.index[partition.index,]>=pars.index[partition.index,1])])
        }
        if(selac.obj$nuc.model == "GTR"){
            c.phi.q.name <- paste("C.phi.q", partition.index, sep="_")
            base.freqs.names <- paste(c("freqA", "freqC", "freqG"), partition.index, sep="_")
            rate.names <- paste(c("C_A", "G_A", "T_A", "G_C", "T_C"), partition.index, sep="_")
            tmp.names <- c(c.phi.q.name, alpha.name, beta.name, Ne.name, base.freqs.names, rate.names, edge.length.names)
            par.names <- c(par.names, tmp.names[which(pars.index[partition.index,]>=pars.index[partition.index,1])])
            pars <- c(pars, pars.mat[partition.index,which(pars.index[partition.index,]>=pars.index[partition.index,1])])
        }
        if(selac.obj$nuc.model == "UNREST"){
            c.phi.q.name <- paste("C.phi.q", partition.index, sep="_")
            base.freqs.names <- paste(c("freqA", "freqC", "freqG"), partition.index, sep="_")
            rate.names <- paste(c("C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T"), partition.index, sep="_")
            tmp.names <- c(c.phi.q.name, alpha.name, beta.name, Ne.name, base.freqs.names, rate.names, edge.length.names)
            par.names <- c(par.names, tmp.names[which(pars.index[partition.index,]>=pars.index[partition.index,1])])
            pars <- c(pars, pars.mat[partition.index,which(pars.index[partition.index,]>=pars.index[partition.index,1])])
        }
    }
 
    site.pattern.data.list <- as.list(numeric(n.partitions))
    site.pattern.count.list <- as.list(numeric(n.partitions))
    nsites.vector <- c()
    if(selac.obj$aa.optim.type == "none"){
        empirical.base.freq.list <- as.list(numeric(n.partitions))
        for (partition.index in sequence(n.partitions)) {
            gene.tmp <- read.dna(partitions[partition.index], format='fasta')
            if(!is.null(fasta.rows.to.keep)){
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
            }else{
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
            }
            nucleotide.data <- DNAbinToNucleotideNumeric(gene.tmp)
            nucleotide.data <- nucleotide.data[phy$tip.label,]
            nsites.vector = c(nsites.vector, dim(nucleotide.data)[2] - 1)
            empirical.base.freq <- as.matrix(nucleotide.data[,-1])
            empirical.base.freq <- table(empirical.base.freq, deparse.level = 0)/sum(table(empirical.base.freq, deparse.level = 0))
            empirical.base.freq.list[[partition.index]] <- as.vector(empirical.base.freq[1:4])
            nucleotide.data <- SitePattern(nucleotide.data, includes.optimal.aa=FALSE)
            site.pattern.data.list[[partition.index]] = nucleotide.data$unique.site.patterns
            site.pattern.count.list[[partition.index]] = nucleotide.data$site.pattern.counts
        }
    }else{
        empirical.codon.freq.list <- as.list(numeric(n.partitions))
        aa.optim.list <- as.list(numeric(n.partitions))
        aa.optim.full.list <- as.list(numeric(n.partitions))
        for (partition.index in sequence(n.partitions)) {
            gene.tmp <- read.dna(partitions[partition.index], format='fasta')
            if(!is.null(fasta.rows.to.keep)){
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
            }else{
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
            }
            codon.data <- DNAbinToCodonNumeric(gene.tmp)
            codon.data <- codon.data[phy$tip.label,]
            nsites.vector = c(nsites.vector, dim(codon.data)[2] - 1)
            aa.optim <- selac.obj$aa.optim[[partition.index]]
            aa.optim.full.list[[partition.index]] <- aa.optim
            empirical.codon.freq.list[[partition.index]] <- CodonEquilibriumFrequencies(codon.data[,-1], aa.optim, numcode=selac.obj$numcode)
            aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
            colnames(aa.optim.frame.to.add) <- colnames(codon.data)
            codon.data <- rbind(codon.data, aa.optim.frame.to.add)
            codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
            site.pattern.data.list[[partition.index]] = codon.data$unique.site.patterns
            site.pattern.count.list[[partition.index]] = codon.data$site.pattern.counts
            aa.optim.list[[partition.index]] = codon.data$optimal.aa		
        }
    }
    codon.index.matrix = CreateCodonMutationMatrixIndex()
	lower <- rep(exp(-21), length(pars))
	upper <- rep(exp(21), length(pars))
    
    interval.results <- AdaptiveConfidenceIntervalSampling(x=pars, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=pars.index, phy=phy, aa.optim_array=aa.optim.list, root.p_array=empirical.codon.freq.list, numcode=selac.obj$numcode, diploid=selac.obj$diploid, aa.properties=selac.obj$aa.properties, volume.fixed.value=0.0003990333, nuc.model=selac.obj$nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=selac.obj$include.gamma, ncats=selac.obj$ncats, k.levels=selac.obj$k.levels, logspace=TRUE, verbose=verbose, parallel.type=parallel.type, n.cores=n.cores, neglnl=TRUE, lower=lower, upper=upper, desired.delta=desired.delta, n.points=n.points, scale.int=scale.int)
	
    par.names <- c("lnLik", par.names)
    interval.results.in <- interval.results[which(interval.results[,1] - min(interval.results[,1])<=2),]
	interval.results.out <- interval.results[which(interval.results[,1] - min(interval.results[,1])>2),]
	
    colnames(interval.results.in) <- colnames(interval.results) <- par.names
    
	obj = NULL	
	obj$ci = apply(interval.results.in, 2, quantile)
	obj$points.within.region = interval.results.in
	obj$all.points = interval.results
	class(obj) = "selac.support"	
	return(obj)
}


AdaptiveConfidenceIntervalSampling <- function(x, codon.site.data, codon.site.counts, n.partitions, nsites.vector, index.matrix, phy, aa.optim_array, root.p_array, numcode, diploid, aa.properties, volume.fixed.value, nuc.model, codon.index.matrix, include.gamma, ncats, k.levels, logspace, verbose, parallel.type, n.cores, neglnl, lower, upper, desired.delta, n.points, scale.int) {
	
	phy$node.label <- NULL
	
    #Now assess the likelihood at the MLE:
	starting <- OptimizeEdgeLengthsGlobal(x=log(x), codon.site.data=codon.site.data, codon.site.counts=codon.site.counts, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim_array, root.p_array=root.p_array, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, parallel.type=parallel.type, n.cores=n.cores, neglnl=neglnl)
    
    #Generate the multipliers for feeling the boundaries:
	min.multipliers <- rep(1, length(x))
	max.multipliers <- rep(1, length(x))
	
    results <- data.frame(data.frame(matrix(nrow=n.points+1, ncol=1+length(x))))
	results[1,] <- unname(c(starting, x))
	
    for (i in sequence(n.points)) {
        sum.vals <- NA
		while(is.na(sum.vals)) {
			sim.points <- GenerateValues(par=x, lower=lower, upper=upper, scale.int=scale.int, examined.max=max.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, max, na.rm=TRUE), examined.min=min.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, min, na.rm=TRUE))
            sum.vals <- sum(sim.points)
        }
        second <- OptimizeEdgeLengthsGlobal(x = log(sim.points), codon.site.data=codon.site.data, codon.site.counts=codon.site.counts, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim_array, root.p_array=root.p_array, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, parallel.type=parallel.type, n.cores=n.cores, neglnl=neglnl)
        results[i+1,] <- c(second, sim.points)
		
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
			examined.max[i] <- max(1e-10, examined.max[i])
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


