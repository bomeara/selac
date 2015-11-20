
######################################################################################################################################
######################################################################################################################################
### Sella and Hirsh model for nucleotides
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu

#source("selac.R")

CreateNucleotideDistanceMatrix <- function(position.multiplier) {
    n.states <- 4
    nucleotide.distances <- matrix(1,nrow=n.states,ncol=n.states)
    diag(nucleotide.distances) <- 0
    nucleotide.distances <- nucleotide.distances * position.multiplier
    rownames(nucleotide.distances) <- colnames(nucleotide.distances) <- n2s(0:3)
    return(nucleotide.distances)
}


GetNucleotideNucleotideDistance <- function(n1, n2, nucleotide.distances){
    site_d <- function(k){
        return(nucleotide.distances[n1[k], n2[k]])
    }
    d <- sapply(c(1:length(n1)),site_d,simplify=TRUE)
    return(d)
}


GetPairwiseNucleotideWeightSingleSite <- function(d1, d2, Ne, diploid){
    if(diploid==TRUE){
        b = 1
    }else{
        b = 2
    }
    if(d1==d2){ #When the fitnesses are the same, neutral case, pure drift
        return(1/(2*Ne))
    }else{
        fit_ratio <- exp(-(d1-d2)) #f1/f2
        if(fit_ratio==Inf) #1 is much better than 2 (the mutant)
        return(0)
        else if(fit_ratio==1)
        return(1/(2*Ne))
        else
        return((1-fit_ratio^b)/(1-fit_ratio^(2*Ne)))
    }
}


GetNucleotideFixationMatrix <- function(site.number, position.multiplier, optimal.nucleotide, Ne, diploid=TRUE){
    nucleotide.set <- 0:3
    nucleotide.distances <- CreateNucleotideDistanceMatrix(position.multiplier[site.number])
    nucleotide.fitness.ratios <- matrix(data=0,4,4)
    unique.nucs <- n2s(0:3)
    for (i in sequence(4)) {
        for (j in sequence(4)) {
            nuc1 <- n2s(nucleotide.set[i])
            nuc2 <- n2s(nucleotide.set[j])
            if(!nuc1 == nuc2){
                d1 <- GetProteinProteinDistance(protein1=nuc1, protein2=unique.nucs[optimal.nucleotide], aa.distances=nucleotide.distances)
                d2 <- GetProteinProteinDistance(protein1=nuc2, protein2=unique.nucs[optimal.nucleotide], aa.distances=nucleotide.distances)
                nucleotide.fitness.ratios[i,j] <- GetPairwiseNucleotideWeightSingleSite(d1=d1, d2=d2, Ne=Ne, diploid=diploid)
            }
        }
    }
    return(nucleotide.fitness.ratios)
}


GetLikelihoodUCEForSingleCharGivenOptimum <- function(charnum=1, nuc.data, phy, Q_position, root.p=NULL, scale.factor, return.all=FALSE) {
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    nl <- nrow(Q_position)
    #Now we need to build the matrix of likelihoods to pass to dev.raydisc:
    liks <- matrix(0, nb.tip + nb.node, nl)
    #Now loop through the tips.
    for(i in 1:nb.tip){
        #The codon at a site for a species is not NA, then just put a 1 in the appropriate column.
        #Note: We add charnum+1, because the first column in the data is the species labels:
        if(nuc.data[i,charnum+1] < 65){
            liks[i,nuc.data[i,charnum+1]] <- 1
        }else{
            #If here, then the site has no data, so we treat it as ambiguous for all possible codons. Likely things might be more complicated, but this can be modified later:
            liks[i,] <- 1
        }
    }
    #The result here is just the likelihood:
    result <- -GetLikelihood(phy=phy, liks=liks, Q=Q_position, scale.factor=scale.factor, root.p=root.p)
    ifelse(return.all, stop("return all not currently implemented"), return(result))
}


GetLikelihoodUCEForManyCharVaryingBySite <- function(nuc.data, phy, nuc.mutation.rates, position.multiplier.vector, Ne, nuc.optim_array=NULL, root.p_array=NULL, diploid=TRUE) {
    nsites <- dim(nuc.data)[2]-1
    final.likelihood.vector <- rep(NA, nsites)
    if(is.null(root.p_array)) {
        #Generate matrix of equal frequencies for each site:
        root.p_array <- rep(0.25, 4)
    }
    if(diploid == TRUE){
        ploidy = 2
    }else{
        ploidy = 1
    }
    for(site.index in sequence(nsites)) {
        weight.matrix <- GetNucleotideFixationMatrix(1, position.multiplier=position.multiplier.vector[site.index], optimal.nucleotide=nuc.optim_array[site.index], Ne=Ne, diploid=diploid)
        Q_position <- (ploidy * Ne) * nuc.mutation.rates * weight.matrix
        #Rescaling Q matrix in order to have a 1 nucleotide change per site if the branch length was 1:
        diag(Q_position) = 0
        Q_position <- t(Q_position * root.p_array)
        diag(Q_position) <- -rowSums(Q_position)
        scale.factor <- -sum(diag(Q_position) * root.p_array)
        final.likelihood.vector[site.index] <- GetLikelihoodUCEForSingleCharGivenOptimum(charnum=site.index, nuc.data=nuc.data, phy=phy, Q_position=Q_position, root.p=root.p_array, scale.factor=scale.factor, return.all=FALSE)
    }
    return(final.likelihood.vector)
}


PositionSensitivityMultiplier <- function(a0, a1, a2, site.index){
    #At the moment assumes a standard normal distribution
    sensitivity.vector <- a0 * exp(-((site.index - a1)^2)/ (2*(a2^2)))
    return(sensitivity.vector)
}


GetLikelihoodUCEForManyCharGivenAllParams <- function(x, nuc.data, phy, nuc.optim_array=NULL, nuc.model, diploid=TRUE, logspace=FALSE, verbose=TRUE, neglnl=FALSE) {
    if(logspace) {
        x = exp(x)
    }
    
    if(nuc.model == "JC") {
        base.freqs=c(x[5:7], 1-sum(x[5:7]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model)
    }
    if(nuc.model == "GTR") {
        base.freqs=c(x[5:7], 1-sum(x[5:7]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[8:length(x)], model=nuc.model)
    }
    if(nuc.model == "UNREST") {
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[8:length(x)], model=nuc.model)
    }
    nsites <- dim(nuc.data)[2]-1
    site.mid.point <- ceiling(nsites/2)
    site.index <- 1:nsites
    site.index <- site.index - site.mid.point
    #Note that I am rescaling x2. The reason being that I can optimize in log space and not have to deal with negatives.
    position.multiplier.vector <- PositionSensitivityMultiplier(x[1], x[2]-site.mid.point, x[3], site.index)
    final.likelihood = GetLikelihoodUCEForManyCharVaryingBySite(nuc.data=nuc.data, phy=phy, nuc.mutation.rates=nuc.mutation.rates, position.multiplier.vector=position.multiplier.vector, Ne=x[4], nuc.optim_array=nuc.optim_array, root.p_array=base.freqs, diploid=diploid)
    likelihood <- sum(final.likelihood)
    
    if(neglnl) {
        likelihood <- -1 * likelihood
    }
    if(is.na(likelihood)){
        return(10000000000)
    }
    if(verbose) {
        results.vector <- c(likelihood)
        names(results.vector) <- c("likelihood")
        print(results.vector)
    }
    return(likelihood)
}


OptimizeEdgeLengthsGlobalUCE <- function(x, site.pattern.data.list, n.partitions, nsites.vector, index.matrix, phy, nuc.optim.list=NULL, diploid=TRUE, nuc.model, logspace=FALSE, verbose=TRUE, n.cores=NULL, neglnl=FALSE) {
    if(logspace) {
        x <- exp(x)
    }
    par.mat <- index.matrix
    par.mat[] <- c(x, 0)[index.matrix]
    #sums the total number of parameters: 4 is the general shape pars, 3 are the base pars, and finally, the transition rates.
    if(nuc.model == "JC"){
        max.par = 4 + 3 + 0
    }
    if(nuc.model == "GTR"){
        max.par = 4 + 3 + 5
    }
    if(nuc.model == "UNREST"){
        max.par = 4 + 3 + 11
    }
    if(is.null(n.cores)){
        likelihood.vector <- c()
        for(partition.index in sequence(n.partitions)){
            phy$edge.length = par.mat[partition.index,(max.par+1):ncol(par.mat)]
            nuc.data = NULL
            nuc.data = site.pattern.data.list[[partition.index]]
            likelihood.vector = c(likelihood.vector, GetLikelihoodUCEForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, nuc.optim_array=nuc.optim.list[[partition.index]], nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl))
        }
        likelihood = sum(likelihood.vector)
    }else{
        MultiCoreLikelihood <- function(partition.index){
            phy$edge.length = par.mat[partition.index,(max.par+1):ncol(par.mat)]
            nuc.data = NULL
            nuc.data = site.pattern.data.list[[partition.index]]
            likelihood.tmp = GetLikelihoodUCEForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, nuc.optim_array=nuc.optim.list[[partition.index]], nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl)
            return(likelihood.tmp)
        }
        #This orders the nsites per partition in decreasing order (to increase efficiency):
        partition.order <- 1:n.partitions
        likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores)))
    }
    print(likelihood)
    return(likelihood)
}


######################################################################################################################################
######################################################################################################################################
### Likelihood calculator -- One step process because at each site Q changes
######################################################################################################################################
######################################################################################################################################

GetLikelihood <- function(phy, liks, Q, scale.factor, root.p){
    Q.scaled = Q * (1/scale.factor)
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    TIPS <- 1:nb.tip
    comp <- numeric(nb.tip + nb.node)
    phy <- reorder(phy, "pruningwise")
    # Obtain an object of all the unique ancestors
    anc <- unique(phy$edge[,1])
    for (i  in seq(from = 1, length.out = nb.node)) {
        # The ancestral node at row i is called focal
        focal <- anc[i]
        # Get descendant information of focal
        desRows<-which(phy$edge[,1]==focal)
        desNodes<-phy$edge[desRows,2]
        v <- 1
        for (desIndex in sequence(length(desRows))){
            v <- v*expm(Q.scaled * phy$edge.length[desRows[desIndex]]) %*% liks[desNodes[desIndex],]
        }
        comp[focal] <- sum(v)
        liks[focal, ] <- v/comp[focal]
    }
    # Specifies the root:
    root <- nb.tip + 1L
    # If any of the logs have NAs restart search:
    if(is.nan(sum(log(comp[-TIPS]))) || is.na(sum(log(comp[-TIPS])))){
        return(10000000000)
    }
    else{
        loglik<- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
        if(is.infinite(loglik)){
            return(10000000000)
        }
    }
    loglik
}

GetMaxNameUCE <- function(x) {
    x = unname(x)
    x = x[which(x!=65)]
    return(names(table(x))[(which.is.max(table(x)))]) #note that this breaks ties at random
}


######################################################################################################################################
######################################################################################################################################
### Organizes the optimization routine
######################################################################################################################################
######################################################################################################################################

#' @title Optimize parameters under the SELON model
#'
#' @description 
#' Optimizes model parameters under the SELON model
#'
#' @param nuc.data.path Provides the path to the directory containing the gene specific fasta files that contains the nucleotide data.
#' @param n.partitions The number of partitions to analyze. The order is based on the Unix order of the fasta files in the directory.
#' @param phy The phylogenetic tree to optimize the model parameters.
#' @param edge.length A logical indicating whether or not edge lengths should be optimized.
#' @param edge.linked A logical indicating whether or not edge lengths should be optimized separately for each gene. By default, a single set of each lengths is optimized for all genes.
#' @param optimal.nuc Indicates what type of optimal.nuc should be used. At the moment there is only a single option: "majrule".
#' @param nuc.model Indicates what type nucleotide model to use. There are three options: "JC", "GTR", or "UNREST".
#' @param include.gamma A logical indicating whether or not to include a discrete gamma model.
#' @param diploid A logical indicating whether or not the organism is diploid or not.
#' @param verbose Logical indicating whether each iteration be printed to the screen.
#' @param n.cores The number of cores to run the analyses over.
#' @param max.tol Supplies the relative optimization tolerance.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.
#'
#' @details 
#' SELON stands for SELection On Nucleotides. This function takes a user supplied topology and a set of fasta formatted sequences and optimizes the parameters in the SELON model. Selection is based on selection towards an optimal nucleotide at each site, which is based simply on the majority rule of the observed data. The strength of selection is then varied along sites based on a Taylor series, which scales the substitution rates. NOTE THIS IS NOT WORKING PROPERLY AT THIS TIME. SO PLEASE DO NOT USE YET.
SelonOptimize <- function(nuc.data.path, n.partitions=NULL, phy, edge.length="optimize", edge.linked=TRUE, optimal.nuc="majrule", nuc.model="GTR", diploid=TRUE, verbose=FALSE, n.cores=NULL, max.tol=.Machine$double.eps^0.25, fasta.rows.to.keep=NULL) {
    
    cat("Initializing data and model parameters...", "\n")
    
    partitions <- system(paste("ls -1 ", nuc.data.path, "*.fasta", sep=""), intern=TRUE)
    
    if(is.null(n.partitions)){
        n.partitions <- length(partitions)
    }else{
        n.partitions = n.partitions
    }
    site.pattern.data.list <- as.list(numeric(n.partitions))
    nuc.optim.list <- as.list(numeric(n.partitions))
    nsites.vector <- c()
    empirical.base.freq.list <- as.list(numeric(n.partitions))
    starting.branch.lengths <- matrix(0, n.partitions, length(phy$edge[,1]))
    for (partition.index in sequence(n.partitions)) {
        gene.tmp <- read.dna(partitions[partition.index], format='fasta')
        if(!is.null(fasta.rows.to.keep)){
            gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
        }else{
            gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
        }
        starting.branch.lengths[partition.index,] <- ComputeStartingBranchLengths(phy, gene.tmp)$edge.length
        nucleotide.data <- DNAbinToNucleotideNumeric(gene.tmp)
        nucleotide.data <- nucleotide.data[phy$tip.label,]
        site.pattern.data.list[[partition.index]] = nucleotide.data
        nsites.vector = c(nsites.vector, dim(nucleotide.data)[2] - 1)
        empirical.base.freq <- as.matrix(nucleotide.data[,-1])
        empirical.base.freq <- table(empirical.base.freq, deparse.level = 0) / sum(table(empirical.base.freq, deparse.level = 0))
        empirical.base.freq.list[[partition.index]] <- as.vector(empirical.base.freq[1:4])
        nuc.optim <- as.numeric(apply(nucleotide.data[,-1], 2, GetMaxNameUCE)) #starting values for all, final values for majrule
        nuc.optim.list[[partition.index]] <- nuc.optim
    }
    opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = "100000", "ftol_rel" = max.tol)
    results.final <- c()
    if(nuc.model == "JC"){
        nuc.ip = NULL
        ip = c(1e-7, ceiling(nsites.vector[1]/2), exp(-21), 5e6, 0.25, 0.25, 0.25)
        upper = c(21, nsites.vector[1], 21, 21, 0, 0, 0)
        lower = rep(-21, length(ip))
        max.par.model.count = 4 + 3 + 0
    }
    if(nuc.model == "GTR"){
        nuc.ip = rep(1, 5)
        ip = c(1e-7, ceiling(nsites.vector[1]/2), exp(-21), 5e6, 0.25, 0.25, 0.25, nuc.ip)
        upper = c(21, nsites.vector[1], 21, 21, 0, 0, 0, rep(21, length(nuc.ip)))
        lower = rep(-21, length(ip))
        max.par.model.count = 4 + 3 + 5
    }
    if(nuc.model == "UNREST"){
        nuc.ip = rep(1, 11)
        ip = c(1e-7, ceiling(nsites.vector[1]/2), exp(-21), 5e6, 0.25, 0.25, 0.25, nuc.ip)
        upper = c(21, nsites.vector[1], 21, 21, 0, 0, 0, rep(21, length(nuc.ip)))
        lower = rep(-21, length(ip))
        max.par.model.count = 4 + 3 + 11
    }
    if(edge.length == "optimize"){
        if(edge.linked == TRUE){
            phy$edge.length <- colMeans(starting.branch.lengths)
            index.matrix = matrix(0, n.partitions, max.par.model.count+length(phy$edge.length))
            index.matrix[1,] = 1:ncol(index.matrix)
            ip.vector = c(ip, phy$edge.length)
            upper.vector = c(upper, rep(log(5), length(phy$edge.length)))
            lower.vector = c(lower, rep(-21, length(phy$edge.length)))
            for(partition.index in 2:n.partitions){
                if(nuc.model == "JC"){
                    ip[2] = ceiling(nsites.vector[partition.index]/2)
                    upper[2] = nsites.vector[partition.index]
                    ip.vector = c(ip.vector, ip[1:3], ip[5], ip[6], ip[7])
                    upper.vector = c(upper.vector, c(upper[1:3], upper[5], upper[6], upper[7]))
                    lower.vector = c(lower.vector, c(lower[1:3], lower[5], lower[6], lower[7]))
                }else{
                    ip[2] = ceiling(nsites.vector[partition.index]/2)
                    upper[2] = nsites.vector[partition.index]
                    ip.vector = c(ip.vector, ip[1:3], ip[5], ip[6], ip[7], nuc.ip)
                    upper.vector = c(upper.vector, c(upper[1:3], upper[5], upper[6], upper[7], rep(21, length(nuc.ip))))
                    lower.vector = c(lower.vector, c(lower[1:3], lower[5], lower[6], lower[7], rep(-21, length(nuc.ip))))
                }
                index.matrix.tmp = numeric(max.par.model.count + length(phy$edge.length))
                #This fixes Ne, and the edge lengths across all partitions:
                index.matrix.tmp[4] = 4
                index.matrix.tmp[(max.par.model.count+1):ncol(index.matrix)] = index.matrix[1,(max.par.model.count+1):ncol(index.matrix)]
                index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
                index.matrix[partition.index,] <- index.matrix.tmp
            }
        }else{
            phy$edge.length <- starting.branch.lengths[1,]
            index.matrix = matrix(0, n.partitions, max.par.model.count+length(phy$edge.length))
            index.matrix[1,] = 1:ncol(index.matrix)
            ip.vector = c(ip, phy$edge.length)
            upper.vector = c(upper, rep(log(5), length(phy$edge.length)))
            lower.vector = c(lower, rep(-21, length(phy$edge.length)))
            for(partition.index in 2:n.partitions){
                phy$edge.length <- starting.branch.lengths[partition.index,]
                if(nuc.model == "JC"){
                    ip[2] = ceiling(nsites.vector[partition.index]/2)
                    upper[2] = nsites.vector[partition.index]
                    ip.vector = c(ip.vector, ip[1:3], ip[5], ip[6], ip[7], phy$edge.length)
                    upper.vector = c(upper.vector, c(upper[1:3], upper[5], upper[6], upper[7], rep(log(5), length(phy$edge.length))))
                    lower.vector = c(lower.vector, c(lower[1:3], lower[5], lower[6], lower[7], rep(log(-21), length(phy$edge.length))))
                }else{
                    ip[2] = ceiling(nsites.vector[partition.index]/2)
                    upper[2] = nsites.vector[partition.index]
                    ip.vector = c(ip.vector, ip[1:3], ip[5], ip[6], ip[7], nuc.ip, phy$edge.length)
                    upper.vector = c(upper.vector, c(upper[1:3], upper[5], upper[6], upper[7], rep(21, length(nuc.ip)), rep(log(5), length(phy$edge.length))))
                    lower.vector = c(lower.vector, c(lower[1:3], lower[5], lower[6], lower[7], rep(-21, length(nuc.ip)), rep(-21, length(phy$edge.length))))
                }
                index.matrix.tmp = numeric(max.par.model.count+length(phy$edge.length))
                #This fixes Ne across all partitions:
                index.matrix.tmp[4] = 4
                index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
                index.matrix[partition.index,] <- index.matrix.tmp							
            }
        }
        cat("Finished. Optimizing model parameters...", "\n")
        results.final <- nloptr(x0=log(ip.vector), eval_f = OptimizeEdgeLengthsGlobalUCE, ub=upper.vector, lb=lower.vector, opts=opts, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=nuc.optim.list, diploid=diploid, nuc.model=nuc.model, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
        cat("Finished. Summarizing results...", "\n")
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
        phy$edge.length <- apply(mle.pars.mat[,7:ncol(mle.pars.mat)], 2, weighted.mean, w=nsites.vector)
        loglik <- -(results.final$objective) #to go from neglnl to lnl
        np = max(index.matrix)
        obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, AICc = NULL, mle.pars=mle.pars.mat, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, nuc.optim=nuc.optim_list, nuc.model=nuc.model, diploid=diploid, empirical.base.freqs=empirical.base.freq.list, max.tol=max.tol)
        class(obj) = "selac"
    }
    return(obj)
}


######################################################################################################################################
######################################################################################################################################
### TESTING CODE
######################################################################################################################################
######################################################################################################################################

#kk <- PositionSensitivityMultiplier(1e-7, 0,10, site.index<-seq(0, 100, by=1))
#ss <- matrix(1, 4,4)
#diag(ss) = 0
#GetNucleotideFixationMatrix(1, position.multiplier=kk, 1, 5e6)*(2*5e6)*matrix(1,4,4)







