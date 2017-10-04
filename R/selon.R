
######################################################################################################################################
######################################################################################################################################
### Sella and Hirsh model for nucleotides
######################################################################################################################################
######################################################################################################################################


CreateNucleotideDistanceMatrix <- function() {
    n.states <- 4
    nucleotide.distances <- matrix(1,nrow=n.states,ncol=n.states)
    diag(nucleotide.distances) <- 0
    rownames(nucleotide.distances) <- colnames(nucleotide.distances) <- n2s(0:3)
    return(nucleotide.distances)
}


CreateNucleotideMutationMatrixSpecial <- function(rates) {
    index <- matrix(NA, 4, 4)
    np <- 12
    index[col(index) != row(index)] <- 1:np
    nuc.mutation.rates <- matrix(0, nrow=4, ncol=4)
    nuc.mutation.rates<-matrix(rates[index], dim(index))
    rownames(nuc.mutation.rates) <- n2s(0:3)
    colnames(nuc.mutation.rates) <- n2s(0:3)
    nuc.mutation.rates[3,4] = 1
    diag(nuc.mutation.rates) <- 0
    diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
    #Next we take our rates and find the homogeneous solution to Q*pi=0 to determine the base freqs:
    base.freqs <- Null(nuc.mutation.rates)
    #Rescale base.freqs so that they sum to 1:
    base.freqs.scaled <- c(base.freqs/sum(base.freqs))
    base.freqs.scaled.matrix <- rbind(base.freqs.scaled, base.freqs.scaled, base.freqs.scaled, base.freqs.scaled)
    diag(nuc.mutation.rates) <- 0
    #Rescale Q to account for base.freqs:
    nuc.mutation.rates <- nuc.mutation.rates * base.freqs.scaled.matrix
    diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
    obj <- NULL
    obj$base.freq <- base.freqs.scaled
    obj$nuc.mutation.rates <- nuc.mutation.rates
    return(obj)
}


GetNucleotideNucleotideDistance <- function(n1, n2, nucleotide.distances){
    site_d <- function(k){
        return(nucleotide.distances[n1[k], n2[k]])
    }
    d <- sapply(c(1:length(n1)), site_d, simplify=TRUE)
    return(d)
}


GetPairwiseNucleotideWeightSingleSite <- function(d1, d2, Ne, ci, diploid){
    if(diploid==TRUE){
        b = 1
    }else{
        b = 2
    }
    if(d1==d2){ #When the fitnesses are the same, neutral case, pure drift
        return(1/(2*Ne))
    }else{
        fit_ratio <- exp(-(d1-d2)*ci*Ne) #f1/f2
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
    nucleotide.distances <- CreateNucleotideDistanceMatrix()
    nucleotide.fitness.ratios <- matrix(data=0,4,4)
    unique.nucs <- n2s(0:3)
    for (i in sequence(4)) {
        for (j in sequence(4)) {
            nuc1 <- n2s(nucleotide.set[i])
            nuc2 <- n2s(nucleotide.set[j])
            if(!nuc1 == nuc2){
                d1 <- GetProteinProteinDistance(protein1=nuc1, protein2=unique.nucs[optimal.nucleotide], aa.distances=nucleotide.distances)
                d2 <- GetProteinProteinDistance(protein1=nuc2, protein2=unique.nucs[optimal.nucleotide], aa.distances=nucleotide.distances)
                nucleotide.fitness.ratios[i,j] <- GetPairwiseNucleotideWeightSingleSite(d1=d1, d2=d2, Ne=Ne, ci=position.multiplier, diploid=diploid)
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
    phy <- reorder(phy, "pruningwise")
    diag(nuc.mutation.rates) = 0
    diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
    scale.factor <- -sum(diag(nuc.mutation.rates) * root.p_array)
    nuc.mutation.rates_scaled <- nuc.mutation.rates_scaled * (1/scale.factor)

    for(site.index in sequence(nsites)) {
        weight.matrix <- GetNucleotideFixationMatrix(site.index, position.multiplier=position.multiplier.vector[site.index], optimal.nucleotide=nuc.optim_array[site.index], Ne=Ne, diploid=diploid)
        Q_position <- (ploidy * Ne) * nuc.mutation.rates_scaled * weight.matrix
        diag(Q_position) = 0
        diag(Q_position) <- -rowSums(Q_position)
        #scale.factor <- -sum(diag(Q_position) * root.p_array)
        final.likelihood.vector[site.index] <- GetLikelihoodUCEForSingleCharGivenOptimum(charnum=site.index, nuc.data=nuc.data, phy=phy, Q_position=Q_position, root.p=root.p_array, scale.factor=NULL, return.all=FALSE)
    }
    return(final.likelihood.vector)
}


PositionSensitivityMultiplierNormal <- function(a0, a1, a2, site.index){
    #At the moment assumes a standard normal distribution
    #sensitivity.vector <- (1/(a2*sqrt(2*pi))) * exp(-((site.index - a1)^2)/ (2*(a2^2)))
    sensitivity.vector <- a0 * exp(-((site.index - a1)^2)/ (2*(a2^2)))
    #sensitivity.vector <- a0 + a1*(site.index) + a2*((site.index)^2)
    return(sensitivity.vector)
}


PositionSensitivityMultiplierSigmoid <- function(slope.left, slope.right, midpoint, gene.length) {
    site.index.left <- 1:midpoint
    inflection.left <- midpoint/2
    sensitivity.vector.left <- 1/(1+exp(slope.left*inflection.left-slope.left*site.index.left))
    sensitivity.vector.left <- sensitivity.vector.left/max(sensitivity.vector.left)
    site.index.right <- gene.length:midpoint
    inflection.right <- gene.length - ((gene.length-midpoint)/2)
    sensitivity.vector.right <- 1/(1+exp(slope.right*inflection.right-slope.right*site.index.right))
    sensitivity.vector.right <- sensitivity.vector.right/max(sensitivity.vector.right)
    return(c(sensitivity.vector.left, sensitivity.vector.right))
}


PositionSensitivityMultiplierQuadratic <- function(a0, a1, a2, site.index){
    #At the moment assumes a standard normal distribution
    #sensitivity.vector <- (1/(a2*sqrt(2*pi))) * exp(-((site.index - a1)^2)/ (2*(a2^2)))
    sensitivity.vector <- a0 + a1*site.index + a2*(site.index^2)
    return(sensitivity.vector)
}


GetLikelihoodUCEForManyCharGivenAllParams <- function(x, nuc.data, phy, nuc.optim_array=NULL, nuc.model, diploid=TRUE, logspace=FALSE, verbose=TRUE, neglnl=FALSE) {
    if(logspace) {
        x = exp(x)
    }
    
    Ne=5e6
    x[1] <- x[1]/Ne
    if(nuc.model == "JC") {
        base.freqs=c(x[4:6], 1-sum(x[4:6]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "GTR") {
        base.freqs=c(x[4:6], 1-sum(x[4:6]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[7:length(x)], model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "UNREST") {
        tmp <- CreateNucleotideMutationMatrixSpecial(x[4:length(x)])
        base.freqs <- tmp$base.freqs
        nuc.mutation.rates <- tmp$nuc.mutation.rates
    }
    
    nsites <- dim(nuc.data)[2]-1
    site.index <- 1:nsites
    #Note that I am rescaling x[2] and x[3] so that I can optimize in log space, but also have negative slopes.
    #position.multiplier.vector <- x[1] * PositionSensitivityMultiplierSigmoid(x[2]+(-5), x[3]+(-5), x[4], nsites)
    position.multiplier.vector <- PositionSensitivityMultiplierNormal(x[1], x[2], x[3], site.index)
    final.likelihood = GetLikelihoodUCEForManyCharVaryingBySite(nuc.data=nuc.data, phy=phy, nuc.mutation.rates=nuc.mutation.rates, position.multiplier.vector=position.multiplier.vector, Ne=Ne, nuc.optim_array=nuc.optim_array, root.p_array=base.freqs, diploid=diploid)
    likelihood <- sum(final.likelihood)
    
    if(neglnl) {
        likelihood <- -1 * likelihood
    }
    if(is.na(likelihood)){
        return(1000000)
    }
    if(verbose) {
        results.vector <- c(likelihood)
        names(results.vector) <- c("likelihood")
        print(results.vector)
    }
    return(likelihood)
}


OptimizeEdgeLengthsUCE <- function(x, par.mat, site.pattern.data.list, n.partitions, nsites.vector, index.matrix, phy, nuc.optim.list=NULL, diploid=TRUE, nuc.model, logspace=FALSE, verbose=TRUE, n.cores=NULL, neglnl=FALSE) {
    if(logspace) {
        x <- exp(x)
    }
    
    phy$edge.length = x
    
    #sums the total number of parameters: 4 is the general shape pars, 3 are the base pars, and finally, the transition rates.
    if(nuc.model == "JC"){
        max.par = 3 + 3 + 0
    }
    if(nuc.model == "GTR"){
        max.par = 3 + 3 + 5
    }
    if(nuc.model == "UNREST"){
        max.par = 3 + 11
    }
    if(is.null(n.cores)){
        likelihood.vector <- c()
        for(partition.index in sequence(n.partitions)){
            nuc.data = NULL
            nuc.data = site.pattern.data.list[[partition.index]]
            likelihood.vector = c(likelihood.vector, GetLikelihoodUCEForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, nuc.optim_array=nuc.optim.list[[partition.index]], nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl))
        }
        likelihood = sum(likelihood.vector)
    }else{
        MultiCoreLikelihood <- function(partition.index){
            nuc.data = NULL
            nuc.data = site.pattern.data.list[[partition.index]]
            likelihood.tmp = GetLikelihoodUCEForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, nuc.optim_array=nuc.optim.list[[partition.index]], nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl)
            return(likelihood.tmp)
        }
        #This orders the nsites per partition in decreasing order (to increase efficiency):
        partition.order <- 1:n.partitions
        likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores)))
    }
    return(likelihood)
}



OptimizeModelParsUCE <- function(x, fixed.pars, site.pattern.data.list, n.partitions, nsites.vector, index.matrix, phy, nuc.optim.list=NULL, diploid=TRUE, nuc.model, logspace=FALSE, verbose=TRUE, n.cores=NULL, neglnl=FALSE, all.pars=FALSE) {
    
    if(logspace) {
        x <- exp(x)
    }
    
    #sums the total number of parameters: 4 is the general shape pars, 3 are the base pars, and finally, the transition rates.
    if(nuc.model == "JC"){
        max.par = 3 + 3 + 0
    }
    if(nuc.model == "GTR"){
        max.par = 3 + 3 + 5
    }
    if(nuc.model == "UNREST"){
        max.par = 3 + 11
    }
    
    if(all.pars == TRUE){
        if(class(index.matrix)=="numeric"){
            index.matrix <- matrix(index.matrix, 1, length(index.matrix))
        }
        par.mat <- index.matrix
        par.mat[] <- c(x, 0)[index.matrix]
    }else{
        par.mat <- matrix(c(x, fixed.pars), 1, max.par)
    }
    
    nuc.data = NULL
    nuc.data = site.pattern.data.list
    likelihood.vector = GetLikelihoodUCEForManyCharGivenAllParams(x=log(par.mat), nuc.data=nuc.data, phy=phy, nuc.optim_array=nuc.optim.list, nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl)
    likelihood = sum(likelihood.vector)
    return(likelihood)
}



OptimizeNucAllGenesUCE <- function(x, fixed.pars, site.pattern.data.list, n.partitions, nsites.vector, index.matrix, phy, nuc.optim.list=NULL, diploid=TRUE, nuc.model, logspace=FALSE, verbose=TRUE, n.cores=NULL, neglnl=FALSE) {
    if(logspace) {
        x <- exp(x)
    }
    #sums the total number of parameters: 4 is the general shape pars, 3 are the base pars, and finally, the transition rates.
    if(nuc.model == "JC"){
        max.par = 3 + 3 + 0
    }
    if(nuc.model == "GTR"){
        max.par = 3 + 3 + 5
    }
    if(nuc.model == "UNREST"){
        max.par = 3 + 11
    }
    
    par.mat <- c()
    for(row.index in 1:dim(fixed.pars)[1]){
        par.mat <- rbind(par.mat, c(fixed.pars[row.index,1:3], x))
    }
    MultiCoreLikelihood <- function(partition.index){
        nuc.data = NULL
        nuc.data = site.pattern.data.list[[partition.index]]
        likelihood.tmp = GetLikelihoodUCEForManyCharGivenAllParams(x=log(par.mat[partition.index,]), nuc.data=nuc.data, phy=phy, nuc.optim_array=nuc.optim.list[[partition.index]], nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl)
        return(likelihood.tmp)
    }
    #This orders the nsites per partition in decreasing order (to increase efficiency):
    partition.order <- 1:n.partitions
    likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores)))
    return(likelihood)
    
}


GetOptimalNucPerSite <- function(x, logspace=TRUE, verbose=FALSE, neglnl=TRUE){
    
    return(1)
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
        return(1000000)
    }
    else{
        loglik<- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
        if(is.infinite(loglik)){
            return(1000000)
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
#' @param edge.length Indicates whether or not edge lengths should be optimized. By default it is set to "optimize", other option is "fixed", which user-supplied branch lengths.
#' @param edge.linked A logical indicating whether or not edge lengths should be optimized separately for each gene. By default, a single set of each lengths is optimized for all genes.
#' @param optimal.nuc Indicates what type of optimal.nuc should be used. At the moment there is only a single option: "majrule".
#' @param nuc.model Indicates what type nucleotide model to use. There are three options: "JC", "GTR", or "UNREST".
#' @param global.nucleotide.model assumes nucleotide model is shared among all partitions
#' @param diploid A logical indicating whether or not the organism is diploid or not.
#' @param verbose Logical indicating whether each iteration be printed to the screen.
#' @param n.cores The number of cores to run the analyses over.
#' @param max.tol Supplies the relative optimization tolerance.
#' @param max.evals Supplies the max number of iterations tried during optimization.
#' @param max.restarts Supplies the number of random restarts.
#' @param output.by.restart Logical indicating whether or not each restart is saved to a file. Default is TRUE.
#' @param output.restart.filename Designates the file name for each random restart.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.
#'
#' @details
#' SELON stands for SELection On Nucleotides. This function takes a user supplied topology and a set of fasta formatted sequences and optimizes the parameters in the SELON model. Selection is based on selection towards an optimal nucleotide at each site, which is based simply on the majority rule of the observed data. The strength of selection is then varied along sites based on a Taylor series, which scales the substitution rates. Still a work in development, but so far, seems very promising.
SelonOptimize <- function(nuc.data.path, n.partitions=NULL, phy, edge.length="optimize", edge.linked=TRUE, optimal.nuc="majrule", nuc.model="GTR", global.nucleotide.model=TRUE, diploid=TRUE, verbose=FALSE, n.cores=1, max.tol=.Machine$double.eps^0.25, max.evals=1000000, max.restarts=10, output.by.restart=TRUE, output.restart.filename="restartResult", fasta.rows.to.keep=NULL) {
    
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
        starting.branch.lengths[partition.index,] <- ComputeStartingBranchLengths(phy, gene.tmp, data.type="dna",recalculate.starting.brlen=TRUE)$edge.length
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
    opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = max.evals, "ftol_rel" = max.tol)
    if(max.restarts > 1){
        selon.starting.vals <- matrix(0, max.restarts+1, 2)
        selon.starting.vals[,1] <- runif(n = max.restarts+1, min = (10^-20)*5e6, max = (10^-7)*5e6)
        #selon.starting.vals[,2] <- runif(n = max.restarts+1, min = 0.01, max = 10)
        selon.starting.vals[,2] <- runif(n = max.restarts+1, min = 0.01, max = 500)
    }else{
        selon.starting.vals <- matrix(c(1e-13*5e6, 100),1,2)
        selon.starting.vals <- rbind(selon.starting.vals, selon.starting.vals)
    }
    if(nuc.model == "JC"){
        ip = c(selon.starting.vals[1,1], ceiling(nsites.vector[1]/2), selon.starting.vals[1,2], 0.25, 0.25, 0.25)
        parameter.column.names <- c("s.Ne", "midpoint", "width", "freqA", "freqC", "freqG")
        upper = c(log(50), log(nsites.vector[1]), log(500), 0, 0, 0)
        lower = rep(-21, length(ip))
        max.par.model.count = 3 + 3 + 0
    }
    if(nuc.model == "GTR"){
        nuc.ip = rep(1, 5)
        ip = c(selon.starting.vals[1,1], ceiling(nsites.vector[1]/2), selon.starting.vals[1,2], 0.25, 0.25, 0.25, nuc.ip)
        parameter.column.names <- c("s.Ne", "midpoint", "width", "freqA", "freqC", "freqG", "C_A", "G_A", "T_A", "G_C", "T_C")
        upper = c(log(50), log(nsites.vector[1]), log(500), 0, 0, 0, rep(21, length(nuc.ip)))
        lower = rep(-21, length(ip))
        max.par.model.count = 3 + 3 + 5
    }
    if(nuc.model == "UNREST"){
        nuc.ip = rep(1, 11)
        ip = c(selon.starting.vals[1,1], ceiling(nsites.vector[1]/2), selon.starting.vals[1,2], nuc.ip)
        parameter.column.names <- c("s.Ne", "midpoint", "width", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T")
        upper = c(log(50), log(nsites.vector[1]), log(500), rep(21, length(nuc.ip)))
        lower = rep(-21, length(ip))
        max.par.model.count = 3 + 11
    }
    index.matrix = matrix(0, n.partitions, max.par.model.count)
    index.matrix[1,] = 1:ncol(index.matrix)
    ip.vector = ip
    if(n.partitions > 1){
        if(global.nucleotide.model == TRUE) {
            upper.mat = upper
            lower.mat = lower
            for(partition.index in 2:n.partitions){
                if(nuc.model == "JC"){
                    ip[2] = ceiling(nsites.vector[partition.index]/2)
                    upper[2] = log(nsites.vector[partition.index])
                    ip.vector = c(ip.vector, ip[1:3])
                    upper.mat = c(upper.mat, upper[1:3])
                    lower.mat = c(lower.mat, lower[1:3])
                }else{
                    if(nuc.model == "GTR"){
                        index.matrix.tmp = numeric(max.par.model.count)
                        ip[2] = ceiling(nsites.vector[partition.index]/2)
                        upper[2] = log(nsites.vector[partition.index])
                        index.matrix.tmp[c(4:11)] = c(4:11)
                        ip.vector = c(ip.vector, ip[1:3])
                        upper.mat = rbind(upper.mat, upper)
                        lower.mat = rbind(lower.mat, lower)
                        
                    }else{
                        index.matrix.tmp = numeric(max.par.model.count)
                        ip[2] = ceiling(nsites.vector[partition.index]/2)
                        upper[2] = log(nsites.vector[partition.index])
                        index.matrix.tmp[c(4:14)] = c(4:14)
                        ip.vector = c(ip.vector, ip[1:3])
                        upper.mat = rbind(upper.mat, upper)
                        lower.mat = rbind(lower.mat, lower)
                    }
                }
                index.matrix.tmp[index.matrix.tmp==0] <- seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
                index.matrix[partition.index,] <- index.matrix.tmp
            }
        }else{
            upper.vector = upper
            lower.vector = lower
            for(partition.index in 2:n.partitions){
                if(nuc.model == "JC"){
                    ip[2] = ceiling(nsites.vector[partition.index]/2)
                    upper[2] = log(nsites.vector[partition.index])
                    ip.vector = c(ip.vector, ip[1:3], ip[4], ip[5], ip[6])
                    upper.vector = c(upper.vector, c(upper[1:3], upper[4], upper[5], upper[6]))
                    lower.vector = c(lower.vector, c(lower[1:3], lower[4], lower[5], lower[6]))
                }else{
                    ip[2] = ceiling(nsites.vector[partition.index]/2)
                    upper[2] = log(nsites.vector[partition.index])
                    ip.vector = c(ip.vector, ip[1:3], ip[4], ip[5], ip[6], nuc.ip)
                    upper.vector = c(upper.vector, c(upper[1:3], upper[4], upper[5], upper[6], rep(21, length(nuc.ip))))
                    lower.vector = c(lower.vector, c(lower[1:3], lower[4], lower[5], lower[6], rep(-21, length(nuc.ip))))
                }
                index.matrix.tmp = numeric(max.par.model.count)
                index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
                index.matrix[partition.index,] <- index.matrix.tmp
            }
        }
    }
    number.of.current.restarts <- 1
    nuc.optim.original <- nuc.optim.list
    best.lik <- 1000000
    while(number.of.current.restarts < (max.restarts+1)){
        cat(paste("Finished. Performing random restart ", number.of.current.restarts,"...", sep=""), "\n")
        if(edge.length == "optimize"){
            phy$edge.length <- colMeans(starting.branch.lengths)
        }
        nuc.optim.list <- nuc.optim.list
        cat("       Doing first pass using majority rule optimal amino acid...", "\n")
        if(edge.length == "optimize"){
            cat("              Optimizing edge lengths", "\n")
            mle.pars.mat <- index.matrix
            mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
            plot(phy)
            print(phy$edge.length)
            print(mle.pars.mat)
            opts.edge <- opts
            upper.edge <- rep(log(5), length(phy$edge.length))
            lower.edge <- rep(log(1e-8), length(phy$edge.length))
            results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengthsUCE, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=nuc.optim.list, diploid=diploid, nuc.model=nuc.model, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
            phy$edge.length <- exp(results.edge.final$solution)
        }
        if(global.nucleotide.model == TRUE) {
            cat("              Optimizing model parameters", "\n")
            substitution.pars <- mle.pars.mat[1,c(4:max.par.model.count)]
            #Part 1: Optimize the shape function:
            index.matrix.red <- t(matrix(1:(n.partitions*3), 3, n.partitions))
            ParallelizedOptimizedByGene <- function(n.partition){
                tmp.par.mat <- as.matrix(mle.pars.mat[,1:3])
                upper.bounds.gene <- upper.mat[n.partition, 1:3]
                lower.bounds.gene <- lower.mat[n.partition, 1:3]
                optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, fixed.pars=substitution.pars, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=nuc.optim.list[[n.partition]], diploid=diploid, nuc.model=nuc.model, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=FALSE)
                tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
                return(tmp.pars)
            }
            results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores)
            parallelized.parameters <- t(matrix(unlist(results.set), 4, n.partitions))
            results.final <- NULL
            results.final$objective <- sum(parallelized.parameters[,1])
            results.final$solution <- c(t(parallelized.parameters[,-1]))
            mle.pars.mat.red <- index.matrix.red
            mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]
            upper.bounds.shared <- upper[c(4:max.par.model.count)]
            lower.bounds.shared <- lower[c(4:max.par.model.count)]
            optim.substitution.pars.by.gene <- nloptr(x0=log(substitution.pars), eval_f = OptimizeNucAllGenesUCE, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=nuc.optim.list, diploid=diploid, nuc.model=nuc.model, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
            results.final$objective <- optim.substitution.pars.by.gene$objective
            substitution.pars <- exp(optim.substitution.pars.by.gene$solution)
            mle.pars.mat <- c()
            for(row.index in 1:dim(mle.pars.mat.red)[1]){
                mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1:3], substitution.pars))
            }
        }else{
            cat("              Optimizing model parameters", "\n")
            # Optimize it all!
            ParallelizedOptimizedByGene <- function(n.partition){
                optim.by.gene <- nloptr(x0=log(mle.pars.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.vector, lb=lower.vector, opts=opts, fixed.pars=NULL, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=nuc.optim.list[[n.partition]], diploid=diploid, nuc.model=nuc.model, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=TRUE)
                tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
                return(tmp.pars)
            }
            
            results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores)
            parallelized.parameters <- t(matrix(unlist(results.set), 2, n.partitions))
            results.final <- NULL
            results.final$objective <- sum(parallelized.parameters[,1])
            results.final$solution <- c(t(parallelized.parameters[,-1]))
            mle.pars.mat <- index.matrix
            mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
            lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
            current.likelihood <- results.final$objective
            cat(paste("       Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
            iteration.number <- iteration.number + 1
        }
        
        current.likelihood <- results.final$objective
        cat(paste("       Current likelihood", current.likelihood, sep=" "), "\n")
        lik.diff <- 10
        iteration.number <- 1
        while(lik.diff != 0 & iteration.number<7){
            cat(paste("       Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
            
            if(optimal.nuc == "optimize"){
                cat("              Optimizing nucleotide", "\n")
                nuc.optim.list <- as.list(numeric(n.partitions))
                for(partition.index in sequence(n.partitions)) {
                    gene.tmp <- read.dna(partitions[partition.index], format='fasta')
                    if(!is.null(fasta.rows.to.keep)){
                        gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
                    }else{
                        gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
                    }
                    nucleotide.data <- DNAbinToNucleotideNumeric(gene.tmp)
                    nucleotide.data <- nucleotide.data[phy$tip.label,]
                    nuc.optim.list[[partition.index]] = GetOptimalNucPerSite(x=log(mle.pars.mat[partition.index,]), logspace=TRUE, verbose=verbose, neglnl=TRUE)
                }
            }
            if(edge.length == "optimize"){
                cat("              Optimizing edge lengths", "\n")
                results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengthsUCE, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=nuc.optim.list, diploid=diploid, nuc.model=nuc.model, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
                phy$edge.length <- exp(results.edge.final$solution)
            }
            if(global.nucleotide.model == TRUE) {
                cat("              Optimizing model parameters", "\n")
                substitution.pars <- mle.pars.mat[1,c(4:max.par.model.count)]
                #Part 1: Optimize the shape function:
                index.matrix.red <- t(matrix(1:(n.partitions*3), 3, n.partitions))
                ParallelizedOptimizedByGene <- function(n.partition){
                    tmp.par.mat <- as.matrix(mle.pars.mat[,1:3])
                    upper.bounds.gene <- upper.mat[n.partition, 1:3]
                    lower.bounds.gene <- lower.mat[n.partition, 1:3]
                    optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, fixed.pars=substitution.pars, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=nuc.optim.list[[n.partition]], diploid=diploid, nuc.model=nuc.model, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=FALSE)
                    tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
                    return(tmp.pars)
                }
                results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores)
                parallelized.parameters <- t(matrix(unlist(results.set), 4, n.partitions))
                results.final <- NULL
                results.final$objective <- sum(parallelized.parameters[,1])
                results.final$solution <- c(t(parallelized.parameters[,-1]))
                mle.pars.mat.red <- index.matrix.red
                mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]
                upper.bounds.shared <- upper[c(4:max.par.model.count)]
                lower.bounds.shared <- lower[c(4:max.par.model.count)]
                optim.substitution.pars.by.gene <- nloptr(x0=log(substitution.pars), eval_f = OptimizeNucAllGenesUCE, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=nuc.optim.list, diploid=diploid, nuc.model=nuc.model, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
                results.final$objective <- optim.substitution.pars.by.gene$objective
                substitution.pars <- exp(optim.substitution.pars.by.gene$solution)
                mle.pars.mat <- c()
                for(row.index in 1:dim(mle.pars.mat.red)[1]){
                    mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1:3], substitution.pars))
                }
                lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
                current.likelihood <- results.final$objective
                cat(paste("       Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
                iteration.number <- iteration.number + 1
            }else{
                cat("              Optimizing model parameters", "\n")
                # Optimize it all!
                ParallelizedOptimizedByGene <- function(n.partition){
                    optim.by.gene <- nloptr(x0=log(mle.pars.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.vector, lb=lower.vector, opts=opts, fixed.pars=NULL, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=nuc.optim.list[[n.partition]], diploid=diploid, nuc.model=nuc.model, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=TRUE)
                    tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
                    return(tmp.pars)
                }
                
                results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores)
                parallelized.parameters <- t(matrix(unlist(results.set), 2, n.partitions))
                results.final <- NULL
                results.final$objective <- sum(parallelized.parameters[,1])
                results.final$solution <- c(t(parallelized.parameters[,-1]))
                mle.pars.mat <- index.matrix
                mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
                lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
                current.likelihood <- results.final$objective
                cat(paste("       Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
                iteration.number <- iteration.number + 1
            }
            
        }
        
        if(output.by.restart == TRUE){
            obj.tmp = list(np=max(index.matrix) + length(phy$edge.lengths) + sum(nsites.vector), loglik = results.final$objective, AIC = -2*results.final$objective+2*(max(index.matrix) + length(phy$edge.lengths) + sum(nsites.vector)), AICc = NULL, mle.pars=mle.pars.mat, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, nuc.optim=nuc.optim.list, nuc.optim.type=optimal.nuc, nuc.model=nuc.model, diploid=diploid, empirical.base.freqs=empirical.base.freq.list, max.tol=max.tol, max.tol=max.tol, max.evals=max.evals, selon.starting.vals=ip.vector)
            class(obj.tmp) = "selac"
            save(obj.tmp,file=paste(paste(nuc.data.path, output.restart.filename, sep=""), number.of.current.restarts, "Rsave", sep="."))
        }
        if(results.final$objective < best.lik){
            best.ip <- ip.vector
            best.lik <- results.final$objective
            best.solution <- mle.pars.mat
            best.edge.lengths <- phy$edge.length
            best.nuc.optim.list <- nuc.optim.list
        }
        number.of.current.restarts <- number.of.current.restarts + 1
        ip.vector[c(index.matrix[,1])] <- selon.starting.vals[number.of.current.restarts, 1]
        ip.vector[2] <- c(selon.starting.vals[number.of.current.restarts, 2])
        nuc.optim.list <- nuc.optim.original
    }
    selon.starting.vals <- best.ip
    loglik <- -(best.lik) #to go from neglnl to lnl
    mle.pars.mat <- best.solution
    nuc.optim.list <- best.nuc.optim.list
    if(edge.length == "optimize"){
        phy$edge.length <- best.edge.lengths
    }
    cat("Finished. Summarizing results...", "\n")
    colnames(mle.pars.mat) <- parameter.column.names
    
    if(edge.length == "optimize"){
        np <- max(index.matrix) + length(phy$edge.lengths) + sum(nsites.vector)
    }else{
        np <- max(index.matrix) + sum(nsites.vector)
    }
    obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, AICc = NULL, mle.pars=mle.pars.mat, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, nuc.optim=nuc.optim.list, nuc.optim.type=optimal.nuc, nuc.model=nuc.model, diploid=diploid, empirical.base.freqs=empirical.base.freq.list, max.tol=max.tol, max.tol=max.tol, max.evals=max.evals, selon.starting.vals=selon.starting.vals)
    class(obj) = "selon"
    return(obj)
}






######################################################################################################################################
######################################################################################################################################
### Print function for the selon class:
######################################################################################################################################
######################################################################################################################################

print.selon <- function(x,...){
    ntips=Ntip(x$phy)
    output<-data.frame(x$loglik,x$AIC,ntips,sum(x$nsites), row.names="")
    names(output)<-c("-lnL","AIC", "ntax", "nsites")
    cat("\nFit\n")
    print(output)
}


######################################################################################################################################
######################################################################################################################################
### TESTING CODE
######################################################################################################################################
######################################################################################################################################

#Stuff to do:
#1. Expand matrix to allow optimal nucleotide to vary along branches.
#2. Add in 1/100 observations for when you have invariants sites.
#3. Topological inference.
#4. Add in normal nucleotide models -- remove from selac essentially.

#tree <- read.tree("/Users/jbeaulieu/Desktop/selonSims5genes/rokasYeast.tre")
#phy <- drop.tip(tree, "Calb")
#gene.tmp <- read.dna("/Users/jbeaulieu/selac/R/fastaSet_1/gene1.fasta", format='fasta')
#gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
#nucleotide.data <- DNAbinToNucleotideNumeric(gene.tmp)
#nucleotide.data <- nucleotide.data[phy$tip.label,]
#empirical.base.freq <- as.matrix(nucleotide.data[,-1])
#empirical.base.freq <- table(empirical.base.freq, deparse.level = 0) / sum(table(empirical.base.freq, deparse.level = 0))
#nuc.optim <- as.numeric(apply(nucleotide.data[,-1], 2, GetMaxNameUCE)) #starting values for all, final values for majrule
#GetLikelihoodUCEForManyCharGivenAllParams(x=log(c(4.26e-6, ceiling(235/2), 10, .25, .25, .25)), nuc.data=nucleotide.data, phy=phy, nuc.optim_array=nuc.optim, nuc.model="JC", Ne=5e6, diploid=TRUE, logspace=TRUE)
#kk <- PositionSensitivityMultiplierNormal(4.26e-6, ceiling(235/2), 10, site.index<-seq(0, 235, by=1))
#nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model="JC", base.freqs=rep(.25,4))
#res <- c()
#nucleotide.data.red <- nucleotide.data[,c(1:20)]
#nuc.optim.red <- nuc.optim[1:19]
#for(i in 1:length(kk)){
#kk.red <- rep(kk[i], length(kk))
#   res <- rbind(res, GetLikelihoodUCEForManyCharVaryingBySite(nuc.data=nucleotide.data.red, phy=phy, nuc.mutation.rates=nuc.mutation.rates, position.multiplier.vector=kk.red, Ne=5e6, nuc.optim_array=nuc.optim.red, root.p_array=NULL, diploid=TRUE))
#}

#pdf("ends.pdf")
#plot(res[,2])
#dev.off()
#nsite.length=300
#optimal.nuc <- sample(1:4, nsite.length, replace=TRUE)
#save(optimal.nuc.list, file="optimal.nuc.list.Rsave")
#gtr <- rep(1,5)
#base.freq <- c(.25,.25,.25)
#par.vec = c(1e-7,1,1, 150, base.freq, gtr)
#tmp <- SelonSimulator(phy=phy, pars=par.vec, nuc.optim_array=optimal.nuc, nuc.model="GTR", diploid=TRUE)
#system(paste("mkdir", paste("fastaSet",1, sep="_"), sep=" "))
#write.dna(tmp, file=paste(paste("fastaSet",1, sep="_"), "/gene",  1, ".fasta", sep=""), format="fasta", colw=1000000)

#pp <- SelonOptimize("", n.partitions=NULL, phy=phy2, edge.length="optimize", edge.linked=TRUE, optimal.nuc="majrule", nuc.model="JC", diploid=TRUE, n.cores=3)




