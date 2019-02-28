
######################################################################################################################################
######################################################################################################################################
### Sella and Hirsh model for nucleotides
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu


.nucleotide.name <- c("a", "c", "g", "t")


rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
}

######################################################################################################################################
######################################################################################################################################
### Functions for canonical implementation
######################################################################################################################################
######################################################################################################################################


CreateNucleotideDistanceMatrix <- function() {
    n.states <- 4
    nucleotide.distances <- matrix(1,nrow=n.states,ncol=n.states)
    diag(nucleotide.distances) <- 0
    rownames(nucleotide.distances) <- colnames(nucleotide.distances) <- .nucleotide.name
    return(nucleotide.distances)
}


CreateNucleotideMutationMatrixSpecial <- function(rates) {
    index <- matrix(NA, 4, 4)
    np <- 12
    index[col(index) != row(index)] <- 1:np
    nuc.mutation.rates <- matrix(0, nrow=4, ncol=4)
    nuc.mutation.rates<-matrix(rates[index], dim(index))
    rownames(nuc.mutation.rates) <- .nucleotide.name
    colnames(nuc.mutation.rates) <- .nucleotide.name
    nuc.mutation.rates[3,4] = 1
    diag(nuc.mutation.rates) <- 0
    diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
    #Next we take our rates and find the homogeneous solution to Q*pi=0 to determine the base freqs:
    base.freqs <- Null(nuc.mutation.rates)
    #Rescale base.freqs so that they sum to 1:
    base.freqs.scaled <- c(base.freqs/sum(base.freqs))
    base.freqs.scaled.matrix <- rep.row(base.freqs.scaled, 4)
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
    unique.nucs <- .nucleotide.name
    for (i in sequence(4)) {
        for (j in sequence(4)) {
            nuc1 <- .nucleotide.name[i]
            nuc2 <- .nucleotide.name[j]
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
        state <- nuc.data[i,charnum+1]
        if(state < 65){
            liks[i,state] <- 1
        }else{
            #If here, then the site has no data, so we treat it as ambiguous for all possible codons. Likely things might be more complicated, but this can be modified later:
            liks[i,] <- 1
        }
    }
    #The result here is just the likelihood:
    result <- -GetLikelihood(phy=phy, liks=liks, Q=Q_position, root.p=root.p)
    #ODE way is commented out
    #Q_position_vectored <- c(t(Q_position)) # has to be transposed
    #result <- -TreeTraversalSelonODE(phy=phy, Q_codon_array_vectored=Q_position_vectored, liks.HMM=liks, bad.likelihood=-100000, root.p=root.p)
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
    nuc.mutation.rates_scaled <- nuc.mutation.rates * (1/scale.factor)
    for(site.index in sequence(nsites)) {
        weight.matrix <- GetNucleotideFixationMatrix(site.index, position.multiplier=position.multiplier.vector[site.index], optimal.nucleotide=nuc.optim_array[site.index], Ne=Ne, diploid=diploid)
        Q_position <- (ploidy * Ne) * nuc.mutation.rates_scaled * weight.matrix
        diag(Q_position) <- 0
        diag(Q_position) <- -rowSums(Q_position)
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
        base.freqs <- tmp$base.freq
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


GetOptimalNucPerSite <- function(x, nuc.data, phy, nuc.model, diploid=TRUE, logspace=TRUE, verbose=FALSE, neglnl=TRUE){
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
        base.freqs <- tmp$base.freq
        nuc.mutation.rates <- tmp$nuc.mutation.rates
    }
    
    nsites <- dim(nuc.data)[2]-1
    site.index <- 1:nsites
    optimal.vector.by.site <- rep(NA, nsites)
    #unique.aa <- GetMatrixAANames(numcode)
    optimal.nuc.likelihood.mat <- matrix(0, nrow=4, ncol=nsites)
    position.multiplier.vector <- PositionSensitivityMultiplierNormal(x[1], x[2], x[3], site.index)
    for(i in 1:4){
        nuc.optim_array = rep(i, nsites)
        tmp = GetLikelihoodUCEForManyCharVaryingBySite(nuc.data=nuc.data, phy=phy, nuc.mutation.rates=nuc.mutation.rates, position.multiplier.vector=position.multiplier.vector, Ne=Ne, nuc.optim_array=nuc.optim_array, root.p_array=base.freqs, diploid=diploid)
        tmp[is.na(tmp)] = -1000000
        final.likelihood = tmp
        optimal.nuc.likelihood.mat[i,] <- final.likelihood
    }
    for(j in 1:nsites){
        optimal.vector.by.site[j] <- which.is.max(optimal.nuc.likelihood.mat[,j])
    }
    return(optimal.vector.by.site)
}



######################################################################################################################################
######################################################################################################################################
### Functions for HMM implementation
######################################################################################################################################
######################################################################################################################################

CreateHMMNucleotideMutationMatrix <- function(model="UNREST", rates, base.freqs=NULL, opt.nucleotide.transition) {
    mat.dim <- 4*4
    evolv.nucleotide.mutation <- matrix(data=0, nrow=mat.dim, ncol=mat.dim)

    if(model=="UNREST"){
        individual.matrix <- CreateNucleotideMutationMatrixHMMSpecial(rates=rates, model=model)
        for(i in 1:4) {
            index.vec.diag.i <- (1+(i-1)*4):(4+(i-1)*4)
            evolv.nucleotide.mutation[index.vec.diag.i, index.vec.diag.i] <- individual.matrix
            for(j in 1:4){
                index.vec.diag.j <- (1+(j-1)*4):(4+(j-1)*4)
                diag(evolv.nucleotide.mutation[index.vec.diag.j, index.vec.diag.i]) <- opt.nucleotide.transition
            }
        }
        diag(evolv.nucleotide.mutation) <- 0
        diag(evolv.nucleotide.mutation) <- -rowSums(evolv.nucleotide.mutation)
        #Next we take our rates and find the homogeneous solution to Q*pi=0 to determine the base freqs:
        base.freqs <- Null(evolv.nucleotide.mutation)
        #Rescale base.freqs so that they sum to 1:
        base.freqs.scaled <- c(base.freqs/sum(base.freqs))
        base.freqs.scaled.matrix <- rep.row(base.freqs, mat.dim)
        diag(evolv.nucleotide.mutation) <- 0
        #Rescale Q to account for base.freqs:
        evolv.nucleotide.mutation <- evolv.nucleotide.mutation * base.freqs.scaled.matrix
        diag(evolv.nucleotide.mutation) <- -rowSums(evolv.nucleotide.mutation)

        rownames(evolv.nucleotide.mutation) <- paste(rep(.nucleotide.name, times=4), rep(.nucleotide.name, each=4), sep="")
        colnames(evolv.nucleotide.mutation) <- paste(rep(.nucleotide.name, times=4), rep(.nucleotide.name, each=4), sep="")
        
        obj <- NULL
        obj$base.freq <- base.freqs.scaled
        obj$nuc.mutation.rates <- evolv.nucleotide.mutation
        
        return(obj)
    }else{
        individual.matrix <- CreateNucleotideMutationMatrixHMMSpecial(rates=rates, model=model)
        for(i in 1:4) {
            index.vec.diag.i <- (1+(i-1)*4):(4+(i-1)*4)
            evolv.nucleotide.mutation[index.vec.diag.i, index.vec.diag.i] <- individual.matrix
            for(j in 1:4){
                index.vec.diag.j <- (1+(j-1)*4):(4+(j-1)*4)
                diag(evolv.nucleotide.mutation[index.vec.diag.j, index.vec.diag.i]) <- opt.nucleotide.transition
            }
        }
        diag(evolv.nucleotide.mutation) <- 0
        diag(evolv.nucleotide.mutation) <- -rowSums(evolv.nucleotide.mutation)
        if(!is.null(base.freqs)){
            diag(evolv.nucleotide.mutation) = 0
            evolv.nucleotide.mutation = t(evolv.nucleotide.mutation * base.freqs)
            diag(evolv.nucleotide.mutation) = -rowSums(evolv.nucleotide.mutation)
        }
        rownames(evolv.nucleotide.mutation) <- paste(rep(.nucleotide.name, times=4), rep(.nucleotide.name, each=4), sep="")
        colnames(evolv.nucleotide.mutation) <- paste(rep(.nucleotide.name, times=4), rep(.nucleotide.name, each=4), sep="")
        
        return(evolv.nucleotide.mutation)
    }
}


CreateNucleotideMutationMatrixHMMSpecial <- function(rates, model="UNREST") {
    if(model == "JC") {
        nuc.mutation.rates <- matrix(data=rates[1], nrow=4, ncol=4)
        diag(nuc.mutation.rates) <- 0
        diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
        return(nuc.mutation.rates)
    }
    if(model == "GTR") {
        index <- matrix(NA, 4, 4)
        np <- 5
        sel <- col(index) < row(index)
        sel[4,3] = FALSE
        index[sel] <- 1:np
        index <- t(index)
        index[sel] <- 1:np
        nuc.mutation.rates <- matrix(0, nrow=4, ncol=4)
        nuc.mutation.rates<-matrix(rates[index], dim(index))
        nuc.mutation.rates[4,3] <- nuc.mutation.rates[3,4] <- 1
        diag(nuc.mutation.rates) <- 0
        diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
        return(nuc.mutation.rates)
    }
    if(model == "UNREST"){
        index <- matrix(NA, 4, 4)
        np <- 12
        index[col(index) != row(index)] <- 1:np
        nuc.mutation.rates <- matrix(0, nrow=4, ncol=4)
        nuc.mutation.rates <- matrix(rates[index], dim(index))
        nuc.mutation.rates[3,4] = 1
        diag(nuc.mutation.rates) <- 0
        diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
        return(nuc.mutation.rates)
    }
}


GetNucleotideFixationHMMMatrix <- function(site.number, position.multiplier, Ne, diploid=TRUE){
    nucleotide.set <- 0:3
    mat.dim <- 4*4
    evolv.nucleotide.fixation.probs <- matrix(data=0, nrow=mat.dim, ncol=mat.dim)
    
    for(i in 1:4) {
        index.vec.diag.i <- (1+(i-1)*4):(4+(i-1)*4)
        evolv.nucleotide.fixation.probs[index.vec.diag.i, index.vec.diag.i] <- GetNucleotideFixationMatrix(site.number=site.number, position.multiplier=position.multiplier, optimal.nucleotide=i, Ne=Ne, diploid=diploid)
        for(j in 1:4){
            index.vec.diag.j <- (1+(j-1)*4):(4+(j-1)*4)
            diag(evolv.nucleotide.fixation.probs[index.vec.diag.j, index.vec.diag.i]) <- 1/(2*Ne)
        }
    }
    rownames(evolv.nucleotide.fixation.probs) <- paste(rep(.nucleotide.name, times=4), rep(.nucleotide.name, each=4), sep="")
    colnames(evolv.nucleotide.fixation.probs) <- paste(rep(.nucleotide.name, times=4), rep(.nucleotide.name, each=4), sep="")

    return(evolv.nucleotide.fixation.probs)
}


GetLikelihoodUCEHMMForSingleCharGivenOptimum <- function(charnum=1, nuc.data, phy, Q_position, root.p=NULL, scale.factor, return.all=FALSE) {
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
            if(nuc.data[i,charnum+1] == 1){
                liks[i,c(1,5,9,13)] <- 1
            }
            if(nuc.data[i,charnum+1] == 2){
                liks[i,c(2,6,10,14)] <- 1
            }
            if(nuc.data[i,charnum+1] == 3){
                liks[i,c(3,7,11,15)] <- 1
            }
            if(nuc.data[i,charnum+1] == 4){
                liks[i,c(4,8,12,16)] <- 1
            }
        }else{
            #If here, then the site has no data, so we treat it as ambiguous for all possible codons. Likely things might be more complicated, but this can be modified later:
            liks[i,] <- 1
        }
    }
    #print(Q_position)
    #print(liks)
    #print(root.p)
    #The result here is just the likelihood:
    result <- -GetLikelihood(phy=phy, liks=liks, Q=Q_position, root.p=root.p)
    #ODE way is commented out
    #Q_position_vectored <- c(t(Q_position)) # has to be transposed
    #result <- -TreeTraversalSelonODE(phy=phy, Q_codon_array_vectored=Q_position_vectored, liks.HMM=liks, bad.likelihood=-100000, root.p=root.p)
    ifelse(return.all, stop("return all not currently implemented"), return(result))
}


GetLikelihoodHMMUCEForManyCharVaryingBySite <- function(nuc.data, phy, nuc.mutation.rates, position.multiplier.vector, Ne, root.p_array=NULL, diploid=TRUE) {
    nsites <- dim(nuc.data)[2]-1
    final.likelihood.vector <- rep(NA, nsites)
    if(is.null(root.p_array)) {
########REMAINING ISSUE -- not clear on the frequencies under HMM. Recaled to normalize to 1, or not?
        #Generate matrix of equal frequencies for each site:
        root.p_array <- rep(1/16, 16)
        root.p_array <- root.p_array/sum(root.p_array)
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
    nuc.mutation.rates_scaled <- nuc.mutation.rates * (1/scale.factor)
    for(site.index in sequence(nsites)) {
        weight.matrix <- GetNucleotideFixationHMMMatrix(site.index, position.multiplier=position.multiplier.vector[site.index], Ne=Ne, diploid=diploid)
        Q_position <- (ploidy * Ne) * nuc.mutation.rates_scaled * weight.matrix
        diag(Q_position) <- 0
        diag(Q_position) <- -rowSums(Q_position)
        final.likelihood.vector[site.index] <- GetLikelihoodUCEHMMForSingleCharGivenOptimum(charnum=site.index, nuc.data=nuc.data, phy=phy, Q_position=Q_position, root.p=root.p_array, scale.factor=NULL, return.all=FALSE)
    }
    return(final.likelihood.vector)
}


GetLikelihoodUCEHMMForManyCharGivenAllParams <- function(x, nuc.data, phy, nuc.optim_array=NULL, nuc.model, diploid=TRUE, logspace=FALSE, verbose=TRUE, neglnl=FALSE) {
    if(logspace) {
        x = exp(x)
    }
    Ne=5e6
    x[1] <- x[1]/Ne
    if(nuc.model == "JC") {
########REMAINING ISSUE -- not clear on the frequencies under HMM. Recaled to normalize to 1, or not?
        opt.nucleotide.transition = x[7]
        base.freqs=c(x[4:6], 1-sum(x[4:6]))
        base.freqs=rep(c(x[4:6], 1-sum(x[4:6])),4)
        base.freqs=base.freqs/sum(base.freqs)
        nuc.mutation.rates <- CreateHMMNucleotideMutationMatrix(model=nuc.model, rates=1, base.freqs=base.freqs, opt.nucleotide.transition=opt.nucleotide.transition)
    }
    if(nuc.model == "GTR") {
########REMAINING ISSUE -- not clear on the frequencies under HMM. Recaled to normalize to 1, or not?
        base.freqs=rep(c(x[4:6], 1-sum(x[4:6])),4)
        base.freqs=base.freqs/sum(base.freqs)
        opt.nucleotide.transition = x[7]
        nuc.mutation.rates <- CreateHMMNucleotideMutationMatrix(model=nuc.model, rates=x[8:length(x)], base.freqs=base.freqs, opt.nucleotide.transition=opt.nucleotide.transition)
    }
    if(nuc.model == "UNREST") {
        opt.nucleotide.transition = x[4]
        tmp <- CreateHMMNucleotideMutationMatrix(model=nuc.model, rates=x[5:length(x)], base.freqs=NULL, opt.nucleotide.transition=opt.nucleotide.transition)
        base.freqs <- tmp$base.freq
        nuc.mutation.rates <- tmp$nuc.mutation.rates
    }
    
    nsites <- dim(nuc.data)[2]-1
    site.index <- 1:nsites
    #Note that I am rescaling x[2] and x[3] so that I can optimize in log space, but also have negative slopes.
    #position.multiplier.vector <- x[1] * PositionSensitivityMultiplierSigmoid(x[2]+(-5), x[3]+(-5), x[4], nsites)
    position.multiplier.vector <- PositionSensitivityMultiplierNormal(x[1], x[2], x[3], site.index)
    final.likelihood <- GetLikelihoodHMMUCEForManyCharVaryingBySite(nuc.data=nuc.data, phy=phy, nuc.mutation.rates=nuc.mutation.rates, position.multiplier.vector=position.multiplier.vector, Ne=Ne, root.p_array=base.freqs, diploid=diploid)
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



######################################################################################################################################
######################################################################################################################################
### Functions used by everything
######################################################################################################################################
######################################################################################################################################


#OptimizeEdgeLengthsUCE <- function(x, par.mat, site.pattern.data.list, n.partitions, nsites.vector, index.matrix, phy, nuc.optim.list=NULL, diploid=TRUE, nuc.model, hmm=FALSE, logspace=FALSE, verbose=TRUE, n.cores=NULL, neglnl=FALSE) {
    
#    if(logspace) {
#        x <- exp(x)
#    }
#
#    phy$edge.length = x
#    if(hmm==TRUE){
#        #sums the total number of parameters: 4 is the general shape pars, 3 are the base pars, 1 for transition rate among hidden nucleotides, and finally, the transition rates.
#        if(nuc.model == "JC"){
#            max.par = 3 + 3 + 1 + 0
#        }
#        if(nuc.model == "GTR"){
#            max.par = 3 + 3 + 1 + 5
#        }
#        if(nuc.model == "UNREST"){
#            max.par = 3 + 1 + 11
#        }
#        #print("here")
#        if(is.null(n.cores)){
#            likelihood.vector <- c()
#            for(partition.index in sequence(n.partitions)){
#                nuc.data = NULL
#                nuc.data = site.pattern.data.list[[partition.index]]
#                likelihood.vector = c(likelihood.vector, GetLikelihoodUCEHMMForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl))
#            }
#            likelihood = sum(likelihood.vector)
#        }else{
#            MultiCoreLikelihood <- function(partition.index){
#                nuc.data = NULL
#                nuc.data = site.pattern.data.list[[partition.index]]
#                likelihood.tmp = GetLikelihoodUCEHMMForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl)
#                return(likelihood.tmp)
#            }
#            #This orders the nsites per partition in decreasing order (to increase efficiency):
#            partition.order <- 1:n.partitions
#            likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores)))
#        }
#        return(likelihood)
#
#    }else{
#        #sums the total number of parameters: 4 is the general shape pars, 3 are the base pars, and finally, the transition rates.
#        if(nuc.model == "JC"){
#            max.par = 3 + 3 + 0
#        }
#        if(nuc.model == "GTR"){
#            max.par = 3 + 3 + 5
#        }
#        if(nuc.model == "UNREST"){
#            max.par = 3 + 11
#        }
#
#        if(is.null(n.cores)){
#            likelihood.vector <- c()
#            for(partition.index in sequence(n.partitions)){
#                nuc.data = NULL
#                nuc.data = site.pattern.data.list[[partition.index]]
#                likelihood.vector = c(likelihood.vector, GetLikelihoodUCEForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, nuc.optim_array=nuc.optim.list[[partition.index]], nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl))
#            }
#            likelihood = sum(likelihood.vector)
#        }else{
#            MultiCoreLikelihood <- function(partition.index){
#                nuc.data = NULL
#                nuc.data = site.pattern.data.list[[partition.index]]
#                likelihood.tmp = GetLikelihoodUCEForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, nuc.optim_array=nuc.optim.list[[partition.index]], nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl)
#                return(likelihood.tmp)
#            }
#            #This orders the nsites per partition in decreasing order (to increase efficiency):
#            partition.order <- 1:n.partitions
#            likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores)))
#        }
#        return(likelihood)
#    }
#}


OptimizeModelParsUCE <- function(x, fixed.pars, site.pattern.data.list, n.partitions, nsites.vector, index.matrix, phy, nuc.optim.list=NULL, diploid=TRUE, nuc.model, hmm=FALSE, logspace=FALSE, verbose=TRUE, n.cores=NULL, neglnl=FALSE, all.pars=FALSE) {
    
    if(logspace) {
        x <- exp(x)
    }
    
    if(hmm==TRUE){
        #sums the total number of parameters: 4 is the general shape pars, 3 are the base pars, 1 for transition rate among hidden nucleotides, and finally, the transition rates.
        if(nuc.model == "JC"){
            max.par = 3 + 3 + 1 + 0
        }
        if(nuc.model == "GTR"){
            max.par = 3 + 3 + 1 + 5
        }
        if(nuc.model == "UNREST"){
            max.par = 3 + 1 + 11
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
        likelihood.vector = GetLikelihoodUCEHMMForManyCharGivenAllParams(x=log(par.mat), nuc.data=nuc.data, phy=phy, nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl)
        likelihood = sum(likelihood.vector)
        return(likelihood)
    }else{
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
}


OptimizeNucAllGenesUCE <- function(x, fixed.pars, site.pattern.data.list, n.partitions, nsites.vector, index.matrix, phy, nuc.optim.list=NULL, diploid=TRUE, nuc.model, hmm=FALSE, logspace=FALSE, verbose=TRUE, n.cores=NULL, neglnl=FALSE) {
    if(logspace) {
        x <- exp(x)
    }
    
    if(hmm == TRUE){
        #sums the total number of parameters: 4 is the general shape pars, 3 are the base pars, 1 for transition rate among hidden nucleotides, and finally, the transition rates.
        if(nuc.model == "JC"){
            max.par = 3 + 3 + 1 + 0
        }
        if(nuc.model == "GTR"){
            max.par = 3 + 3 + 1 + 5
        }
        if(nuc.model == "UNREST"){
            max.par = 3 + 1 + 11
        }
        
        par.mat <- c()
        for(row.index in 1:dim(fixed.pars)[1]){
            par.mat <- rbind(par.mat, c(fixed.pars[row.index,1:3], x))
        }
        MultiCoreLikelihood <- function(partition.index){
            nuc.data = NULL
            nuc.data = site.pattern.data.list[[partition.index]]
            likelihood.tmp = GetLikelihoodUCEHMMForManyCharGivenAllParams(x=log(par.mat[partition.index,]), nuc.data=nuc.data, phy=phy, nuc.model=nuc.model, diploid=diploid, logspace=logspace, verbose=verbose, neglnl=neglnl)
            return(likelihood.tmp)
        }
        #This orders the nsites per partition in decreasing order (to increase efficiency):
        partition.order <- 1:n.partitions
        likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores)))
        return(likelihood)
    }else{
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
}



######################################################################################################################################
######################################################################################################################################
### Edge Length Optimizer
######################################################################################################################################
######################################################################################################################################

FindBranchGenerations <- function(phy) {
    generation <- list()
    known <- 1:Ntip(phy)
    unknown <- phy$edge[,1]
    needed <- phy$edge[,2]
    root <- min(unknown)
    i <- 1
    repeat{
        knowable <- unknown[needed %in% known]
        knowable <- knowable[duplicated(knowable)]
        generation[[i]] <- phy$edge[which(phy$edge[,1] %in% knowable),2]
        known <- c(known, knowable)
        needed <- needed[!unknown %in% knowable]
        unknown <- unknown[!unknown %in% knowable]
        i <- i + 1
        if (any(root == knowable)) break
    }
    res <- generation
    return(res)
}


#Generates a DataArray for all sites -- slow but no prohibitively so. Unclear if data.table is necessary.
MakeDataArray <- function(site.pattern.data.list, phy, nsites.vector) {
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    #Now we need to build the matrix of likelihoods to pass to dev.raydisc:
    #liks.array <- array(data=0, dim=c(nb.tip+nb.node, nl, sum(nsites.vector)))
    site.pattern.data.frame <- site.pattern.data.list[[1]]
    #Now loop through the tips.
    for(partition.index in 2:length(site.pattern.data.list)){
        site.pattern.data.frame <- cbind(site.pattern.data.frame, site.pattern.data.list[[partition.index]][,2:(nsites.vector[partition.index]+1)])
    }
    site.pattern.data.table <- as.data.table(site.pattern.data.frame[,-1])
    site.pattern.data.table <- t(site.pattern.data.table)
    colnames(site.pattern.data.table) <- phy$tip.label
    return(site.pattern.data.table)
}


#Goal is to make a list with all the things I need for a site.
MakeParameterArray <- function(nuc.optim.list, pars.mat, nsites.vector) {
    Ne <- 5e6
    pars.array <- c()
    for(partition.index in 1:length(nsites.vector)){
        pars.site.tmp <- as.list(1:nsites.vector[partition.index])
        site.index <- 1:nsites.vector[partition.index]
        position.multiplier.vector <- PositionSensitivityMultiplierNormal(pars.mat[partition.index,1]/Ne, pars.mat[partition.index,2], pars.mat[partition.index,3], site.index)
        for(site.index in 1:nsites.vector[partition.index]) {
            pars.site.tmp[[site.index]] <- c(nuc.optim.list[[partition.index]][site.index], position.multiplier.vector[site.index], pars.mat[partition.index,4:dim(pars.mat)[2]])
        }
        pars.array <- append(pars.array, pars.site.tmp)
    }
    return(pars.array)
}


#Goal is to make a list with all the things I need for a site.
MakeParameterArrayGTR <- function(site.pattern.count.list, empirical.base.freq.list, pars.mat, nsites.vector) {
        pars.array <- c()
        for(partition.index in 1:length(nsites.vector)){
            pars.site.tmp <- as.list(1:nsites.vector[partition.index])
            site.index <- 1:nsites.vector[partition.index]
            for(site.index in 1:nsites.vector[partition.index]) {
                pars.site.tmp[[site.index]] <- c(site.pattern.count.list[[partition.index]][site.index], empirical.base.freq.list[[partition.index]], pars.mat[partition.index,1:dim(pars.mat)[2]])
            }
            pars.array <- append(pars.array, pars.site.tmp)
        }
    return(pars.array)
}



#Will conduct single branch calculations
SingleBranchCalculation <- function(Q, init.cond, edge.length, root.p) {
    BranchProbs <- expm(Q * edge.length, method="Ward77") %*% init.cond
    if(is.nan(BranchProbs) || is.na(BranchProbs)){
        return(1000000)
    }
    return(BranchProbs)
}


GetBranchLikeAcrossAllSites <- function(p, edge.number, phy, data.array, pars.array, nuc.model, diploid, n.cores, logspace) {
    
    if(diploid == TRUE){
        ploidy <- 2
        Ne <- 5e6
    }else{
        ploidy <- 1
        Ne <- 5e6
    }
    
    if(logspace == TRUE){
        p <- exp(p)
    }
    if(!is.null(edge.number)){
        phy$edge.length[which(phy$edge[,2]==edge.number)] <- p
    }
    
    phy <- reorder(phy, "pruningwise")
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode

    MultiCoreLikelihood <- function(site.index, phy){
        # Parse parameters #
        x <- pars.array[[site.index]]
        optim.nuc <- x[1]
        x <- x[-1]
        position.multiplier <- x[1]
        x <- x[-1]
        ####################
        
        if(nuc.model == "JC") {
            base.freqs=c(x[1:3], 1-sum(x[1:3]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "GTR") {
            base.freqs=c(x[1:3], 1-sum(x[1:3]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[4:length(x)], model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "UNREST") {
            tmp <- CreateNucleotideMutationMatrixSpecial(x[1:length(x)])
            base.freqs <- tmp$base.freq
            nuc.mutation.rates <- tmp$nuc.mutation.rates
        }
        
        diag(nuc.mutation.rates) <- 0
        diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
        scale.factor <- -sum(diag(nuc.mutation.rates) * base.freqs)
        nuc.mutation.rates_scaled <- nuc.mutation.rates * (1/scale.factor)
        weight.matrix <- GetNucleotideFixationMatrix(site.index, position.multiplier=position.multiplier, optimal.nucleotide=optim.nuc, Ne=Ne, diploid=diploid)
        Q_position <- (ploidy * Ne) * nuc.mutation.rates_scaled * weight.matrix
        diag(Q_position) <- 0
        diag(Q_position) <- -rowSums(Q_position)
        
        liks <- matrix(0, nb.tip + nb.node, dim(Q_position)[1])
        for(i in 1:Ntip(phy)){
            state <- data.array[site.index,phy$tip.label[i]]
            if(state < 65){
                liks[i,state] <- 1
            }else{
                #If here, then the site has no data, so we treat it as ambiguous for all possible codons. Likely things might be more complicated, but this can be modified later:
                liks[i,] <- 1
            }
        }
        
        branchLikPerSite <- GetLikelihood(phy=phy, liks=liks, Q=Q_position, root.p=base.freqs)
        return(branchLikPerSite)
    }
    site.order <- 1:dim(data.array)[1]
    branchLikAllSites <- sum(unlist(mclapply(site.order, MultiCoreLikelihood, phy=phy, mc.cores=n.cores)))
    return(sum(branchLikAllSites))
}



GetBranchLikeAcrossAllSitesGTR <- function(p, edge.number, phy, data.array, pars.array, nuc.model, include.gamma, gamma.type, ncats, n.cores, logspace) {

    if(logspace == TRUE){
        p <- exp(p)
    }
    if(!is.null(edge.number)){
        phy$edge.length[which(phy$edge[,2]==edge.number)] <- p
    }
    
    phy <- reorder(phy, "pruningwise")
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    
    MultiCoreLikelihood <- function(site.index, phy){
        
        # Parse parameters #
        x <- pars.array[[site.index]]
        site.pattern.count <- x[1]
        x <- x[-1]
        base.freqs <- x[1:4]
        x <- x[-c(1:4)]
        if(include.gamma == TRUE){
            shape <- x[1]
            x <- x[-1]
        }
        ####################
        if(nuc.model == "JC") {
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "GTR") {
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[1:length(x)], model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "UNREST") {
            tmp <- CreateNucleotideMutationMatrixSpecial(x[1:length(x)])
            nuc.mutation.rates <- tmp$nuc.mutation.rates
        }
        diag(nuc.mutation.rates) = 0
        diag(nuc.mutation.rates) = -rowSums(nuc.mutation.rates)
        scale.factor <- -sum(diag(nuc.mutation.rates) * base.freqs)
        
        Q <- nuc.mutation.rates * (1/scale.factor)
        liks <- matrix(0, nb.tip + nb.node, dim(Q)[1])
        for(i in 1:Ntip(phy)){
            state <- data.array[site.index,phy$tip.label[i]]
            if(state < 65){
                liks[i,state] <- 1
            }else{
                #If here, then the site has no data, so we treat it as ambiguous for all possible codons. Likely things might be more complicated, but this can be modified later:
                liks[i,] <- 1
            }
        }

        if(include.gamma==TRUE){
            if(gamma.type == "median"){
                rates.k <- DiscreteGamma(shape=shape, ncats=ncats)
                weights.k <- rep(1/ncats, ncats)
            }
            if(gamma.type == "quadrature"){
                rates.and.weights <- LaguerreQuad(shape=shape, ncats=ncats)
                rates.k <- rates.and.weights[1:ncats]
                weights.k <- rates.and.weights[(ncats+1):(ncats*2)]
            }
            if(gamma.type == "lognormal"){
                rates.and.weights <- LogNormalQuad(shape=shape, ncats=ncats)
                rates.k <- rates.and.weights[1:ncats]
                weights.k <- rates.and.weights[(ncats+1):(ncats*2)]
            }
            tmp <- c()
            for(k in sequence(ncats)){
                tmp <- c(tmp, GetLikelihood(phy=phy, liks=liks, Q=(Q * rates.k[k]), root.p=base.freqs))
            }
            branchLikPerSite <- -log(sum(exp(-tmp) * weights.k)) * site.pattern.count
        }else{
            branchLikPerSite <- GetLikelihood(phy=phy, liks=liks, Q=Q, root.p=base.freqs) * site.pattern.count
        }
        return(branchLikPerSite)
    }
    site.order <- 1:dim(data.array)[1]
    branchLikAllSites <- sum(unlist(mclapply(site.order, MultiCoreLikelihood, phy=phy, mc.cores=n.cores)))
    return(sum(branchLikAllSites))
}

#out <- optimize(GetBranchLikeAcrossAllSites, edge.number=generations[[gen.index]][index], phy=phy, data.array=data.array, pars.array=pars.array, nuc.model=nuc.model, diploid=TRUE, lower=log(1e-8), upper=log(10), maximum=FALSE, tol = .Machine$double.eps^0.25)


## Go by independent generations. As we get deeper and deeper in the tree, we have to do less of the traversal. Needs: To update data matrix as we go down and to ignore edges we have already ML'd.
## Step 1: Send appropriate info to SingleBranch calculation to get right info based on new MLE of branch we just evaluated
## Step 2: Replace row info, across each site. Issue though is that we'd have to regenerate data.array after we're done? Actually no because basically once we done a single round we're done here.
OptimizeEdgeLengthsUCENew <- function(phy, pars.mat, site.pattern.data.list, nuc.optim.list, nuc.model, nsites.vector, diploid, logspace, n.cores, neglnl=FALSE) {
    
    maxit <- 11
    tol <- .Machine$double.eps^0.25
    nb.tip <- Ntip(phy)
    nb.node <- Nnode(phy)
    TIPS <- 1:nb.tip
    generations <- FindBranchGenerations(phy)
    data.array <- MakeDataArray(site.pattern.data.list=site.pattern.data.list, phy=phy, nsites.vector=nsites.vector)
    pars.array <- MakeParameterArray(nuc.optim.list=nuc.optim.list, pars.mat=pars.mat, nsites.vector=nsites.vector)
    are_we_there_yet <- 1
    iteration.number <- 1
    old.likelihood <- GetBranchLikeAcrossAllSites(p=phy$edge.length, edge.number=NULL, phy=phy, data.array=data.array, pars.array=pars.array, nuc.model=nuc.model, diploid=diploid, n.cores=n.cores, logspace=logspace)
    while (are_we_there_yet > tol && iteration.number < maxit) {
        cat("                   Round number",  iteration.number, "\n")
        for(gen.index in 1:length(generations)){
            for(index in 1:length(generations[[gen.index]])){
                cat("                        Optimizing edge number",  generations[[gen.index]][index],"\n")
                out <- optimize(GetBranchLikeAcrossAllSites, edge.number=generations[[gen.index]][index], phy=phy, data.array=data.array, pars.array=pars.array, nuc.model=nuc.model, diploid=diploid, n.cores=n.cores, logspace=logspace, lower=1e-8, upper=10, maximum=FALSE, tol=tol)
                phy$edge.length[which(phy$edge[,2]==generations[[gen.index]][index])] <- out$minimum
            }
        }
        new.likelihood <- out$objective
        iteration.number <- iteration.number + 1
        are_we_there_yet <- (old.likelihood - new.likelihood ) / new.likelihood
        old.likelihood <- new.likelihood
    }
    
    final.likelihood <- out$objective
    if(neglnl) {
        final.likelihood <- -1 * final.likelihood
    }
    tree_and_likelihood <- NULL
    tree_and_likelihood$final.likelihood <- final.likelihood
    tree_and_likelihood$phy <- phy
    return(tree_and_likelihood)
}

#ppp <- OptimizeEdgeLengthsUCENew(phy=phy, pars.mat=pars.mat, site.pattern.data.list=site.pattern.data.list, nuc.optim.list=nuc.optim.list, nuc.model=nuc.model, nsites.vector=nsites.vector, diploid=TRUE, logspace=FALSE, n.cores=n.cores, neglnl=TRUE)

## Go by independent generations. As we get deeper and deeper in the tree, we have to do less of the traversal. Needs: To update data matrix as we go down and to ignore edges we have already ML'd.
## Step 1: Send appropriate info to SingleBranch calculation to get right info based on new MLE of branch we just evaluated
## Step 2: Replace row info, across each site. Issue though is that we'd have to regenerate data.array after we're done? Actually no because basically once we done a single round we're done here.
OptimizeEdgeLengthsGTRNew <- function(phy, pars.mat, site.pattern.data.list, site.pattern.count.list, empirical.base.freq.list, nuc.model, include.gamma, gamma.type, ncats, nsites.vector, logspace, n.cores, neglnl=FALSE) {
    
    maxit <- 11
    tol <- .Machine$double.eps^0.25
    nb.tip <- Ntip(phy)
    nb.node <- Nnode(phy)
    TIPS <- 1:nb.tip
    generations <- FindBranchGenerations(phy)

    nsites.vector.update <- c()
    for(partition.index in 1:length(site.pattern.data.list)){
        nsites.vector.update <- c(nsites.vector.update, length(site.pattern.count.list[[partition.index]]))
    }
    
    data.array <- MakeDataArray(site.pattern.data.list=site.pattern.data.list, phy=phy, nsites.vector=nsites.vector.update)
    pars.array <- MakeParameterArrayGTR(site.pattern.count.list=site.pattern.count.list, empirical.base.freq.list=empirical.base.freq.list, pars.mat=pars.mat, nsites.vector=nsites.vector.update)

    are_we_there_yet <- 1
    iteration.number <- 1
    old.likelihood <- GetBranchLikeAcrossAllSitesGTR(p=phy$edge.length, edge.number=NULL, phy=phy, data.array=data.array, pars.array=pars.array, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, n.cores=n.cores, logspace=logspace)
    while (are_we_there_yet > tol && iteration.number < maxit) {
        cat("                   Round number",  iteration.number, "\n")
        for(gen.index in 1:length(generations)){
            for(index in 1:length(generations[[gen.index]])){
                cat("                        Optimizing edge number",  generations[[gen.index]][index],"\n")
                out <- optimize(GetBranchLikeAcrossAllSitesGTR, edge.number=generations[[gen.index]][index], phy=phy, data.array=data.array, pars.array=pars.array, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, n.cores=n.cores, logspace=logspace, lower=1e-8, upper=10, maximum=FALSE, tol=tol)
                phy$edge.length[which(phy$edge[,2]==generations[[gen.index]][index])] <- out$minimum
            }
        }
        new.likelihood <- out$objective
        iteration.number <- iteration.number + 1
        are_we_there_yet <- (old.likelihood - new.likelihood) / new.likelihood
        old.likelihood <- new.likelihood
    }
    
    final.likelihood <- out$objective
    if(neglnl) {
        final.likelihood <- -1 * final.likelihood
    }
    tree_and_likelihood <- NULL
    tree_and_likelihood$final.likelihood <- final.likelihood
    tree_and_likelihood$phy <- phy
    return(tree_and_likelihood)
}



######################################################################################################################################
######################################################################################################################################
### Likelihood calculator -- One step process because at each site Q changes
######################################################################################################################################
######################################################################################################################################

GetLikelihood <- function(phy, liks, Q, root.p){
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
        desRows <- which(phy$edge[,1] == focal)
        desNodes <- phy$edge[desRows,2]
        v <- 1
        for (desIndex in sequence(length(desRows))){
            v <- v * internal_expmt(Q, phy$edge.length[desRows[desIndex]])[[1]] %*% liks[desNodes[desIndex],]
            #v <- v * expm(Q * phy$edge.length[desRows[desIndex]]) %*% liks[desNodes[desIndex],]
        }
        comp[focal] <- sum(v)
        liks[focal,] <- v/comp[focal]
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



######################################################################################################################################
######################################################################################################################################
### Likelihood calculator -- ODE solver
######################################################################################################################################

#TreeTraversalSelonODE <- function(phy, Q_codon_array_vectored, liks.HMM, bad.likelihood=-100000, root.p) {
    
    ##start with first method and move to next if problems encountered
    ## when solving ode, such as negative pr values < neg.pr.threshold
#    ode.method.vec <- c("lsoda", "ode45")
#    num.ode.method <- length(ode.method.vec)
    
#    rtol = 1e-6 #default 1e-6 returns a negative value under long branch testing conditions
#    atol = 1e-6 #default 1e-6
    
#    neg.pr.threshold <- -10*atol
    
#    nb.tip <- length(phy$tip.label)
#    nb.node <- phy$Nnode
    
#    anc <- unique(phy$edge[,1])
#    TIPS <- 1:nb.tip
    
#    comp <- numeric(nb.tip + nb.node)
    
#    for (i in seq(from = 1, length.out = nb.node)) {
#        focal <- anc[i]
#        desRows <- which(phy$edge[,1]==focal) ##des = descendant
#        desNodes <- phy$edge[desRows,2]
#        state.pr.vector = rep(1, dim(liks.HMM)[2]) ##
        
#        for (desIndex in sequence(length(desRows))){
#            yini <- liks.HMM[desNodes[desIndex],]
#            times=c(0, phy$edge.length[desRows[desIndex]])

#            ode.not.solved <- TRUE
#            ode.solver.attempt <- 0
            
#            while(ode.not.solved && ode.solver.attempt < num.ode.method){
#                ode.solver.attempt <- ode.solver.attempt+1
#                ode.method <-  ode.method.vec[ode.solver.attempt]
                
#                subtree.pr.ode.obj <- lsoda(
#                y=yini, times=times, func = "selon_ode",
#                parms=Q_codon_array_vectored, initfunc="initmod_selon",
#                dllname = "selonODE",
#                rtol=rtol, atol=atol
#                )
                
                ## CHECK TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
                ## $istate should be = 0 [documentation in doc/deSolve.Rnw indicates
                ## it should be 2]
                ## Values < 0 indicate problems
                ## TODO: take advantage of while() around ode solving created
                ## for when we hit negative values
#               istate <- attributes(subtree.pr.ode.obj)$istate[1]
                
#               if(istate < 0){
                    ## For \code{lsoda, lsodar, lsode, lsodes, vode, rk, rk4, euler} these are
#                   error.text <- switch(as.character(istate),
#                   "-1"="excess work done",
#                   "-2"="excess accuracy requested",
#                   "-3"="illegal input detected",
#                   "-4"="repeated error test failures",
#                   "-5"="repeated convergence failures",
#                   "-6"="error weight became zero",
#                   paste("unknown error. ode() istate value: ", as.character(istate))
#                   )
#                   warning(print(paste("selac.R: Integration of desIndex", desIndex, " ode solver returned istate[1] = ",  istate, " : ", error.text, " returning bad.likelihood")))
#                   return(bad.likelihood)
#               }else{
                    ##no integration issues,
                    ## object consists of pr values at start and end time
                    ## extract final state variable, dropping time entry
#                   subtree.pr.vector <- subtree.pr.ode.obj[dim(subtree.pr.ode.obj)[[1]],-1]
#               }
                
                ## test for negative entries
                ## if encountered and less than neg.pr.threshold
                ## replace the negative values to 0
                ## if there are values less than neg.pr.threshold, then
                ## resolve equations using more robust method on the list
                ## Alternative: use 'event' option in deSolve as described at
                ## http://stackoverflow.com/questions/34424716/using-events-in-desolve-to-prevent-negative-state-variables-r
#               neg.vector.pos <- which(subtree.pr.vector < 0, arr.ind=TRUE)
#               num.neg.vector.pos <- length(neg.vector.pos)
                
#               if(num.neg.vector.pos > 0){
#                   min.vector.val <- min(subtree.pr.vector[neg.vector.pos])
#                   neg.vector.pos.as.string <- toString(neg.vector.pos)
                    
#                   warning.message <- paste("WARNING: subtree.pr.vector solved with ode method ", ode.method, " contains ", num.neg.vector.pos, " negative values at positions ", neg.vector.pos.as.string ,  "of a ", length(subtree.pr.vector), " vector." )

#                   if(min.vector.val > neg.pr.threshold){
#                       warning.message <- paste(warning.message, "\nMinimum value ", min.vector.val, " >  ", neg.pr.threshold, " the neg.pr.threshold.\nSetting all negative values to 0.")
#                       warning(warning.message)
#                       subtree.pr.vector[neg.vector.pos] <- 0
                        
#                   }else{
#                       warning.message <- paste(warning.message, "selon.R: minimum value ", min.vector.val, " <  ", neg.pr.threshold, " the neg.pr.threshold.")
#
#                       if(ode.solver.attempt < num.ode.method){
#                           warning.message <- paste(warning.message, " Trying ode method ", ode.method.vec[ode.solver.attempt+1])
#                           warning(warning.message)
                            
#                       }else{
#                           warning.message <- paste(warning.message, "No additional ode methods available. Returning bad.likelihood: ", bad.likelihood)
#                           warning(warning.message)
#                           return(bad.likelihood)
#                       }
#                   }
#               }else{
                    ## no negative values in pr.vs.time.matrix
#                   ode.not.solved <- FALSE
#               }
                
#           } ##end while() for ode solver
#           state.pr.vector <- state.pr.vector * subtree.pr.vector
#       }
#       comp[focal] <- sum(state.pr.vector)
#       liks.HMM[focal,] <- state.pr.vector/comp[focal]
#   }
#   root.node <- nb.tip + 1L
    
    ##Check for negative transition rates
    ##mikeg:  For now, just issue warning
    
#   neg.nodes <- which(liks.HMM[root.node,] <0)
#   if(length(neg.nodes)>0){
#       warning(paste("selac.R: encountered " , length(neg.nodes), " negatives values in liks.HMM[", root.node, ", ", neg.nodes, " ] =  ",  liks.HMM[root.node, neg.nodes], " at position ", i, " , desIndex ", desIndex))
#   }
    
#   loglik <- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks.HMM[root.node,])))
    
    ##return bad.likelihood if loglik is bad
#   if(!is.finite(loglik)) return(bad.likelihood)

#    return(loglik)
#}


GetMaxNameUCE <- function(x) {
    x = unname(x)
    x = x[which(x!=65)]
    return(names(table(x))[(which.is.max(table(x)))]) #note that this breaks ties at random
}




######################################################################################################################################
######################################################################################################################################
### Organizes the optimization routine for canonical model
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
#' @param cycle.stage Specifies the number of cycles per restart. Default is 12.
#' @param max.restarts Supplies the number of random restarts.
#' @param output.by.restart Logical indicating whether or not each restart is saved to a file. Default is TRUE.
#' @param output.restart.filename Designates the file name for each random restart.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.
#'
#' @details
#' SELON stands for SELection On Nucleotides. This function takes a user supplied topology and a set of fasta formatted sequences and optimizes the parameters in the SELON model. Selection is based on selection towards an optimal nucleotide at each site, which is based simply on the majority rule of the observed data. The strength of selection is then varied along sites based on a Taylor series, which scales the substitution rates. Still a work in development, but so far, seems very promising.
SelonOptimize <- function(nuc.data.path, n.partitions=NULL, phy, edge.length="optimize", edge.linked=TRUE, optimal.nuc="majrule", nuc.model="GTR", global.nucleotide.model=TRUE, diploid=TRUE, verbose=FALSE, n.cores=1, max.tol=.Machine$double.eps^0.25, max.evals=1000000, cycle.stage=12, max.restarts=3, output.by.restart=TRUE, output.restart.filename="restartResult", fasta.rows.to.keep=NULL) {
    
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
        selon.starting.vals[,1] <- runif(n = max.restarts+1, min = (10^-20)*5e6, max = (10^-12)*5e6)
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
            phy$edge.length <- apply(starting.branch.lengths, 2, weighted.mean, nsites.vector)
        }
        nuc.optim.list <- nuc.optim.list
        cat("       Doing first pass using majority rule optimal nucleotide...", "\n")
        if(edge.length == "optimize"){
            cat("              Optimizing edge lengths", "\n")
            mle.pars.mat <- index.matrix
            mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
            phy$edge.length[phy$edge.length < 1e-08] <- 1e-08
            results.edge.final <- OptimizeEdgeLengthsUCENew(phy=phy, pars.mat=mle.pars.mat, site.pattern.data.list=site.pattern.data.list, nuc.optim.list=nuc.optim.list, nuc.model=nuc.model, nsites.vector=nsites.vector, diploid=diploid, logspace=FALSE, n.cores=n.cores, neglnl=TRUE)
            #results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengthsUCE, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=nuc.optim.list, diploid=diploid, nuc.model=nuc.model, hmm=FALSE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
            print(results.edge.final$final.likelihood)
            phy <- results.edge.final$phy
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
                optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, fixed.pars=substitution.pars, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=nuc.optim.list[[n.partition]], diploid=diploid, nuc.model=nuc.model, hmm=FALSE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=FALSE)
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
            optim.substitution.pars.by.gene <- nloptr(x0=log(substitution.pars), eval_f = OptimizeNucAllGenesUCE, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=nuc.optim.list, diploid=diploid, nuc.model=nuc.model, hmm=FALSE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
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
                optim.by.gene <- nloptr(x0=log(mle.pars.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.vector, lb=lower.vector, opts=opts, fixed.pars=NULL, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=nuc.optim.list[[n.partition]], diploid=diploid, nuc.model=nuc.model, hmm=FALSE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=TRUE)
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
        }
        
        current.likelihood <- results.final$objective
        cat(paste("       Current likelihood", current.likelihood, sep=" "), "\n")
        are_we_there_yet <- 1
        iteration.number <- 1
        while(are_we_there_yet > max.tol && iteration.number < (cycle.stage+1)){
            cat(paste("       Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
            
            if(optimal.nuc == "optimize"){
                cat("              Optimizing nucleotide", "\n")
                nuc.optim.list <- as.list(numeric(n.partitions))
                ParallelizedOptimizeNucByGene <- function(n.partition){
                    nuc.optim.list[[partition.index]] = GetOptimalNucPerSite(x=log(mle.pars.mat[n.partition,]), nuc.data=site.pattern.data.list[[n.partition]], phy=phy, nuc.model=nuc.model, diploid=diploid, logspace=TRUE, verbose=verbose, neglnl=TRUE)
                }
                nuc.optim.list <- mclapply(1:n.partitions, ParallelizedOptimizeNucByGene, mc.cores=n.cores)
            }
            if(edge.length == "optimize"){
                cat("              Optimizing edge lengths", "\n")
                phy$edge.length[phy$edge.length < 1e-08] <- 1e-08
                results.edge.final <- OptimizeEdgeLengthsUCENew(phy=phy, pars.mat=mle.pars.mat, site.pattern.data.list=site.pattern.data.list, nuc.optim.list=nuc.optim.list, nuc.model=nuc.model, nsites.vector=nsites.vector, diploid=diploid, logspace=FALSE, n.cores=n.cores, neglnl=TRUE)
                phy <- results.edge.final$phy
                #results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengthsUCE, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=nuc.optim.list, diploid=diploid, nuc.model=nuc.model, hmm=FALSE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
                #phy$edge.length <- exp(results.edge.final$solution)
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
                    optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, fixed.pars=substitution.pars, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=nuc.optim.list[[n.partition]], diploid=diploid, nuc.model=nuc.model, hmm=FALSE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=FALSE)
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
                optim.substitution.pars.by.gene <- nloptr(x0=log(substitution.pars), eval_f = OptimizeNucAllGenesUCE, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=nuc.optim.list, diploid=diploid, nuc.model=nuc.model, hmm=FALSE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
                results.final$objective <- optim.substitution.pars.by.gene$objective
                substitution.pars <- exp(optim.substitution.pars.by.gene$solution)
                mle.pars.mat <- c()
                for(row.index in 1:dim(mle.pars.mat.red)[1]){
                    mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1:3], substitution.pars))
                }
                are_we_there_yet <- (current.likelihood - results.final$objective ) / results.final$objective
                current.likelihood <- results.final$objective
                cat(paste("       Current likelihood", current.likelihood, sep=" "), paste("% difference from previous round", are_we_there_yet, sep=" "), "\n")
                iteration.number <- iteration.number + 1
            }else{
                cat("              Optimizing model parameters", "\n")
                # Optimize it all!
                ParallelizedOptimizedByGene <- function(n.partition){
                    optim.by.gene <- nloptr(x0=log(mle.pars.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.vector, lb=lower.vector, opts=opts, fixed.pars=NULL, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=nuc.optim.list[[n.partition]], diploid=diploid, nuc.model=nuc.model, hmm=FALSE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=TRUE)
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
                are_we_there_yet <- (current.likelihood - results.final$objective ) / results.final$objective
                current.likelihood <- results.final$objective
                cat(paste("       Current likelihood", current.likelihood, sep=" "), paste("% difference from previous round", are_we_there_yet, sep=" "), "\n")
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
        ip.vector[3] <- c(selon.starting.vals[number.of.current.restarts, 2])
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
### Organizes the optimization routine
######################################################################################################################################
######################################################################################################################################

#' @title Optimize parameters under the HMM SELON model
#'
#' @description
#' Optimizes model parameters under the HMM SELON model
#'
#' @param nuc.data.path Provides the path to the directory containing the gene specific fasta files that contains the nucleotide data.
#' @param n.partitions The number of partitions to analyze. The order is based on the Unix order of the fasta files in the directory.
#' @param phy The phylogenetic tree to optimize the model parameters.
#' @param edge.length Indicates whether or not edge lengths should be optimized. By default it is set to "optimize", other option is "fixed", which user-supplied branch lengths.
#' @param edge.linked A logical indicating whether or not edge lengths should be optimized separately for each gene. By default, a single set of each lengths is optimized for all genes.
#' @param nuc.model Indicates what type nucleotide model to use. There are three options: "JC", "GTR", or "UNREST".
#' @param global.nucleotide.model assumes nucleotide model is shared among all partitions
#' @param diploid A logical indicating whether or not the organism is diploid or not.
#' @param verbose Logical indicating whether each iteration be printed to the screen.
#' @param n.cores The number of cores to run the analyses over.
#' @param max.tol Supplies the relative optimization tolerance.
#' @param max.evals Supplies the max number of iterations tried during optimization.
#' @param cycle.stage Specifies the number of cycles per restart. Default is 12.
#' @param max.restarts Supplies the number of random restarts.
#' @param output.by.restart Logical indicating whether or not each restart is saved to a file. Default is TRUE.
#' @param output.restart.filename Designates the file name for each random restart.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.
#'
#' @details
#' SELON stands for SELection On Nucleotides. This function takes a user supplied topology and a set of fasta formatted sequences and optimizes the parameters in the SELON model. Selection is based on selection towards an optimal nucleotide at each site, which is based simply on the majority rule of the observed data. The strength of selection is then varied along sites based on a Taylor series, which scales the substitution rates. Still a work in development, but so far, seems very promising.
SelonHMMOptimize <- function(nuc.data.path, n.partitions=NULL, phy, edge.length="optimize", edge.linked=TRUE, nuc.model="GTR", global.nucleotide.model=TRUE, diploid=TRUE, verbose=FALSE, n.cores=1, max.tol=.Machine$double.eps^0.5, max.evals=1000000, cycle.stage=12, max.restarts=10, output.by.restart=TRUE, output.restart.filename="restartResult", fasta.rows.to.keep=NULL) {
    
    cat("Initializing data and model parameters...", "\n")
    
    partitions <- system(paste("ls -1 ", nuc.data.path, "*.fasta", sep=""), intern=TRUE)
    if(is.null(n.partitions)){
        n.partitions <- length(partitions)
    }else{
        n.partitions = n.partitions
    }
    site.pattern.data.list <- as.list(numeric(n.partitions))
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
    }
    opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = max.evals, "ftol_rel" = max.tol)
    
    if(max.restarts > 1){
        selon.starting.vals <- matrix(0, max.restarts+1, 2)
        selon.starting.vals[,1] <- runif(n = max.restarts+1, min = (10^-20)*5e6, max = (10^-12)*5e6)
        #selon.starting.vals[,2] <- runif(n = max.restarts+1, min = 0.01, max = 10)
        selon.starting.vals[,2] <- runif(n = max.restarts+1, min = 0.01, max = 500)
    }else{
        selon.starting.vals <- matrix(c(1e-13*5e6, 100),1,2)
        selon.starting.vals <- rbind(selon.starting.vals, selon.starting.vals)
    }

    if(nuc.model == "JC"){
        ip = c(selon.starting.vals[1,1], ceiling(nsites.vector[1]/2), selon.starting.vals[1,2], 0.25, 0.25, 0.25, 1)
        parameter.column.names <- c("s.Ne", "midpoint", "width", "freqA", "freqC", "freqG", "opt.rate")
        upper = c(log(50), log(nsites.vector[1]), log(500), 0, 0, 0, 21)
        lower = rep(-21, length(ip))
        max.par.model.count = 3 + 3 + 1 + 0
    }
    if(nuc.model == "GTR"){
        nuc.ip = rep(1, 5)
        ip = c(selon.starting.vals[1,1], ceiling(nsites.vector[1]/2), selon.starting.vals[1,2], 0.25, 0.25, 0.25, 1, nuc.ip)
        parameter.column.names <- c("s.Ne", "midpoint", "width", "freqA", "freqC", "freqG", "opt.rate", "C_A", "G_A", "T_A", "G_C", "T_C")
        upper = c(log(50), log(nsites.vector[1]), log(500), 0, 0, 0, 21, rep(21, length(nuc.ip)))
        lower = rep(-21, length(ip))
        max.par.model.count = 3 + 3 + 1 + 5
    }
    if(nuc.model == "UNREST"){
        nuc.ip = rep(1, 11)
        ip = c(selon.starting.vals[1,1], ceiling(nsites.vector[1]/2), selon.starting.vals[1,2], 1, nuc.ip)
        parameter.column.names <- c("s.Ne", "midpoint", "width", "opt.rate", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T")
        upper = c(log(50), log(nsites.vector[1]), log(500), 21, rep(21, length(nuc.ip)))
        lower = rep(-21, length(ip))
        max.par.model.count = 3 + 1 + 11
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
                        index.matrix.tmp[c(4:12)] = c(4:12)
                        ip.vector = c(ip.vector, ip[1:3])
                        upper.mat = rbind(upper.mat, upper)
                        lower.mat = rbind(lower.mat, lower)
                        
                    }else{
                        index.matrix.tmp = numeric(max.par.model.count)
                        ip[2] = ceiling(nsites.vector[partition.index]/2)
                        upper[2] = log(nsites.vector[partition.index])
                        index.matrix.tmp[c(4:15)] = c(4:15)
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
                    ip.vector = c(ip.vector, ip[1:3], ip[4], ip[5], ip[6], ip[7])
                    upper.vector = c(upper.vector, c(upper[1:3], upper[4], upper[5], upper[6], upper[7]))
                    lower.vector = c(lower.vector, c(lower[1:3], lower[4], lower[5], lower[6], lower[7]))
                }else{
                    ip[2] = ceiling(nsites.vector[partition.index]/2)
                    upper[2] = log(nsites.vector[partition.index])
                    ip.vector = c(ip.vector, ip[1:3], ip[4], ip[5], ip[6], ip[7], nuc.ip)
                    upper.vector = c(upper.vector, c(upper[1:3], upper[4], upper[5], upper[6], 21,  rep(21, length(nuc.ip))))
                    lower.vector = c(lower.vector, c(lower[1:3], lower[4], lower[5], lower[6], -21, rep(-21, length(nuc.ip))))
                }
                index.matrix.tmp = numeric(max.par.model.count)
                index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
                index.matrix[partition.index,] <- index.matrix.tmp
            }
        }
    }
    number.of.current.restarts <- 1
    best.lik <- 1000000
    while(number.of.current.restarts < (max.restarts+1)){
        cat(paste("Finished. Performing random restart ", number.of.current.restarts,"...", sep=""), "\n")
        if(edge.length == "optimize"){
            phy$edge.length <- colMeans(starting.branch.lengths)
        }
        cat("       Doing first pass HMM...", "\n")
        if(edge.length == "optimize"){
            cat("              Optimizing edge lengths", "\n")
            mle.pars.mat <- index.matrix
            mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
            phy$edge.length[phy$edge.length < 1e-08] <- 1e-08
            results.edge.final <- OptimizeEdgeLengthsUCENew(phy=phy, pars.mat=mle.pars.mat, site.pattern.data.list=site.pattern.data.list, nuc.optim.list=NULL, nuc.model=nuc.model, nsites.vector=nsites.vector, diploid=diploid, logspace=FALSE, n.cores=n.cores, neglnl=TRUE)
            #results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengthsUCE, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=nuc.optim.list, diploid=diploid, nuc.model=nuc.model, hmm=FALSE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
            print(results.edge.final$final.likelihood)
            phy <- results.edge.final$phy
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
                optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, fixed.pars=substitution.pars, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=NULL, diploid=diploid, nuc.model=nuc.model, hmm=TRUE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=FALSE)
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
            optim.substitution.pars.by.gene <- nloptr(x0=log(substitution.pars), eval_f = OptimizeNucAllGenesUCE, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=NULL, diploid=diploid, nuc.model=nuc.model, hmm=TRUE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
            results.final$objective <- optim.substitution.pars.by.gene$objective
            substitution.pars <- exp(optim.substitution.pars.by.gene$solution)
            mle.pars.mat <- c()
            for(row.index in 1:dim(mle.pars.mat.red)[1]){
                mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1:3], substitution.pars))
            }
            print(mle.pars.mat)
        }else{
            cat("              Optimizing model parameters", "\n")
            # Optimize it all!
            ParallelizedOptimizedByGene <- function(n.partition){
                optim.by.gene <- nloptr(x0=log(mle.pars.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.vector, lb=lower.vector, opts=opts, fixed.pars=NULL, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=NULL, diploid=diploid, nuc.model=nuc.model, hmm=TRUE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=TRUE)
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
        while(lik.diff != 0 & iteration.number < (cycle.stage+1)){
            cat(paste("       Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
            if(edge.length == "optimize"){
                cat("              Optimizing edge lengths", "\n")
                phy$edge.length[phy$edge.length < 1e-08] <- 1e-08
                results.edge.final <- OptimizeEdgeLengthsUCENew(phy=phy, pars.mat=mle.pars.mat, site.pattern.data.list=site.pattern.data.list, nuc.optim.list=NULL, nuc.model=nuc.model, nsites.vector=nsites.vector, diploid=diploid, logspace=FALSE, n.cores=n.cores, neglnl=TRUE)
                #results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengthsUCE, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=nuc.optim.list, diploid=diploid, nuc.model=nuc.model, hmm=FALSE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
                print(results.edge.final$final.likelihood)
                phy <- results.edge.final$phy
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
                    optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, fixed.pars=substitution.pars, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=NULL, diploid=diploid, nuc.model=nuc.model, hmm=TRUE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=FALSE)
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
                optim.substitution.pars.by.gene <- nloptr(x0=log(substitution.pars), eval_f = OptimizeNucAllGenesUCE, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, site.pattern.data.list=site.pattern.data.list, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, nuc.optim.list=NULL, diploid=diploid, nuc.model=nuc.model, hmm=TRUE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
                results.final$objective <- optim.substitution.pars.by.gene$objective
                substitution.pars <- exp(optim.substitution.pars.by.gene$solution)
                mle.pars.mat <- c()
                for(row.index in 1:dim(mle.pars.mat.red)[1]){
                    mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1:3], substitution.pars))
                }
                print(mle.pars.mat)
                lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
                current.likelihood <- results.final$objective
                cat(paste("       Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
                iteration.number <- iteration.number + 1
            }else{
                cat("              Optimizing model parameters", "\n")
                # Optimize it all!
                ParallelizedOptimizedByGene <- function(n.partition){
                    optim.by.gene <- nloptr(x0=log(mle.pars.mat[n.partition,]), eval_f = OptimizeModelParsUCE, ub=upper.vector, lb=lower.vector, opts=opts, fixed.pars=NULL, site.pattern.data.list=site.pattern.data.list[[n.partition]], n.partitions=n.partitions, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, nuc.optim.list=NULL, diploid=diploid, nuc.model=nuc.model, hmm=TRUE, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE, all.pars=TRUE)
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
            obj.tmp = list(np=max(index.matrix) + length(phy$edge.lengths) + sum(nsites.vector), loglik = results.final$objective, AIC = -2*results.final$objective+2*(max(index.matrix) + length(phy$edge.lengths) + sum(nsites.vector)), AICc = NULL, mle.pars=mle.pars.mat, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, nuc.model=nuc.model, diploid=diploid, empirical.base.freqs=empirical.base.freq.list, max.tol=max.tol, max.tol=max.tol, max.evals=max.evals, selon.starting.vals=ip.vector)
            class(obj.tmp) = "selac"
            save(obj.tmp,file=paste(paste(nuc.data.path, output.restart.filename, sep=""), number.of.current.restarts, "Rsave", sep="."))
        }
        if(results.final$objective < best.lik){
            best.ip <- ip.vector
            best.lik <- results.final$objective
            best.solution <- mle.pars.mat
            best.edge.lengths <- phy$edge.length
        }
        number.of.current.restarts <- number.of.current.restarts + 1
        ip.vector[c(index.matrix[,1])] <- selon.starting.vals[number.of.current.restarts, 1]
        ip.vector[3] <- c(selon.starting.vals[number.of.current.restarts, 2])
    }
    selon.starting.vals <- best.ip
    loglik <- -(best.lik) #to go from neglnl to lnl
    mle.pars.mat <- best.solution
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
    obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, AICc = NULL, mle.pars=mle.pars.mat, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, nuc.model=nuc.model, diploid=diploid, empirical.base.freqs=empirical.base.freq.list, max.tol=max.tol, max.tol=max.tol, max.evals=max.evals, selon.starting.vals=selon.starting.vals)
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




