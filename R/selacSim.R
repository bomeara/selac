
######################################################################################################################################
######################################################################################################################################
### SELAC simulator -- simulator for SELection on Amino acids and/or Codons model
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu

#source("selon.R")
#source("selac.R")

SingleSiteUpPass <- function(phy, Q_codon, root.value){
    #Randomly choose starting state at root using the root.values as the probability:
    root.value = sample.int(dim(Q_codon)[2], 1, FALSE, prob=root.value)
    #Reorder the phy:
    phy <- reorder(phy, "postorder")
    ntips <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- ntips + 1 #perhaps use an accessor to get the root node id
    #Generate vector that contains the simulated states:
    sim.codon.data.site <- integer(ntips + phy$Nnode)
    sim.codon.data.site[ROOT] <- as.integer(root.value)
    anc <- phy$edge[,1]
    des <- phy$edge[,2]
    edge.length <- phy$edge.length
    for (i in N:1) {
        p <- expm(Q_codon * edge.length[i], method="Ward77")[sim.codon.data.site[anc[i]], ]
        sim.codon.data.site[des[i]] <- sample.int(dim(Q_codon)[2], size = 1, FALSE, prob = p)
    }
    sim.codon.data.site <- sim.codon.data.site[1:ntips]
    return(sim.codon.data.site)
}


SingleSiteUpPassEvolvingRates <- function(phy, root.value, aa.distances, codon_mutation_matrix, aa.optim, nsites, C, Phi, q, Ne, numcode, diploid, pars.to.evolve="phi"){
    #Randomly choose starting state at root using the root.values as the probability:
    root.value = sample.int(dim(codon_mutation_matrix)[2], 1, FALSE, prob=root.value)
    #Reorder the phy:
    phy <- reorder(phy, "postorder")
    ntips <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- ntips + 1 #perhaps use an accessor to get the root node id
    unique.aa <- GetMatrixAANames(numcode=numcode)
    #Generate vector that contains the simulated states:
    sim.codon.data.site <- integer(ntips + phy$Nnode)
    sim.codon.data.site[ROOT] <- as.integer(root.value)
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    edge.length <- phy$edge.length
    for (i in N:1) {
        if(pars.to.evolve=="phi"){
            Phi.current <- Phi[des[i]]
            Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi.current, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
            for(k in 1:21){
                if(diploid == TRUE){
                    Q_codon_array[,,unique.aa[k]] = 2 * Ne * (codon_mutation_matrix * Q_codon_array[,,unique.aa[k]])
                }else{
                    Q_codon_array[,,unique.aa[k]] = Ne * (codon_mutation_matrix * Q_codon_array[,,unique.aa[k]])
                }
                diag(Q_codon_array[,,unique.aa[k]]) = 0
                diag(Q_codon_array[,,unique.aa[k]]) = -rowSums(Q_codon_array[,,unique.aa[k]])
            }
        }
        if(pars.to.evolve=="Ne"){
            Ne.current <- Ne[des[i]]
            Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne.current, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
            for(k in 1:21){
                if(diploid == TRUE){
                    Q_codon_array[,,unique.aa[k]] = 2 * Ne.current * (codon_mutation_matrix * Q_codon_array[,,unique.aa[k]])
                }else{
                    Q_codon_array[,,unique.aa[k]] = Ne.current * (codon_mutation_matrix * Q_codon_array[,,unique.aa[k]])
                }
                diag(Q_codon_array[,,unique.aa[k]]) = 0
                diag(Q_codon_array[,,unique.aa[k]]) = -rowSums(Q_codon_array[,,unique.aa[k]])
            }
        }
        Q_codon <- Q_codon_array[,,aa.optim]
        p <- expm(Q_codon * edge.length[i], method="Ward77")[sim.codon.data.site[anc[i]], ]
        sim.codon.data.site[des[i]] <- sample.int(dim(Q_codon)[2], size = 1, FALSE, prob = p)
    }
    sim.codon.data.site <- sim.codon.data.site[1:ntips]
    return(sim.codon.data.site)
}


BrownianEvolveParameters <- function(phy, start.value, rate, logspace=TRUE){
    phy <- reorder(phy, "postorder")
    ntips <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- ntips + 1
    evolved.parameter <- integer(ntips + phy$Nnode)
    evolved.parameter[ROOT] <- start.value
    anc <- phy$edge[,1]
    des <- phy$edge[,2]
    edge.length <- phy$edge.length
    for(i in N:1) {
        evolved.parameter[des[i]] <- rnorm(1, mean=evolved.parameter[anc[i]], sd=sqrt(edge.length[i] * rate))
    }
    if(logspace==TRUE){
        return(exp(evolved.parameter))
    }else{
        return(evolved.parameter)
    }
}


OUEvolveParameters <- function(phy, alpha, sigma.sq, mean, logspace=TRUE){
    phy <- reorder(phy, "postorder")
    ntips <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- ntips + 1
    evolved.parameter <- integer(ntips + phy$Nnode)
    evolved.parameter[ROOT] <- mean
    anc <- phy$edge[,1]
    des <- phy$edge[,2]
    edge.length <- phy$edge.length
    for(i in N:1) {
        evolved.parameter[des[i]] = evolved.parameter[anc[i]] * exp(-alpha*edge.length[i])+(mean)*(1-exp(-alpha*edge.length[i]))+sigma.sq*rnorm(1,0,1)*sqrt((1-exp(-2*alpha*edge.length[i]))/(2*alpha))
    }
    if(logspace==TRUE){
        return(exp(evolved.parameter))
    }else{
        return(evolved.parameter)
    }
}


#' @title Simulate DNA under the SELAC model
#'
#' @description
#' Simulates nucleotide data based on parameters under the SELAC model
#'
#' @param phy The phylogenetic tree with branch lengths.
#' @param pars A vector of parameters used for the simulation. They are ordered as follows: C.q.phi, alpha, beta, Ne, base.freqs for A C G, and the rates for the nucleotide model.
#' @param aa.optim_array A vector of optimal amino acids for each site to be simulated.
#' @param root.codon.frequencies A vector of codon frequencies for each possible optimal amino acid. Thus, the vector is of length 64x21.
#' @param root.codon.array A matrix of codon frequencies for each possible optimal amino acid. Rows are aa (including stop codon), cols are codons.
#' @param numcode The ncbi genetic code number for translation. By default the standard (numcode=1) genetic code is used.
#' @param aa.properties User-supplied amino acid distance properties. By default we assume Grantham (1974) properties.
#' @param nuc.model Indicates what type nucleotide model to use. There are three options: "JC", "GTR", or "UNREST".
#' @param k.levels Provides how many levels in the polynomial. By default we assume a single level (i.e., linear).
#' @param diploid A logical indicating whether or not the organism is diploid or not.
#'
#' @details
#' Simulates a nucleotide matrix using parameters under the SELAC model. Note that the output can be written to a fasta file using the write.dna() function in the \code{ape} package.
SelacSimulator <- function(phy, pars, aa.optim_array, root.codon.frequencies=NULL, root.codon.array=NULL, numcode=1, aa.properties=NULL, nuc.model, k.levels=0, diploid=TRUE){
    nsites <- length(aa.optim_array)
    #Start organizing the user input parameters:
    C.q.phi <- pars[1]
    C=4
    q=4e-7
    Phi.q <- C.q.phi / C
    Phi <- Phi.q / q
    alpha <- pars[2]
    beta <- pars[3]
    gamma <- GetAADistanceStartingParameters(aa.properties)[3]
    Ne <- pars[4]
    
    if(k.levels > 0){
        if(nuc.model == "JC") {
            base.freqs=c(pars[5:7], 1-sum(pars[5:7]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "GTR") {
            base.freqs=c(pars[5:7], 1-sum(pars[5:7]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[10:length(pars)], model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "UNREST") {
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[10:length(pars)], model=nuc.model)
        }
    }else{
        if(nuc.model == "JC") {
            base.freqs=c(pars[5:7], 1-sum(pars[5:7]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "GTR") {
            base.freqs=c(pars[5:7], 1-sum(pars[5:7]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[8:length(pars)], model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "UNREST") {
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[8:length(pars)], model=nuc.model)
        }
    }
    #Generate our codon matrix:
    codon.index.matrix = CreateCodonMutationMatrixIndex()
    codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
    codon_mutation_matrix[is.na(codon_mutation_matrix)]=0
    #Generate our fixation probability array:
    unique.aa <- GetMatrixAANames(numcode)
    aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE)
    Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
    for(k in 1:21){
        if(diploid == TRUE){
            Q_codon_array[,,unique.aa[k]] = 2 * Ne * (codon_mutation_matrix * Q_codon_array[,,unique.aa[k]])
        }else{
            Q_codon_array[,,unique.aa[k]] = Ne * (codon_mutation_matrix * Q_codon_array[,,unique.aa[k]])
        }
        diag(Q_codon_array[,,unique.aa[k]]) = 0
        diag(Q_codon_array[,,unique.aa[k]]) = -rowSums(Q_codon_array[,,unique.aa[k]])
    }
    root.p_array <- NA
    if(is.null(root.codon.array)) {
	    #Generate matrix of root frequencies for each optimal AA:
	    root.p_array <- matrix(root.codon.frequencies, nrow=dim(Q_codon_array)[2], ncol=21) #make sure user gives you root.codon.frequencies in the right order
	    root.p_array <- t(root.p_array)
	    root.p_array <- root.p_array / rowSums(root.p_array)
	    rownames(root.p_array) <- unique.aa
    } else {
    	root.p_array <- root.codon.array
    }
	root.p_array[is.na(root.p_array)] = 0

    #Perform simulation by looping over desired number of sites. The optimal aa for any given site is based on the user input vector of optimal AA:
    sim.codon.data <- matrix(0, nrow=Ntip(phy), ncol=nsites)
    for(site in 1:nsites){
        Q_codon = Q_codon_array[,,aa.optim_array[site]]
        sim.codon.data[,site] = SingleSiteUpPass(phy, Q_codon=Q_codon, root.value=root.p_array[aa.optim_array[site],])
    }
    codon.names <- rownames(Q_codon_array[,,aa.optim_array[site]])
    #Finally, translate this information into a matrix of nucleotides -- this format allows for write.dna() to write a fasta formatted file:
    nucleotide.data <- c()
    for(codon.sequence in seq(1, nsites, by=1)){
        nucleotide.data <- cbind(nucleotide.data, t(matrix(unlist(strsplit(codon.names[sim.codon.data[,codon.sequence]], split="")), 3,length(phy$tip.label))))
    }
    rownames(nucleotide.data) <- phy$tip.label
    
    #Done.
    return(nucleotide.data)
}


#' @title Simulate DNA under the SELON model
#'
#' @description
#' Simulates nucleotide data based on parameters under the SELAC model
#'
#' @param phy The phylogenetic tree with branch lengths.
#' @param pars A vector of parameters used for the simulation. They are ordered as follows: a0, a1, a2, Ne, base.freqs for A C G, and the nucleotide rates.
#' @param nuc.optim_array A vector of optimal nucleotide for each site to be simulated.
#' @param nuc.model Indicates what type nucleotide model to use. There are three options: "JC", "GTR", or "UNREST".
#' @param diploid A logical indicating whether or not the organism is diploid or not.
#'
#' @details
#' Simulates a nucleotide matrix using parameters under the SELON model. Note that the output can be written to a fasta file using the write.dna() function in the \code{ape} package.
SelonSimulator <- function(phy, pars, nuc.optim_array, nuc.model, diploid=TRUE){
    nsites <- length(nuc.optim_array)
    
    #Start organizing the user input parameters:
    scalor <- pars[1]
    left.slope <- pars[2]
    right.slope <- pars[3]
    mid.point <- pars[4]
    Ne = 5e6
	scalor = scalor/Ne
    site.index <- 1:nsites
    position.multiplier.vector <- scalor * PositionSensitivityMultiplierSigmoid(left.slope, right.slope, mid.point, nsites)

    if(nuc.model == "JC") {
        base.freqs=c(pars[5:7], 1-sum(pars[5:7]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "GTR") {
        base.freqs=c(pars[5:7], 1-sum(pars[5:7]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[7:length(pars)], model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "UNREST") {
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[7:length(pars)], model=nuc.model)
    }
    if(diploid == TRUE){
        ploidy = 2
    }else{
        ploidy = 1
    }
    #Perform simulation by looping over desired number of sites. The optimal aa for any given site is based on the user input vector of optimal AA:
    sim.nuc.data <- matrix(0, nrow=Ntip(phy), ncol=nsites)
    for(site in 1:nsites){
        weight.matrix <- GetNucleotideFixationMatrix(1, position.multiplier=position.multiplier.vector[site], optimal.nucleotide=nuc.optim_array[site], Ne=Ne, diploid=diploid)
        Q_position <- (ploidy * Ne) * nuc.mutation.rates * weight.matrix
        #Rescaling Q matrix in order to have a 1 nucleotide change per site if the branch length was 1:
        diag(Q_position) = 0
        diag(Q_position) <- -rowSums(Q_position)
        sim.nuc.data[,site] = SingleSiteUpPass(phy, Q_codon=Q_position, root.value=base.freqs)
    }
    nuc.names <- n2s(0:3)
    #Finally, translate this information into a matrix of nucleotides -- this format allows for write.dna() to write a fasta formatted file:
    nucleotide.data <- c()
    for(nuc.sequence in seq(1, nsites, by=1)){
        nucleotide.data <- cbind(nucleotide.data, nuc.names[sim.nuc.data[,nuc.sequence]])
    }
    rownames(nucleotide.data) <- phy$tip.label
    
    #Done.
    return(nucleotide.data)
}


#' @title Simulate DNA under General-Time Reversible model
#'
#' @description
#' Simulates nucleotide data based on parameters under the GTR+G model
#'
#' @param phy The phylogenetic tree with branch lengths.
#' @param pars A vector of parameters used for the simulation. They are ordered as follows: gamma shape and the rates for the nucleotide model.
#' @param nsites The number of sites to simulate.
#' @param nuc.model Indicates what type nucleotide model to use. There are three options: "JC", "GTR", or "UNREST".
#' @param base.freqs The base frequencies for A C G T (in that order).
#' @param ncats The number of discrete gamma categories.
#'
#' @details
#' Simulates a nucleotide matrix using parameters under the GTR+G model. Note that the output can be written to a fasta file using the write.dna() function in the \code{ape} package.
NucSimulator <- function(phy, pars, nsites, nuc.model, base.freqs, ncats){
    
    if(nuc.model == "JC") {
        base.freqs=base.freqs
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "GTR") {
        base.freqs=base.freqs
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[2:6], model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "UNREST") {
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars, model=nuc.model)
    }
    
    Q_mat <- nuc.mutation.rates
    diag(Q_mat) = 0
    diag(Q_mat) <- -rowSums(Q_mat)
    
    rate.vector <- DiscreteGamma(pars[1], ncats)
    rate.indicator <- sample.int(dim(Q_mat)[2], nsites, TRUE, prob=rep(1/ncats, ncats))
    # Perform simulation by looping over desired number of sites. The optimal aa for any given site is based on the user input vector of optimal AA:
    sim.nuc.data <- matrix(0, nrow=Ntip(phy), ncol=nsites)
    for(site in 1:nsites){
        Q_tmp <- Q_mat * rate.vector[rate.indicator[site]]
        sim.nuc.data[,site] = SingleSiteUpPass(phy, Q_codon=Q_tmp, root.value=base.freqs)
    }
    nuc.names <- n2s(0:3)
    # Finally, translate this information into a matrix of nucleotides -- this format allows for write.dna() to write a fasta formatted file:
    nucleotide.data <- c()
    for(nuc.sequence in seq(1, nsites, by=1)){
        nucleotide.data <- cbind(nucleotide.data, nuc.names[sim.nuc.data[,nuc.sequence]])
    }
    rownames(nucleotide.data) <- phy$tip.label
    # Done.
    return(nucleotide.data)
}


#' @title Simulate DNA under the SELAC model and evolving rates
#'
#' @description
#' Simulates nucleotide data based on parameters under the SELAC model but assumes either Phi or Ne evolves along the tree.
#'
#' @param phy The phylogenetic tree with branch lengths.
#' @param pars A vector of parameters used for the simulation. They are ordered as follows: C.q.phi, alpha, beta, and Ne.
#' @param aa.optim_array A vector of optimal amino acids for each site to be simulated.
#' @param root.codon.frequencies A vector of codon frequencies for each possible optimal amino acid. Thus, the vector is of length 64x64.
#' @param numcode The The ncbi genetic code number for translation. By default the standard (numcode=1) genetic code is used.
#' @param aa.properties User-supplied amino acid distance properties. By default we assume Grantham (1974) properties.
#' @param nuc.model Indicates what type nucleotide model to use. There are three options: "JC", "GTR", or "UNREST".
#' @param k.levels Provides how many levels in the polynomial. By default we assume a single level (i.e., linear).
#' @param diploid A logical indicating whether or not the organism is diploid or not.
#' @param pars.to.evolve Indicates which parameters to assume evolve along the tree. Only two options: "phi" or "Ne".
#' @param evolve.type The process by which the focal parameter evovles. There are two options: Brownian motion ("BM") or Ornstein-Uhlenbeck ("OU").
#' @param evolve.pars The process parameters used to simulate focal parameter evolution. Under "BM", the order is root.state, rate; under "OU", the order is alpha, sigma.sq, and the mean.
#' @param Ne.vals.evolved Under selac we assume a global Ne for all genes. Thus, when the focal parameter to evolve is "Ne", then a user specified vector of simulated Ne values are provided here.
#'
#' @details
#' Simulates a nucleotide matrix using parameters under the SELAC model, but allows either Phi or Ne to evolve along the tree. Note that the output can be written to a fasta file using the write.dna() function in the \code{ape} package.
SelacSimulatorEvolvingRates <- function(phy, pars, aa.optim_array, root.codon.frequencies, numcode=1, aa.properties=NULL, nuc.model, k.levels=0, diploid=TRUE, pars.to.evolve="phi", evolve.type="BM", evolve.pars=c(1,0), Ne.vals.evolved=NULL){
    nsites <- length(aa.optim_array)
    # Start organizing the user input parameters:
    C.q.phi <- pars[1]
    C=4
    q=4e-7
    Phi.q <- C.q.phi / C
    Phi.values <- Phi.q / q
    alpha <- pars[2]
    beta <- pars[3]
    gamma <- GetAADistanceStartingParameters(aa.properties)[3]
    Ne.values <- pars[4]
    if(pars.to.evolve == "phi"){
        if(evolve.type == "BM"){
            Phi.values <- BrownianEvolveParameters(phy=phy, start.value = evolve.pars[1], rate=evolve.pars[2])
            print(Phi.values)
        }
        if(evolve.type == "OU"){
            Phi.values <- OUEvolveParameters(phy=phy, alpha=evolve.pars[1], sigma.sq=evolve.pars[2], mean=evolve.pars[3])
            print(Phi.values)
        }
    }
    if(pars.to.evolve == "Ne"){
        Ne.values = Ne.vals.evolved
        print(Ne.values)
        print("Ne")
    }
    if(k.levels > 0){
        if(nuc.model == "JC") {
            base.freqs=c(pars[5:7], 1-sum(pars[5:7]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "GTR") {
            base.freqs=c(pars[5:7], 1-sum(pars[5:7]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[10:length(pars)], model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "UNREST") {
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[10:length(pars)], model=nuc.model)
        }
    }else{
        if(nuc.model == "JC") {
            base.freqs=c(pars[5:7], 1-sum(pars[5:7]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "GTR") {
            base.freqs=c(pars[5:7], 1-sum(pars[5:7]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[8:length(pars)], model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "UNREST") {
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[8:length(pars)], model=nuc.model)
        }
    }
    #Generate our codon matrix:
    codon.index.matrix = CreateCodonMutationMatrixIndex()
    codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
    codon_mutation_matrix[is.na(codon_mutation_matrix)]=0
    #Generate our fixation probability array:
    unique.aa <- GetMatrixAANames(numcode)
    aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE)
    #Generate matrix of root frequencies for each optimal AA:
    root.p_array <- matrix(root.codon.frequencies, nrow=dim(codon.index.matrix)[2], ncol=21) #make sure user gives you root.codon.frequencies in the right order
    root.p_array <- t(root.p_array)
    root.p_array <- root.p_array / rowSums(root.p_array)
    rownames(root.p_array) <- unique.aa
    root.p_array[is.na(root.p_array)] = 0
    #Perform simulation by looping over desired number of sites. The optimal aa for any given site is based on the user input vector of optimal AA:
    sim.codon.data <- matrix(0, nrow=Ntip(phy), ncol=nsites)
    for(site in 1:nsites){
        sim.codon.data[,site] = SingleSiteUpPassEvolvingRates(phy, root.value=root.p_array[aa.optim_array[site],], aa.distances=aa.distances, codon_mutation_matrix=codon_mutation_matrix, aa.optim=aa.optim_array[site], nsites=nsites, C=C, Phi=Phi.values, q=q, Ne=Ne.values, numcode=numcode, diploid=diploid, pars.to.evolve=pars.to.evolve)
    }
    codon.names <- rownames(codon.index.matrix)
    #Finally, translate this information into a matrix of nucleotides -- this format allows for write.dna() to write a fasta formatted file:
    nucleotide.data <- c()
    for(codon.sequence in seq(1, nsites, by=1)){
        nucleotide.data <- cbind(nucleotide.data, t(matrix(unlist(strsplit(codon.names[sim.codon.data[,codon.sequence]], split="")), 3,length(phy$tip.label))))
    }
    rownames(nucleotide.data) <- phy$tip.label
    
    #Done.
    return(nucleotide.data)
}


# Simulate using an individual based model
#
# @param phy A phylogram
# @param root.states The starting sequence
# @param phi.mean The mean value of phi across sites
# @param phi.sd The standard deviation of phi across sites
# @param Ne Effective population size. Must be an integer.
# @param GTR.rates rates in the usual order for selac
# @param recombination.rate The rate at which there is recombination
# @param optimal.aa.switch.rate The rate at which the optimal amino acid changes
# @param diploid If TRUE, simulate using a diploid population
# @param aa.properties The aa.properties
# @param population Rather than a single sequence to start the sim, allows use of a population
# @param brlen.conversion How many generations correspond to one unit of branch length on a tree
# @param burnin.generations How many generations to simulate on before actually starting the sim.
#FullWrightFisherSimulator <- function(phy, root.states=NULL, phi.mean, phi.sd=0, Ne=100, GTR.rates, recombination.rate, optimal.aa.switch.rate, diploid=TRUE, aa.properties=NULL, population=NULL, brlen.conversion=1e6, burnin.generations=1e3) {
#    if (is.null(root.states) & is.null(population)) {
#        stop("You must specify the root state or the root population")
#    }
#    if(is.null(population)) {
#        population <- replicate(n=Ne, root.states)
#    }
#}



###TESTING CODE:
#phy <- read.tree("../data/rokasYeast.tre")
#phy <- drop.tip(phy, "Calb")
#gene.tmp <- read.dna("../data/gene1Yeast.fasta", format="fasta")
#gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[1:7,])
#codon.data <- DNAbinToCodonNumeric(gene.tmp)
#codon.data <- codon.data[phy$tip.label,]
#aa.data <- ConvertCodonNumericDataToAAData(codon.data, numcode=1)
#aa.optim <- apply(aa.data[, -1], 2, GetMaxName) #starting values for all, final values for majrule
#root.codon.freq <- CodonEquilibriumFrequencies(codon.data[,-1], aa.optim, numcode=1)

#tmp <- SelacSimulatorEvolvingRates(phy=phy, pars=c(1e-9,.4,.1,5e6,.25,.25,.25), aa.optim_array=aa.optim, root.codon.frequencies=root.codon.freq, numcode=1, aa.properties=NULL, nuc.model="JC")





