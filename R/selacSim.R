
######################################################################################################################################
######################################################################################################################################
### SELAC simulator -- simulator for SELection on Amino acids and/or Codons model
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu

source("selac.R")

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
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    edge.length <- phy$edge.length
	for (i in N:1) {
		p <- expm(Q_codon * edge.length[i], method="Ward77")[sim.codon.data.site[anc[i]], ]
		sim.codon.data.site[des[i]] <- sample.int(dim(Q_codon)[2], size = 1, FALSE, prob = p)
	}
	sim.codon.data.site <- sim.codon.data.site[1:ntips]
	return(sim.codon.data.site)
}


SelacSimulator <- function(phy, pars, aa.optim_array, root.codon.frequencies, numcode=1, aa.properties=NULL, nuc.model, k.levels=0, diploid=TRUE){
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
	codon_mutation_matrix = c(as.vector(nuc.mutation.rates), 0)[codon.index.matrix]
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
	#Generate matrix of root frequencies for each optimal AA:
	root.p_array <- matrix(root.codon.frequencies, nrow=dim(Q_codon_array)[2], ncol=21) #make sure user gives you root.codon.frequencies in the right order
	root.p_array <- t(root.p_array)
	root.p_array <- root.p_array / rowSums(root.p_array)
	rownames(root.p_array) <- unique.aa
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


SelonSimulator <- function(phy, pars, nuc.optim_array, nuc.model, diploid=TRUE){

	nsites <- length(nuc.optim_array)
	#Start organizing the user input parameters:
	a0 <- pars[1]
	a1 <- pars[2]
	a2 <- pars[3]
	Ne <- pars[4]
	
	site.mid.point <- ceiling(nsites/2)
    site.index <- 1:nsites
    site.index <- site.index - site.mid.point
	position.multiplier.vector <- PositionSensitivityMultiplier(a0, a1-site.mid.point, a2, site.index)

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
	if(is.null(root.p_array)) {
		#Generate matrix of equal frequencies for each site:
        root.p_array <- rep(0.25, 4)
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
		nucleotide.data <- cbind(nucleotide.data, t(matrix(nuc.names[sim.nu.data[,nuc.sequence]], split="")), 3,length(phy$tip.label))))
	}
	rownames(nucleotide.data) <- phy$tip.label
	
	#Done.
	return(nucleotide.data)
}

#' Simulate using an individual based model
#'
#' @param phy A phylogram
#' @param root.states The starting sequence
#' @param phi.mean The mean value of phi across sites
#' @param phi.sd The standard deviation of phi across sites
#' @param Ne Effective population size. Must be an integer.
#' @param GTR rates in the usual order for selac
#' @param recombination.rate The rate at which there is recombination
#' @param optimal.aa.switch.rate The rate at which the optimal amino acid changes
#' @param diplod If TRUE, simulate using a diploid population
#' @param aa.properties The aa.properties
#' @param population Rather than a single sequence to start the sim, allows use of a population
#' @param brlen.conversion How many generations correspond to one unit of branch length on a tree
#' @param burnin.generations How many generations to simulate on before actually starting the sim.
FullWrightFisherSimulator <- function(phy, root.states=NULL, phi.mean, phi.sd=0, Ne=100, GTR.rates, recombination.rate, optimal.aa.switch.rate, diploid=TRUE, aa.properties=NULL, population=NULL, brlen.conversion=1e6, burnin.generations=1e3) {
	if (is.null(root.states) & is.null(population)) {
		stop("You must specify the root state or the root population")
	}
	if(is.null(population)) {
		population <- replicate(n=Ne, root.states)
	}
		
	
}

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

#tmp <- SelacSimulator(nsites=10, phy=phy, pars=c(1e-9,.4,.1,5e6,.25,.25,.25), aa.distances, aa.optim_array=aa.optim, root.codon.frequencies=root.codon.freq, numcode=1, aa.properties=NULL, nuc.model="JC")





