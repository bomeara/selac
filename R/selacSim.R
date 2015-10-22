
######################################################################################################################################
######################################################################################################################################
### SELAC simulator -- simulator for SELection on Amino acids and/or Codons model
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu

source("selac.R")

SingleSiteUpPass <- function(phy, Q_codon, root.value){	
	#Randomly choose starting state at root using the root.values as the probability:
	root.value = sample.int(64, 1, FALSE, prob=root.value)
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
		sim.codon.data.site[des[i]] <- sample.int(64, size = 1, FALSE, prob = p)
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





