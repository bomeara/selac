
######################################################################################################################################
######################################################################################################################################
### Tests for model adequacy based on functionality
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu


GetTipIntervalStateSingleSite <- function(charnum=1, codon.data, phy, root.p=NULL, taxon.to.drop=1, Q.to.reconstruct, Q.to.simulate,  model.to.reconstruct.under="selac", model.to.simulate.under="selac"){

    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    nl <- nrow(Q.to.reconstruct)
    #Now we need to build the matrix of likelihoods to pass to dev.raydisc:
    liks <- matrix(0, nb.tip + nb.node, nl)
    #Now loop through the tips.
    for(i in 1:nb.tip){
        #The codon at a site for a species is not NA, then just put a 1 in the appropriate column.
        #Note: We add charnum+1, because the first column in the data is the species labels:
        if(codon.data[i,charnum+1] < 65){
            liks[i,codon.data[i,charnum+1]] <- 1
        }else{
            #If here, then the site has no data, so we treat it as ambiguous for all possible codons. Likely things might be more complicated, but this can be modified later:
            liks[i,] <- 1
            #This is to deal with stop codons, which are effectively removed from the model at this point. Someday they might be allowed back in. If so, the following line needs to be dealt with.
        }
    }

    TIPS <- 1:nb.tip
    anc <- unique(phy$edge[,1])

    diag(Q.to.reconstruct) <- 0
    diag(Q.to.reconstruct) <- -rowSums(Q.to.reconstruct)

    #A temporary likelihood matrix so that the original does not get written over:
    liks.down <- liks
    #A transpose of Q for assessing probability of j to i, rather than i to j:
    tranQ <- t(Q.to.reconstruct)
    comp <- matrix(0, nb.tip+nb.node, ncol(liks))
    #The first down-pass: The same algorithm as in the main function to calculate the conditional likelihood at each node:
    for (i  in seq(from = 1, length.out = nb.node)) {
        #the ancestral node at row i is called focal
        focal <- anc[i]
        #Get descendant information of focal
        desRows<-which(phy$edge[,1]==focal)
        desNodes<-phy$edge[desRows,2]
        v <- 1
        for (desIndex in sequence(length(desRows))){
            if(!desNodes[desIndex] == taxon.to.drop){
                v <- v*expm(Q.to.reconstruct * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks.down[desNodes[desIndex],]
            }
        }
        comp[focal] <- sum(v)
        liks.down[focal, ] <- v/comp[focal]
    }
    root <- nb.tip + 1L
    #Enter the root defined root probabilities if they are supplied by the user:
    equil.root <- NULL
    for(i in 1:ncol(Q.to.reconstruct)){
        posrows <- which(Q.to.reconstruct[,i] >= 0)
        rowsum <- sum(Q.to.reconstruct[posrows,i])
        poscols <- which(Q.to.reconstruct[i,] >= 0)
        colsum <- sum(Q.to.reconstruct[i,poscols])
        equil.root <- c(equil.root,rowsum/(rowsum+colsum))
    }
    if (is.null(root.p)){
        flat.root = equil.root
        k.rates <- 1/length(which(!is.na(equil.root)))
        flat.root[!is.na(flat.root)] = k.rates
        flat.root[is.na(flat.root)] = 0
        liks.down[root, ] <- flat.root * liks.down[root, ]
        liks.down[root, ] <- liks.down[root,] / sum(liks.down[root, ])
        root.p = flat.root
    }

    #The up-pass
    liks.up <- liks
    states <- apply(liks,1,which.max)
    N <- dim(phy$edge)[1]
    comp <- numeric(nb.tip + nb.node)
    for(i in length(anc):1){
        focal <- anc[i]
        if(!focal==root){
            #Gets mother and sister information of focal:
            focalRow <- which(phy$edge[,2]==focal)
            motherRow <- which(phy$edge[,1]==phy$edge[focalRow,1])
            motherNode <- phy$edge[focalRow,1]
            desNodes <- phy$edge[motherRow,2]
            sisterNodes <- desNodes[(which(!desNodes==focal))]
            sisterRows <- which(phy$edge[,2]%in%sisterNodes==TRUE)
            #If the mother is not the root then you are calculating the probability of being in either state.
            #But note we are assessing the reverse transition, j to i, rather than i to j, so we transpose Q to carry out this calculation:
            if(motherNode!=root){
                v <- expm(tranQ * phy$edge.length[which(phy$edge[,2]==motherNode)], method=c("Ward77")) %*% liks.up[motherNode,]
            }
            #If the mother is the root then just use the marginal. This can also be the prior, which I think is the equilibrium frequency.
            #But for now we are just going to use the marginal at the root -- it is unclear what Mesquite does.
            else{
                v <- root.p
            }
            #Now calculate the probability that each sister is in either state. Sister can be more than 1 when the node is a polytomy.
            #This is essentially calculating the product of the mothers probability and the sisters probability:
            for (sisterIndex in sequence(length(sisterRows))){
                if(!sisterNodes[sisterIndex] == taxon.to.drop){
                    v <- v*expm(Q.to.reconstruct * phy$edge.length[sisterRows[sisterIndex]], method=c("Ward77")) %*% liks.down[sisterNodes[sisterIndex],]
                }
            }
            comp[focal] <- sum(v)
            liks.up[focal,] <- v/comp[focal]
        }
    }
    #The final pass
    liks.final <- liks
    comp <- numeric(nb.tip + nb.node)

    #In this final pass, root is never encountered. But its OK, because root likelihoods are set after the loop:
    for (i in seq(from = 1, length.out = nb.node-1)) {
        #the ancestral node at row i is called focal
        focal <- anc[i]
        focalRows <- which(phy$edge[,2]==focal)
        #Now you are assessing the change along the branch subtending the focal by multiplying the probability of
        #everything at and above focal by the probability of the mother and all the sisters given time t:
        v <- liks.down[focal,]*expm(tranQ * phy$edge.length[focalRows], method=c("Ward77")) %*% liks.up[focal,]
        comp[focal] <- sum(v)
        liks.final[focal, ] <- v/comp[focal]
    }

    if(model.to.reconstruct.under == "gtr" & model.to.simulate.under == "selac"){
        return(liks.final)
    }else{
        tot.interval <- phy$edge.length[phy$edge[,2]==taxon.to.drop]
        prop.interval <- seq(0,1 , by=0.05)
        time.interval <- tot.interval * prop.interval
        focal <- taxon.to.drop
        focalRows <- which(phy$edge[,2]==focal)
        starting.probs <- liks.final[phy$edge[focalRows,1],]
        focal.starting.state <- sample(1:dim(Q.to.reconstruct)[2], 1, prob=starting.probs)

        if(model.to.reconstruct.under == "selac"){
            if(dim(Q.to.simulate)[2]==4){
                focal.starting.state.converted <- as.vector(s2n(strsplit(.codon.name[focal.starting.state], "")[[1]]))+1
                reconstructed.sequence <- c()
                for(site.index in 1:3){
                    site.by.site.codon.reconstruction <- c()
                    focal.starting.state <- focal.starting.state.converted[site.index]
                    for (time.index in 2:length(time.interval)){
                        p <- expm(Q.to.simulate * time.interval[time.index], method="Ward77")[focal.starting.state, ]
                        focal.starting.state <- sample.int(dim(Q.to.simulate)[2], size = 1, FALSE, prob = p)
                        site.by.site.codon.reconstruction <- c(site.by.site.codon.reconstruction, focal.starting.state)
                    }
                    reconstructed.sequence <- cbind(reconstructed.sequence, site.by.site.codon.reconstruction)
                }
                reconstructed.sequence <- rbind(focal.starting.state.converted, reconstructed.sequence)
                reconstructed.sequence.codon <- c()
                for(time.index in 1:length(time.interval)){
                    reconstructed.sequence.codon <- c(reconstructed.sequence.codon, which(.codon.name==paste(as.vector(n2s(reconstructed.sequence[time.index,]-1)), collapse="")))
                }
            }else{
                reconstructed.sequence.codon <- c(focal.starting.state)
                for (time.index in 2:length(time.interval)){
                    p <- expm(Q.to.simulate * time.interval[time.index], method="Ward77")[focal.starting.state, ]
                    focal.starting.state <- sample.int(dim(Q.to.simulate)[2], size = 1, FALSE, prob = p)
                    reconstructed.sequence.codon <- c(reconstructed.sequence.codon, focal.starting.state)
                }
            }
        }else{
            reconstructed.sequence.codon <- c(focal.starting.state)
            for (time.index in 2:length(time.interval)){
                p <- expm(Q.to.simulate * time.interval[time.index], method="Ward77")[focal.starting.state, ]
                focal.starting.state <- sample.int(dim(Q.to.simulate)[2], size = 1, FALSE, prob = p)
                reconstructed.sequence.codon <- c(reconstructed.sequence.codon, focal.starting.state)
            }

        }
        return(reconstructed.sequence.codon)
    }
}



GetKnownFunctionality <- function(selac.obj, partition.number, fasta.rows.to.keep=NULL, taxon.to.do, numcode=1){
    partitions <- selac.obj$partitions
    gene.tmp <- read.dna(partitions[partition.number], format='fasta')
    if(!is.null(fasta.rows.to.keep)){
        gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
    }else{
        gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
    }
    codon.data <- DNAbinToCodonNumeric(gene.tmp)
    codon.data <- codon.data[selac.obj$phy$tip.label,]

    aa.optim_array <- selac.obj$aa.optim[[partition.number]]
    aa.data <- ConvertCodonNumericDataToAAData(codon.data, numcode=numcode)

    actual.sequence <- as.vector(aa.data[taxon.to.do,-1])
    functionality.taxon <- GetFunctionality(gene.length=selac.obj$nsites[partition.number], aa.data=actual.sequence, optimal.aa=selac.obj$aa.optim[[partition.number]], alpha=selac.obj$mle.pars[partition.number,2], beta=selac.obj$mle.pars[partition.number,3], gamma=selac.obj$volume.fixed.value, aa.properties=selac.obj$aa.properties)

    return(functionality.taxon)
}


GetFunctionalityModelAdequacy <- function(gene.length, aa.data, optimal.aa, alpha, beta, gamma, gp=NULL, aa.properties=NULL, fmutsel=FALSE){
    if(is.null(aa.properties)) {
        #     aa.properties <- structure(c(0, 2.75, 1.38, 0.92, 0, 0.74, 0.58, 0, 0.33, 0, 0,
        # 1.33, 0.39, 0.89, 0.65, 1.42, 0.71, 0, 0.13, 0.2, 8.1, 5.5, 13,
        # 12.3, 5.2, 9, 10.4, 5.2, 11.3, 4.9, 5.7, 11.6, 8, 10.5, 10.5,
        # 9.2, 8.6, 5.9, 5.4, 6.2, 31, 55, 54, 83, 132, 3, 96, 111, 119,
        # 111, 105, 56, 32.5, 85, 124, 32, 61, 84, 170, 136), .Dim = c(20L,
        # 3L), .Dimnames = list(c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly",
        # "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg",
        # "Ser", "Thr", "Val", "Trp", "Tyr"), c("c", "p", "v"))) #properties from Grantham paper
        aa.properties <- structure(c(0, 2.75, 1.38, 0.92, 0, 0.74, 0.58, 0, 0.33, 0, 0,
        1.33, 0.39, 0.89, 0.65, 1.42, 0.71, 0, 0.13, 0.2, 8.1, 5.5, 13,
        12.3, 5.2, 9, 10.4, 5.2, 11.3, 4.9, 5.7, 11.6, 8, 10.5, 10.5,
        9.2, 8.6, 5.9, 5.4, 6.2, 31, 55, 54, 83, 132, 3, 96, 111, 119,
        111, 105, 56, 32.5, 85, 124, 32, 61, 84, 170, 136), .Dim = c(20L,
        3L), .Dimnames = list(c("A", "C", "D", "E", "F", "G",
        "H", "I", "K", "L", "M", "N", "P", "Q", "R",
        "S", "T", "V", "W", "Y"), c("c", "p", "v"))) #properties from Grantham paper
    }
    if(is.null(gp)){
        gp <- rep(1, gene.length)
    }
    aa.distances <- c()
    #Note using only the second row, because we are comparing empirical S. cervisae rates:
    for(site.index in 1:gene.length){
        if(aa.data[site.index]!="NA" & aa.data[site.index]!="*"){
            #broke this up to make debugging easier:
            distance <- ((alpha*(aa.properties[aa.data[site.index],1] - aa.properties[optimal.aa[site.index],1])^2 + beta*(aa.properties[aa.data[site.index],2]-aa.properties[optimal.aa[site.index],2])^2+gamma*(aa.properties[aa.data[site.index],3]-aa.properties[optimal.aa[site.index],3])^2)^(1/2))
            aa.distances <- c(aa.distances, (1+gp[site.index]*distance))
        }else{
            aa.distances <- c(aa.distances, 0)
            gene.length <- gene.length - 1
        }
    }
    functionality = 1/((1/gene.length) * sum(aa.distances))
    return(functionality)
}


GetGtrSimulateInfo <- function(selac.obj, partition.number){

    pars <- c(selac.obj$mle.pars[partition.number,])
    phy <- selac.obj$phy
    partitions <- selac.obj$partitions
    include.gamma <- selac.obj$include.gamma
    aa.properties <- selac.obj$aa.properties
    diploid <- selac.obj$diploid
    gamma.type <- selac.obj$gamma.type
    ncats <- selac.obj$ncats
    numcode <- selac.obj$numcode
    gamma <- selac.obj$volume.fixed.value
    nuc.model <- selac.obj$nuc.model
    k.levels <- selac.obj$k.levels
    parallel.type <- "by.gene"
    n.cores <- NULL
    codon.freq.by.aa <- NULL
    codon.freq.by.gene <- selac.obj$empirical.base.freqs[[partition.number]]
    nsites <- selac.obj$nsites[partition.number]

    pars <- c(selac.obj$mle.pars[partition.number,])
    if(include.gamma == TRUE){
        shape = pars[1]
        pars = pars[-1]
    }
    if(length(pars)==0){
        transition.rates <- 1
    }else{
        transition.rates <- pars[1:length(pars)]
    }

    nuc.mutation.rates <- CreateNucleotideMutationMatrix(transition.rates, model=nuc.model)

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

        #Rescaling Q matrix in order to have a 1 nucleotide change per site if the branch length was 1:
        if(is.null(codon.freq.by.gene)) {
            #Generate matrix of equal frequencies for each site:
            root.p_array <- rep(0.25, 4)
        }else{
            root.p_array <- codon.freq.by.gene
        }

        diag(nuc.mutation.rates) = 0
        nuc.mutation.rates = t(nuc.mutation.rates * root.p_array)
        diag(nuc.mutation.rates) = -rowSums(nuc.mutation.rates)
        scale.factor <- -sum(diag(nuc.mutation.rates) * root.p_array)

        Q_array <- nuc.mutation.rates * (1/scale.factor)
    }else{
        if(is.null(codon.freq.by.gene)) {
            #Generate matrix of equal frequencies for each site:
            root.p_array <- rep(0.25, 4)
        }else{
            root.p_array <- codon.freq.by.gene
        }

        #Rescaling Q matrix in order to have a 1 nucleotide change per site if the branch length was 1:
        diag(nuc.mutation.rates) = 0
        nuc.mutation.rates = t(nuc.mutation.rates * root.p_array)
        diag(nuc.mutation.rates) = -rowSums(nuc.mutation.rates)
        scale.factor <- -sum(diag(nuc.mutation.rates) * root.p_array)
        Q_array <- nuc.mutation.rates * (1/scale.factor)
    }
    obj <- NULL
    obj$gamma.rates <- rates.k
    obj$gamma.weights <- weights.k
    obj$Q_matrix <- Q_array
    return(obj)
}


GetFMutSelSimulateInfo <- function(selac.obj, partition.number){

    pars <- c(selac.obj$mle.pars[partition.number,])
    phy <- selac.obj$phy
    partitions <- selac.obj$partitions
    include.gamma <- selac.obj$include.gamma
    aa.properties <- selac.obj$aa.properties
    diploid <- selac.obj$diploid
    gamma.type <- selac.obj$gamma.type
    ncats <- selac.obj$ncats
    numcode <- selac.obj$numcode
    gamma <- selac.obj$volume.fixed.value
    nuc.model <- selac.obj$nuc.model
    k.levels <- selac.obj$k.levels
    nuc.model <- selac.obj$nuc.model
    parallel.type <- "by.gene"
    n.cores <- NULL
    codon.freq.by.aa <- NULL
    codon.freq.by.gene <- NULL
    nsites <- selac.obj$nsites[partition.number]
    pars <- c(selac.obj$mle.pars[partition.number,])
    if(nuc.model == "JC") {
        base.freqs <- c(pars[1:3], 1-sum(pars[1:3]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
        pars = pars[-(1:3)]
    }
    if(nuc.model == "GTR") {
        base.freqs <- c(pars[1:3], 1-sum(pars[1:3]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[4:8], model=nuc.model, base.freqs=base.freqs)
        pars = pars[-(1:8)]
    }
    if(nuc.model == "UNREST") {
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[1:11], model=nuc.model, base.freqs=NULL)
        pars = pars[-(1:11)]
    }
    #.codon.sets <- .codon.sets
    n.codons <- dim(.codon.sets)[1]
    codon.eq.freq <- numeric(n.codons)
    fitness.pars <- c(pars[-1],0)
    fitness.pars.ordered <- numeric(n.codons)
    if(length(fitness.pars)>21){
        fitness.pars.ordered = c(fitness.pars[1:48], 0, fitness.pars[49], 0, fitness.pars[50:54], 0, fitness.pars[55:61])
    }else{
        fitness.pars.ordered <- numeric(n.codons)
        #codon.set.translate <- apply(.codon.sets, 2, n2s)
        #codon.name <- apply(.codon.set.translate, 1, paste, collapse="")
        aa.translations <- .aa.translation[[numcode]][.codon.name]
        unique.aa.nostop = .unique.aa[-which(.unique.aa=="*")]
        for(par.index in 1:length(unique.aa.nostop)){
            fitness.pars.ordered[which(aa.translations == unique.aa.nostop[par.index])] <- fitness.pars[par.index]
        }
    }
    for(codon.index in 1:n.codons){
        #In the canonical model stop codons are ignored. We do the same here.
        if(codon.index == 49 | codon.index == 51 | codon.index == 57){
            codon.eq.freq[codon.index] = 0
        }else{
            codon.eq.freq[codon.index] <- base.freqs[unname(.codon.sets[codon.index,1])+1] * base.freqs[unname(.codon.sets[codon.index,2])+1] * base.freqs[unname(.codon.sets[codon.index,3])+1] * exp(fitness.pars.ordered[codon.index])
        }
    }
    codon.eq.freq <- codon.eq.freq[1:64]/sum(codon.eq.freq[1:64])
    Q_codon = CreateCodonMutationMatrixMutSel(omega.par=pars[1], fitness.pars=fitness.pars.ordered, nuc.mutation.rates=nuc.mutation.rates, numcode=numcode)

    #Rescaling Q matrix in order to have a 1 nucleotide change per site if the branch length was 1:
    diag(Q_codon) = 0
    diag(Q_codon) = -rowSums(Q_codon)
    scale.factor <- -sum(diag(Q_codon) * codon.eq.freq, na.rm=TRUE)
    Q_array = Q_codon * (1/scale.factor)

    obj <- NULL
    obj$gamma.rates <- NULL
    obj$gamma.weights <- NULL
    obj$Q_matrix <- Q_array
    return(obj)
}



GetSelacSimulateInfo <- function(selac.obj, partition.number){
    pars <- c(selac.obj$mle.pars[partition.number,])
    phy <- selac.obj$phy
    partitions <- selac.obj$partitions
    include.gamma <- selac.obj$include.gamma
    aa.properties <- selac.obj$aa.properties
    diploid <- selac.obj$diploid
    gamma.type <- selac.obj$gamma.type
    ncats <- selac.obj$ncats
    numcode <- selac.obj$numcode
    gamma <- selac.obj$volume.fixed.value
    nuc.model <- selac.obj$nuc.model
    k.levels <- selac.obj$k.levels
    parallel.type <- "by.gene"
    n.cores <- NULL
    codon.freq.by.aa <- selac.obj$codon.freq.by.aa[[partition.number]]
    codon.freq.by.gene <- selac.obj$codon.freq.by.gene[[partition.number]]
    nsites <- selac.obj$nsites[partition.number]

    codon.index.matrix <- CreateCodonMutationMatrixIndex()

    if(include.gamma == TRUE){
        shape = pars[length(pars)]
        pars = pars[-length(pars)]
    }

    C.Phi.q.Ne <- pars[1]
    C <- 4
    q <- 4e-7
    Ne <- 5e6
    Phi.q.Ne <- C.Phi.q.Ne / C
    Phi.Ne <- Phi.q.Ne / q
    Phi <- Phi.Ne / Ne
    alpha <- pars[2]
    beta <- pars[3]
    gamma <- 0.0003990333

    if(k.levels > 0){
        if(nuc.model == "JC") {
            base.freqs=c(pars[4:6], 1-sum(pars[4:6]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
            poly.params <- pars[7:8]
        }
        if(nuc.model == "GTR") {
            base.freqs=c(pars[4:6], 1-sum(pars[4:6]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[9:length(pars)], model=nuc.model, base.freqs=base.freqs)
            poly.params <- pars[7:8]
        }
        if(nuc.model == "UNREST") {
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[6:length(pars)], model=nuc.model)
            poly.params <- pars[4:5]
        }
    }else{
        if(nuc.model == "JC") {
            base.freqs=c(pars[4:6], 1-sum(pars[4:6]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "GTR") {
            base.freqs=c(pars[4:6], 1-sum(pars[4:6]))
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[7:length(pars)], model=nuc.model, base.freqs=base.freqs)
        }
        if(nuc.model == "UNREST") {
            nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[4:length(pars)], model=nuc.model)
        }
    }

    codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
    codon_mutation_matrix[is.na(codon_mutation_matrix)]=0

    #We rescale the codon matrix only:
    diag(codon_mutation_matrix) = 0
    diag(codon_mutation_matrix) = -rowSums(codon_mutation_matrix)
    scale.factor <- -sum(diag(codon_mutation_matrix) * codon.freq.by.gene, na.rm=TRUE)
    codon_mutation_matrix_scaled = codon_mutation_matrix * (1/scale.factor)

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
        rate.Q_codon.list <- as.list(ncats)
        for(cat.index in 1:ncats){
            if(k.levels > 0){
                aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
            }else{
                aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
            }
            rate.Q_codon.list[[cat.index]] <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[cat.index], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
        }
        for(cat.index in 1:ncats){
            Q_codon_array <- rate.Q_codon.list[[cat.index]]
            for(aa.index in 1:21){
                if(diploid == TRUE){
                    Q_codon_array[,,.unique.aa[aa.index]] <- 2 * Ne * (codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[aa.index]])
                }else{
                    Q_codon_array[,,.unique.aa[aa.index]] <- Ne * (codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[aa.index]])
                }
                diag(Q_codon_array[,,.unique.aa[aa.index]]) <- 0
                diag(Q_codon_array[,,.unique.aa[aa.index]]) <- -rowSums(Q_codon_array[,,.unique.aa[aa.index]])
                Q_codon_array[,,.unique.aa[aa.index]] <- Q_codon_array[,,.unique.aa[aa.index]]
            }
            rate.Q_codon.list[[cat.index]] <- Q_codon_array
        }

        Q_array <- rate.Q_codon.list

    }else{
        rates.k <- NULL
        weights.k <- NULL
        if(k.levels > 0){
            aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
        }else{
            aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
        }
        Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999)
        #Finish the Q_array codon mutation matrix multiplication here:
        for(k in 1:21){
            if(diploid == TRUE){
                Q_codon_array[,,.unique.aa[k]] = (2 * Ne) * codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[k]]
            }else{
                Q_codon_array[,,.unique.aa[k]] = Ne * codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[k]]
            }
            diag(Q_codon_array[,,.unique.aa[k]]) = 0
            diag(Q_codon_array[,,.unique.aa[k]]) = -rowSums(Q_codon_array[,,.unique.aa[k]])
        }
        Q_array <- Q_codon_array
    }
    obj <- NULL
    obj$aa.optim_array <- selac.obj$aa.optim[[partition.number]]
    obj$gamma.rates <- rates.k
    obj$gamma.weights <- weights.k
    obj$Q_matrix <- Q_array
    return(obj)
}



GetIntervalSequencesAllSites <- function(model.to.reconstruct.under, model.to.simulate.under, selac.obj1, selac.obj2, aa.optim.input, fasta.rows.to.keep, taxon.to.drop, partition.number){

    phy <- selac.obj1$phy
    partitions <- selac.obj1$partitions
    include.gamma <- selac.obj1$include.gamma
    aa.properties <- selac.obj1$aa.properties
    diploid <- selac.obj1$diploid
    gamma.type <- selac.obj1$gamma.type
    ncats <- selac.obj1$ncats
    numcode <- selac.obj1$numcode
    gamma <- selac.obj1$volume.fixed.value
    nuc.model <- selac.obj1$nuc.model
    k.levels <- selac.obj1$k.levels
    parallel.type <- "by.gene"
    n.cores <- NULL

    if(model.to.reconstruct.under == "selac"){
        
        for(partition.index in partition.number:partition.number){
            pars <- c(selac.obj1$mle.pars[partition.index,])
            gene.tmp <- read.dna(partitions[partition.index], format='fasta')
            if(!is.null(fasta.rows.to.keep)){
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
            }else{
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
            }
            codon.data <- DNAbinToCodonNumeric(gene.tmp)
            codon.data <- codon.data[phy$tip.label,]
            codon.freq.by.aa <- selac.obj1$codon.freq.by.aa[[partition.index]]
            codon.freq.by.gene <- selac.obj1$codon.freq.by.gene[[partition.index]]
            if(is.null(aa.optim.input)){
                aa.optim_array <- selac.obj1$aa.optim[[partition.index]]
            }
        }

        nsites <- dim(codon.data)[2]-1
        prop.interval <- seq(0, 1, by=0.05)

        codon.index.matrix <- CreateCodonMutationMatrixIndex()
        nsites <- dim(codon.data)[2]-1

        interval.recon_array <- c()
        simulated.site.phi.cat <- c()

        if(include.gamma == TRUE){
            shape = pars[length(pars)]
            pars = pars[-length(pars)]
        }

        C.Phi.q.Ne <- pars[1]
        C <- 4
        q <- 4e-7
        Ne <- 5e6
        Phi.q.Ne <- C.Phi.q.Ne / C
        Phi.Ne <- Phi.q.Ne / q
        Phi <- Phi.Ne / Ne
        alpha <- pars[2]
        beta <- pars[3]
        gamma <- 0.0003990333

        if(k.levels > 0){
            if(nuc.model == "JC") {
                base.freqs=c(pars[4:6], 1-sum(pars[4:6]))
                nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
                poly.params <- pars[7:8]
            }
            if(nuc.model == "GTR") {
                base.freqs=c(pars[4:6], 1-sum(pars[4:6]))
                nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[9:length(pars)], model=nuc.model, base.freqs=base.freqs)
                poly.params <- pars[7:8]
            }
            if(nuc.model == "UNREST") {
                nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[6:length(pars)], model=nuc.model)
                poly.params <- pars[4:5]
            }
        }else{
            if(nuc.model == "JC") {
                base.freqs=c(pars[4:6], 1-sum(pars[4:6]))
                nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
            }
            if(nuc.model == "GTR") {
                base.freqs=c(pars[4:6], 1-sum(pars[4:6]))
                nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[7:length(pars)], model=nuc.model, base.freqs=base.freqs)
            }
            if(nuc.model == "UNREST") {
                nuc.mutation.rates <- CreateNucleotideMutationMatrix(pars[4:length(pars)], model=nuc.model)
            }
        }

        codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
        codon_mutation_matrix[is.na(codon_mutation_matrix)]=0

        #We rescale the codon matrix only:
        diag(codon_mutation_matrix) = 0
        diag(codon_mutation_matrix) = -rowSums(codon_mutation_matrix)
        scale.factor <- -sum(diag(codon_mutation_matrix) * codon.freq.by.gene, na.rm=TRUE)
        codon_mutation_matrix_scaled = codon_mutation_matrix * (1/scale.factor)
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
            rate.Q_codon.list <- as.list(ncats)
            for(cat.index in 1:ncats){
                if(k.levels > 0){
                    aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
                }else{
                    aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
                }
                rate.Q_codon.list[[cat.index]] <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[cat.index], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
            }
            for(cat.index in 1:ncats){
                Q_codon_array <- rate.Q_codon.list[[cat.index]]
                for(aa.index in 1:21){
                    if(diploid == TRUE){
                        Q_codon_array[,,.unique.aa[aa.index]] <- 2 * Ne * (codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[aa.index]])
                    }else{
                        Q_codon_array[,,.unique.aa[aa.index]] <- Ne * (codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[aa.index]])
                    }
                    diag(Q_codon_array[,,.unique.aa[aa.index]]) <- 0
                    diag(Q_codon_array[,,.unique.aa[aa.index]]) <- -rowSums(Q_codon_array[,,.unique.aa[aa.index]])
                    Q_codon_array[,,.unique.aa[aa.index]] <- Q_codon_array[,,.unique.aa[aa.index]]
                }
                rate.Q_codon.list[[cat.index]] <- Q_codon_array
            }

            #Put the na.rm=TRUE bit here just in case -- when the amino acid is a stop codon, there is a bunch of NaNs. Should be fixed now.
            phy.sort <- reorder(phy, "pruningwise")

            #Generate matrix of root frequencies for each optimal AA:
            root.p_array <- matrix(codon.freq.by.aa, nrow=dim(Q_codon_array)[2], ncol=21)
            root.p_array <- t(root.p_array)
            root.p_array <- root.p_array / rowSums(root.p_array)
            rownames(root.p_array) <- .unique.aa

            if(model.to.simulate.under == "selac" | model.to.simulate.under == "fmutsel"){
                if(model.to.simulate.under == "selac"){
                    simulation.model.info <- GetSelacSimulateInfo(selac.obj=selac.obj2, partition.number=partition.number)
                    if(!is.null(simulation.model.info$gamma.rates)){
                        for(i in sequence(nsites)){
                            site.rate <- sample(1:4, 1, prob=weights.k)
                            simulated.site.phi.cat <- c(simulated.site.phi.cat, site.rate)
                            Q_codon_array <- rate.Q_codon.list[[site.rate]]
                            Q_codon_recon <- Q_codon_array[,,aa.optim_array[i]]
                            Q_codon_array <- simulation.model.info$Q_matrix[[site.rate]]
                            Q_codon_sim <- Q_codon_array[,,aa.optim_array[i]]
                            interval.recon_array <- cbind(interval.recon_array, GetTipIntervalStateSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, root.p=root.p_array[aa.optim_array[i],], taxon.to.drop=taxon.to.drop, Q.to.reconstruct=Q_codon_recon, Q.to.simulate=Q_codon_sim,  model.to.reconstruct.under=model.to.reconstruct.under, model.to.simulate.under=model.to.simulate.under))
                        }
                    }else{
                        for(i in sequence(nsites)){
                            site.rate <- sample(1:4, 1, prob=weights.k)
                            Q_codon_array <- rate.Q_codon.list[[site.rate]]
                            Q_codon_recon <- Q_codon_array[,,aa.optim_array[i]]
                            Q_codon_sim <- simulation.model.info$Q_matrix[,,aa.optim_array[i]]
                            interval.recon_array <- cbind(interval.recon_array, GetTipIntervalStateSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, root.p=root.p_array[aa.optim_array[i],], taxon.to.drop=taxon.to.drop, Q.to.reconstruct=Q_codon_recon, Q.to.simulate=Q_codon_sim,  model.to.reconstruct.under=model.to.reconstruct.under, model.to.simulate.under=model.to.simulate.under))
                        }
                    }
                }else{
                    simulation.model.info <- GetFMutSelSimulateInfo(selac.obj=selac.obj2, partition.number=partition.number)
                    for(i in sequence(nsites)){
                        site.rate <- sample(1:4, 1, prob=weights.k)
                        Q_codon_array <- rate.Q_codon.list[[site.rate]]
                        Q_codon_recon <- Q_codon_array[,,aa.optim_array[i]]
                        Q_codon_sim <- simulation.model.info$Q_matrix
                        interval.recon_array <- cbind(interval.recon_array, GetTipIntervalStateSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, root.p=root.p_array[aa.optim_array[i],], taxon.to.drop=taxon.to.drop, Q.to.reconstruct=Q_codon_recon, Q.to.simulate=Q_codon_sim,  model.to.reconstruct.under=model.to.reconstruct.under, model.to.simulate.under=model.to.simulate.under))
                    }
                }
            }else{
                simulation.model.info <- GetGtrSimulateInfo(selac.obj=selac.obj2, partition.number=partition.number)
                for(i in sequence(nsites)){
                    site.rate <- sample(1:4, 1, prob=weights.k)
                    Q_codon_array <- rate.Q_codon.list[[site.rate]]
                    Q_codon <- Q_codon_array[,,aa.optim_array[i]]
                    Q_nuc.mat <- simulation.model.info$Q_mat
                    site.rate <- sample(1:4, 1, prob=simulation.model.info$gamma.weights)
                    interval.recon_array <- cbind(interval.recon_array, GetTipIntervalStateSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, root.p=root.p_array[aa.optim_array[i],], taxon.to.drop=taxon.to.drop, Q.to.reconstruct=Q_codon_array[,,aa.optim_array[i]], Q.to.simulate=Q_nuc.mat*simulation.model.info$gamma.rates[site.rate], model.to.reconstruct.under=model.to.reconstruct.under, model.to.simulate.under=model.to.simulate.under))
                }
            }
        }else{
            if(k.levels > 0){
                aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
            }else{
                aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
            }
            Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999)
            #Finish the Q_array codon mutation matrix multiplication here:
            for(k in 1:21){
                if(diploid == TRUE){
                    Q_codon_array[,,.unique.aa[k]] = (2 * Ne) * codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[k]]
                }else{
                    Q_codon_array[,,.unique.aa[k]] = Ne * codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[k]]
                }
                diag(Q_codon_array[,,.unique.aa[k]]) = 0
                diag(Q_codon_array[,,.unique.aa[k]]) = -rowSums(Q_codon_array[,,.unique.aa[k]])
            }

            #Put the na.rm=TRUE bit here just in case -- when the amino acid is a stop codon, there is a bunch of NaNs. Should be fixed now.
            #scale.factor <- -sum(Q_codon_array[DiagArray(dim(Q_codon_array))] * equilibrium.codon.freq, na.rm=TRUE)
            phy <- reorder(phy, "pruningwise")

            #Generate matrix of root frequencies for each optimal AA:
            root.p_array <- matrix(codon.freq.by.aa, nrow=dim(Q_codon_array)[2], ncol=21)
            root.p_array <- t(root.p_array)
            root.p_array <- root.p_array / rowSums(root.p_array)
            rownames(root.p_array) <- .unique.aa
            phy.sort <- reorder(phy, "pruningwise")

            if(model.to.simulate.under == "selac"){
                simulation.model.info <- GetSelacSimulateInfo(selac.obj=selac.obj2, partition.number=partition.number)
                if(!is.null(simulation.model.info$gamma.rates)){
                    for(i in sequence(nsites)){
                        Q_codon_recon <- Q_codon_array[,,aa.optim_array[i]]
                        site.rate <- sample(1:4, 1, prob=simulation.model.info$gamma.weights)
                        simulated.site.phi.cat <- c(simulated.site.phi.cat, site.rate)
                        Q_codon_array <- simulation.model.info$Q_matrix[[site.rate]]
                        Q_codon_sim <- Q_codon_array[,,aa.optim_array[i]]
                        interval.recon_array <- cbind(interval.recon_array, GetTipIntervalStateSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, root.p=root.p_array[aa.optim_array[i],], taxon.to.drop=taxon.to.drop, Q.to.reconstruct=Q_codon_recon, Q.to.simulate=Q_codon_sim,  model.to.reconstruct.under=model.to.reconstruct.under, model.to.simulate.under=model.to.simulate.under))
                    }
                }else{
                    for(i in sequence(nsites)){
                        Q_codon_recon <- Q_codon_array[,,aa.optim_array[i]]
                        Q_codon_sim <- simulation.model.info$Q_matrix[,,aa.optim_array[i]]
                        interval.recon_array <- cbind(interval.recon_array, GetTipIntervalStateSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, root.p=root.p_array[aa.optim_array[i],], taxon.to.drop=taxon.to.drop, Q.to.reconstruct=Q_codon_recon, Q.to.simulate=Q_codon_sim,  model.to.reconstruct.under=model.to.reconstruct.under, model.to.simulate.under=model.to.simulate.under))
                    }
                }
            }else{
                simulation.model.info <- GetGtrSimulateInfo(selac.obj=selac.obj2, partition.number=partition.number)
                for(i in sequence(nsites)){
                    Q_codon <- Q_codon_array[,,aa.optim_array[i]]
                    Q_nuc.mat <- simulation.model.info$Q_mat
                    site.rate <- sample(1:4, 1, prob=simulation.model.info$gamma.weights)
                    interval.recon_array <- cbind(interval.recon_array, GetTipIntervalStateSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, root.p=root.p_array[aa.optim_array[i],], taxon.to.drop=taxon.to.drop, Q.to.reconstruct=Q_codon_array[,,aa.optim_array[i]], Q.to.simulate=Q_nuc.mat*simulation.model.info$gamma.rates[site.rate], model.to.reconstruct.under=model.to.reconstruct.under, model.to.simulate.under=model.to.simulate.under))
                }
            }
        }
    }else{
        for(partition.index in partition.number:partition.number){
            pars <- c(selac.obj1$mle.pars[partition.index,])
            gene.tmp <- read.dna(partitions[partition.index], format='fasta')
            if(!is.null(fasta.rows.to.keep)){
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
            }else{
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
            }

            codon.data <- DNAbinToNucleotideNumeric(gene.tmp)
            codon.data <- codon.data[phy$tip.label,]
            codon.freq.by.gene <- selac.obj1$empirical.base.freqs[[partition.index]]
            codon.freq.by.aa=NULL
        }

        nsites <- dim(codon.data)[2]-1
        prop.interval <- seq(0, 1, by=0.05)

        interval.recon_array <- c()

        if(include.gamma == TRUE){
            shape = pars[1]
            pars = pars[-1]
        }
        if(length(pars)==0){
            transition.rates <- 1
        }else{
            transition.rates <- pars[1:length(pars)]
        }
        nsites <- dim(codon.data)[2]-1
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(transition.rates, model=nuc.model)
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

            #Rescaling Q matrix in order to have a 1 nucleotide change per site if the branch length was 1:
            if(is.null(codon.freq.by.gene)) {
                #Generate matrix of equal frequencies for each site:
                root.p_array <- rep(0.25, 4)
            }else{
                root.p_array <- codon.freq.by.gene
            }

            phy.sort <- reorder(phy, "pruningwise")

            diag(nuc.mutation.rates) = 0
            nuc.mutation.rates = t(nuc.mutation.rates * root.p_array)
            diag(nuc.mutation.rates) = -rowSums(nuc.mutation.rates)
            scale.factor <- -sum(diag(nuc.mutation.rates) * root.p_array)
            Q_array <- nuc.mutation.rates * (1/scale.factor)
            if(model.to.simulate.under == "selac" | model.to.simulate.under == "fmutsel"){
                if(model.to.simulate.under == "selac"){
                    simulation.model.info <- GetSelacSimulateInfo(selac.obj=selac.obj2, partition.number=partition.number)
                    marginal.recon_array <- array(1, dim=c(Ntip(phy)+Nnode(phy), 4, nsites))
                    for(i in sequence(nsites)){
                        site.rate <- sample(1:4, 1, prob=weights.k)
                        Q_codon_recon <- Q_array * rates.k[site.rate]
                        Q_codon_array <- simulation.model.info$Q_matrix[[site.rate]]
                        Q_codon_sim <- Q_codon_array
                        marginal.recon_array[,,i] <- GetTipIntervalStateSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, root.p=root.p_array, taxon.to.drop=taxon.to.drop, Q.to.reconstruct=Q_codon_recon, Q.to.simulate=Q_codon_sim, model.to.reconstruct.under=model.to.reconstruct.under, model.to.simulate.under=model.to.simulate.under)
                    }
                    if(!is.null(simulation.model.info$gamma.rates)){
                        tot.interval <- phy$edge.length[phy$edge[,2]==taxon.to.drop]
                        prop.interval <- seq(0,1 , by=0.05)
                        time.interval <- tot.interval * prop.interval
                        focal <- taxon.to.drop
                        focalRows <- which(phy$edge[,2]==focal)
                        codon.count <- 1
                        Q_codon_array <- simulation.model.info$Q_matrix
                        end.codon <- seq(3,nsites, by=3)
                        for(i in seq(1, nsites, by=3)){
                            site.rate <- sample(1:4, 1, prob=weights.k)
                            Q_codon_array <- simulation.model.info$Q_matrix[[site.rate]]
                            Q_codon_sim <- Q_codon_array[,,simulation.model.info$aa.optim_array[codon.count]]
                            #Adds a slight bias to this, but cannot have stop codons in the selac case -- nothing holding this back
                            stop.codon <- 1
                            while(stop.codon == 1){
                                focal.starting.state <- c()
                                for(k in i:end.codon[codon.count]){
                                    focal.starting.state <- c(focal.starting.state, sample(1:4, 1, prob=marginal.recon_array[phy$edge[focalRows,1],,k]))
                                }
                                focal.starting.state.converted <- which(.codon.name==paste(as.vector(n2s(focal.starting.state-1)), collapse=""))
                                if(!focal.starting.state.converted == 49 | !focal.starting.state.converted == 51 | !focal.starting.state.converted == 57){
                                    stop.codon = 0
                                }
                            }
                            reconstructed.sequence.codon <- c(focal.starting.state.converted)
                            for (time.index in 2:length(time.interval)){
                                p <- expm(Q_codon_sim * time.interval[time.index], method="Ward77")[focal.starting.state.converted, ]
                                focal.starting.state.converted <- sample.int(64, size = 1, FALSE, prob = p)
                                reconstructed.sequence.codon <- c(reconstructed.sequence.codon, focal.starting.state.converted)
                            }
                            interval.recon_array <- cbind(interval.recon_array, reconstructed.sequence.codon)
                            codon.count <- codon.count+1
                        }
                    }else{
                        tot.interval <- phy$edge.length[phy$edge[,2]==taxon.to.drop]
                        prop.interval <- seq(0,1 , by=0.05)
                        time.interval <- tot.interval * prop.interval
                        focal <- taxon.to.drop
                        focalRows <- which(phy$edge[,2]==focal)
                        codon.count <- 1
                        Q_codon_array <- simulation.model.info$Q_matrix
                        end.codon <- seq(3,nsites, by=3)
                        for(i in seq(1, nsites, by=3)){
                            site.rate <- sample(1:4, 1, prob=weights.k)
                            Q_codon_sim <- Q_codon_array[,,simulation.model.info$aa.optim_array[codon.count]]
                            #Adds a slight bias to this, but cannot have stop codons in the selac case -- nothing holding this back
                            stop.codon <- 1
                            while(stop.codon == 1){
                                focal.starting.state <- c()
                                for(k in i:end.codon[codon.count]){
                                    focal.starting.state <- c(focal.starting.state, sample(1:4, 1, prob=marginal.recon_array[phy$edge[focalRows,1],,k]))
                                }
                                focal.starting.state.converted <- which(.codon.name==paste(as.vector(n2s(focal.starting.state-1)), collapse=""))
                                if(!focal.starting.state.converted == 49 | !focal.starting.state.converted == 51 | !focal.starting.state.converted == 57){
                                    stop.codon = 0
                                }
                            }
                            reconstructed.sequence.codon <- c(focal.starting.state.converted)
                            for (time.index in 2:length(time.interval)){
                                p <- expm(Q_codon_sim * time.interval[time.index], method="Ward77")[focal.starting.state.converted, ]
                                focal.starting.state.converted <- sample.int(64, size = 1, FALSE, prob = p)
                                reconstructed.sequence.codon <- c(reconstructed.sequence.codon, focal.starting.state.converted)
                            }
                            interval.recon_array <- cbind(interval.recon_array, reconstructed.sequence.codon)
                            codon.count <- codon.count+1
                        }
                    }
                }else{
                    print("why are we here?")
                }
            }else{
                simulation.model.info <- GetGtrSimulateInfo(selac.obj=selac.obj2, partition.number=partition.number)
                for(i in sequence(nsites)){
                    site.rate <- sample(1:4, 1, prob=weights.k)
                    Q_codon_recon <- Q_array * rates.k[site.rate]
                    Q_nuc.mat <- simulation.model.info$Q_mat
                    interval.recon_array <- cbind(interval.recon_array, GetTipIntervalStateSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, root.p=root.p_array, taxon.to.drop=taxon.to.drop, Q.to.reconstruct=Q_codon_recon, Q.to.simulate=Q_nuc.mat*simulation.model.info$gamma.rates[site.rate], model.to.reconstruct.under=model.to.reconstruct.under, model.to.simulate.under=model.to.simulate.under))
                }
            }
        }
    }
    obj <- NULL
    obj$interval.recon_array <- interval.recon_array
    if(length(simulated.site.phi.cat)==0){
        obj$site.gamma.indicator <- NULL
    }else{
        obj$site.gamma.indicator <- simulated.site.phi.cat
    }
    return(obj)
}


#' @title Model adequacy simulation
#'
#' @description
#' Performs a single model adequacy simulation
#'
#' @param model.to.reconstruct.under Specifies the model that the internal nodes are to be reconstructed under assuming a single tip is pruned from the tree.
#' @param model.to.simulate.under Specifies the model that the simulation will be conducted along the pruned tip.
#' @param selac.obj.to.reconstruct The selac output object that contains the model parameters to be used in the reconstruction.
#' @param selac.obj.to.simulate The selac output object that contains the model parameters to be used in the simulation.
#' @param aa.optim.input A list of optimal amino acids with each list element designating a character vector for each gene. The optimal amino acids be the MLE from a selac run (default) or a list of user defined optimal A.A.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.
#' @param taxon.to.drop Specifies the tip based on the number in the phy object to be removed and simulated.
#' @param partition.number Specifies the partition number to conduct the model adequacy test.
#' @param numcode The ncbi genetic code number for translation. By default the standard (numcode=1) genetic code is used.
#' @param for.gtr.only A selac object that can be used as the reference optimal AA for when the adequacy of a GTR+G model is tested only.
#'
#' @details
#' Performs a single model adequacy simulation. The test prunes out a user-specified taxon from the tree, performs site data reconstruction for all nodes in the tree under a user-specified model, then simulates the expected data of the pruned taxon according to a user-specified model along uniformly sampled points along the branch. The functionality of the reconstructed sequence is also calculated along the way to see how functionality changes as the simulation reaches the end of the known branch length. The output is a vector with elements containing the functionality of the simulated points along equally spaced sampling points along the known branch length (i.e., edge.length * seq(0, 1, by=0.05))
GetAdequateSelac <- function(model.to.reconstruct.under, model.to.simulate.under, selac.obj.to.reconstruct, selac.obj.to.simulate, aa.optim.input=NULL, fasta.rows.to.keep=NULL, taxon.to.drop=4, partition.number=55, numcode=1, for.gtr.only=NULL){
    prop.intervals <- seq(0,1, by=0.05)

    if(model.to.reconstruct.under == "gtr" & model.to.simulate.under == "selac"){
        simulated.across.intervals.and.sites <- GetIntervalSequencesAllSites(model.to.simulate.under=model.to.simulate.under, model.to.reconstruct.under=model.to.reconstruct.under, selac.obj1=selac.obj.to.reconstruct, selac.obj2=selac.obj.to.simulate, aa.optim.input=aa.optim.input, fasta.rows.to.keep=fasta.rows.to.keep, taxon.to.drop=taxon.to.drop, partition.number=partition.number)
        functionality.taxon <- c()
        for(interval.index in 1:length(prop.intervals)){
            reconstructed.sequence <- c()
            for(site.index in 1:selac.obj.to.simulate$nsites[partition.number]){
                reconstructed.sequence <- c(reconstructed.sequence, .aa.translation[[numcode]][simulated.across.intervals.and.sites$interval.recon_array[interval.index,site.index]])
                reconstructed.sequence <- unname(reconstructed.sequence)
            }
            selac.obj.to.reconstruct <- selac.obj.to.simulate
            if(selac.obj.to.simulate$include.gamma == TRUE){
                rates.cat <- LaguerreQuad(selac.obj.to.simulate$mle.pars[1,length(selac.obj.to.simulate$mle.pars[1,])], ncats=4)[1:4]
                functionality.taxon.interval <- GetFunctionalityModelAdequacy(gene.length=length(reconstructed.sequence), aa.data=reconstructed.sequence, optimal.aa=selac.obj.to.reconstruct$aa.optim[[partition.number]], alpha=selac.obj.to.reconstruct$mle.pars[1,2], beta=selac.obj.to.reconstruct$mle.pars[1,3], gamma=selac.obj.to.reconstruct$volume.fixed.value, gp=rates.cat[simulated.across.intervals.and.sites$site.gamma.indicator], aa.properties=selac.obj.to.reconstruct$aa.properties)
            }else{
                functionality.taxon.interval <- GetFunctionalityModelAdequacy(gene.length=length(reconstructed.sequence), aa.data=reconstructed.sequence, optimal.aa=selac.obj.to.reconstruct$aa.optim[[partition.number]], alpha=selac.obj.to.reconstruct$mle.pars[1,2], beta=selac.obj.to.reconstruct$mle.pars[1,3], gamma=selac.obj.to.reconstruct$volume.fixed.value, gp=NULL, aa.properties=selac.obj.to.reconstruct$aa.properties)
            }

            functionality.taxon <- c(functionality.taxon, functionality.taxon.interval)
        }
    }

    if(model.to.reconstruct.under == "gtr" & model.to.simulate.under == "gtr"){
        simulated.across.intervals.and.sites <- GetIntervalSequencesAllSites(model.to.simulate.under=model.to.simulate.under, model.to.reconstruct.under=model.to.reconstruct.under, selac.obj1=selac.obj.to.reconstruct, selac.obj2=selac.obj.to.simulate, aa.optim.input=aa.optim.input, fasta.rows.to.keep=fasta.rows.to.keep, taxon.to.drop=taxon.to.drop, partition.number=partition.number)
        reconstructed.sequence <- c()
        functionality.taxon <- c()
        for(interval.index in 1:length(prop.intervals)){
            end.site <- seq(3,selac.obj.to.simulate$nsites[partition.number], by=3)
            count <- 1
            reconstructed.sequence.codon <- c()
            for(start.site in seq(1, selac.obj.to.simulate$nsites[partition.number], by=3)){
                reconstructed.sequence.codon <- c(reconstructed.sequence.codon, .aa.translation[[numcode]][which(.codon.name==paste(as.vector(n2s(simulated.across.intervals.and.sites$interval.recon_array[interval.index,start.site:end.site[count]]-1)), collapse=""))])
                count <- count + 1
            }
            functionality.taxon.interval <- GetFunctionalityModelAdequacy(gene.length=length(reconstructed.sequence.codon), aa.data=reconstructed.sequence.codon, optimal.aa=for.gtr.only$aa.optim[[partition.number]], alpha=for.gtr.only$mle.pars[1,2], beta=for.gtr.only$mle.pars[1,3], gamma=for.gtr.only$volume.fixed.value, aa.properties=for.gtr.only$aa.properties)
            functionality.taxon <- c(functionality.taxon, functionality.taxon.interval)
        }
    }

    if(model.to.reconstruct.under == "selac" & model.to.simulate.under == "gtr" | model.to.reconstruct.under == "selac" & model.to.simulate.under == "selac" | model.to.reconstruct.under == "selac" & model.to.simulate.under == "fmutsel"){
        simulated.across.intervals.and.sites <- GetIntervalSequencesAllSites(model.to.simulate.under=model.to.simulate.under, model.to.reconstruct.under=model.to.reconstruct.under, selac.obj1=selac.obj.to.reconstruct, selac.obj2=selac.obj.to.simulate, aa.optim.input=aa.optim.input, fasta.rows.to.keep=fasta.rows.to.keep, taxon.to.drop=taxon.to.drop, partition.number=partition.number)
        functionality.taxon <- c()
        for(interval.index in 1:length(prop.intervals)){
            reconstructed.sequence <- c()
            for(site.index in 1:selac.obj.to.reconstruct$nsites[partition.number]){
                reconstructed.sequence <- c(reconstructed.sequence, .aa.translation[[numcode]][simulated.across.intervals.and.sites$interval.recon_array[interval.index,site.index]])
                reconstructed.sequence <- unname(reconstructed.sequence)
            }
            if(selac.obj.to.simulate$include.gamma == TRUE){
                if(model.to.simulate.under == "selac"){
                    rates.cat <- LaguerreQuad(selac.obj.to.simulate$mle.pars[1,length(selac.obj.to.simulate$mle.pars[1,])], ncats=4)[1:4]
                    functionality.taxon.interval <- GetFunctionalityModelAdequacy(gene.length=length(reconstructed.sequence), aa.data=reconstructed.sequence, optimal.aa=selac.obj.to.reconstruct$aa.optim[[partition.number]], alpha=selac.obj.to.reconstruct$mle.pars[1,2], beta=selac.obj.to.reconstruct$mle.pars[1,3], gamma=selac.obj.to.reconstruct$volume.fixed.value, gp=rates.cat[simulated.across.intervals.and.sites$site.gamma.indicator], aa.properties=selac.obj.to.reconstruct$aa.properties)
                }else{
                    functionality.taxon.interval <- GetFunctionalityModelAdequacy(gene.length=length(reconstructed.sequence), aa.data=reconstructed.sequence, optimal.aa=selac.obj.to.reconstruct$aa.optim[[partition.number]], alpha=selac.obj.to.reconstruct$mle.pars[1,2], beta=selac.obj.to.reconstruct$mle.pars[1,3], gamma=selac.obj.to.reconstruct$volume.fixed.value, gp=NULL, aa.properties=selac.obj.to.reconstruct$aa.properties)
                }
            }else{
                functionality.taxon.interval <- GetFunctionalityModelAdequacy(gene.length=length(reconstructed.sequence), aa.data=reconstructed.sequence, optimal.aa=selac.obj.to.reconstruct$aa.optim[[partition.number]], alpha=selac.obj.to.reconstruct$mle.pars[1,2], beta=selac.obj.to.reconstruct$mle.pars[1,3], gamma=selac.obj.to.reconstruct$volume.fixed.value, gp=NULL, aa.properties=selac.obj.to.reconstruct$aa.properties)
            }
            functionality.taxon <- c(functionality.taxon, functionality.taxon.interval)
        }
    }
    return(functionality.taxon)
}



#' @title Parallel model adequacy test
#'
#' @description
#' Performs model adequacy test using multiple cores
#'
#' @param nreps Specifies the number of repeated model adequact simulations.
#' @param n.cores Specifies the number of cores you want to use.
#' @param model.to.reconstruct.under Specifies the model that the internal nodes are to be reconstructed under assuming a single tip is pruned from the tree.
#' @param model.to.simulate.under Specifies the model that the simulation will be conducted along the pruned tip.
#' @param selac.obj.to.reconstruct The selac output object that contains the model parameters to be used in the reconstruction.
#' @param selac.obj.to.simulate The selac output object that contains the model parameters to be used in the simulation.
#' @param aa.optim.input A list of optimal amino acids with each list element designating a character vector for each gene. The optimal amino acids be the MLE from a selac run (default) or a list of user defined optimal A.A.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.
#' @param taxon.to.drop Specifies the tip based on the number in the phy object to be removed and simulated.
#' @param partition.number Specifies the partition number to conduct the model adequacy test.
#' @param numcode The ncbi genetic code number for translation. By default the standard (numcode=1) genetic code is used.
#' @param for.gtr.only A selac object that can be used as the reference optimal AA for when the adequacy of a GTR+G model is tested only.
#'
#' @details
#' Performs a parallelized analysis of the model adequacy test. The test prunes out a user-specified taxon from the tree, performs site data reconstruction for all nodes in the tree under a user-specified model, then simulates the expected data of the pruned taxon according to a user-specified model along uniformly sampled points along the branch. The functionality of the reconstructed sequence is also calculated along the way to see how functionality changes as the simulation reaches the end of the known branch length. The output is a list with elements equally the number of repititions. Each element contains the functionality of the simulated points along equally spaced sampling points along the known branch length (i.e., edge.length * seq(0, 1, by=0.05))
GetAdequateManyReps <- function(nreps, n.cores, model.to.reconstruct.under="selac", model.to.simulate.under="gtr", selac.obj.to.reconstruct, selac.obj.to.simulate, aa.optim.input=NULL, fasta.rows.to.keep=NULL, taxon.to.drop=2, partition.number=17, numcode=1, for.gtr.only=NULL){
    RunMany <- function(nrep){
        tmp <- GetAdequateSelac(model.to.simulate.under=model.to.simulate.under, model.to.reconstruct.under=model.to.reconstruct.under, selac.obj.to.reconstruct=selac.obj.to.reconstruct, selac.obj.to.simulate=selac.obj.to.simulate, aa.optim.input=aa.optim.input, fasta.rows.to.keep=fasta.rows.to.keep, taxon.to.drop=taxon.to.drop, partition.number=partition.number, numcode=numcode, for.gtr.only=for.gtr.only)
        return(tmp)
    }
    pp <- mclapply(1:nreps, RunMany, mc.cores=n.cores)
    return(pp)
}



#PlotAdequacyResults <- function(adequate.obj.gtr, adequate.obj.selac.nog, adequate.obj.selac.wg, alpha=.1, known.functionality=0.9293696, file.name="modelAdequacy.pdf"){
#    prop.interval <- seq(0,1 , by=0.05)
#    pdf(file.name)
#    plot(prop.interval, adequate.obj.gtr[[1]], ylab="", xlab="", pch=19, xlim=c(0,1), ylim=c(.5, 1), axes=FALSE, col=0)
#
#    for(rep.index in 1:100){
#        lines(prop.interval, adequate.obj.gtr[[rep.index]], col=add.alpha("blue", alpha=alpha))
#    }
#
#    for(rep.index in 1:100){
#        lines(prop.interval, adequate.obj.selac.nog[[rep.index]], col=add.alpha("red", alpha=alpha))
#    }
#
#    for(rep.index in 1:100){
#        lines(prop.interval, adequate.obj.selac.wg[[rep.index]], col=add.alpha("green", alpha=alpha))
#    }
#
#    title(xlab = "Proportion of true edge length", line=2)
#    title(ylab = "Functionality", line=3)
#    par(tck=.01)
#    axis(2, at = seq(.5, 1, by = .1), las =1, lwd=1, labels=TRUE, mgp=c(1,.5,0))
#    axis(1, at = seq(0, 1, by = .2), las =1, lwd=1, labels=TRUE, mgp=c(1,.5,0))
#    abline(h=known.functionality, lty=2)
#    dev.off()
#
#}

#GetKnownFunctionality(result, partition.number=55, fasta.rows.to.keep=1:7, taxon.to.do=2)

## Adequacy of selac + Gamma ##
#load("yeastRokasSelacUNREST.Rdata")
#selac.nog <- result
#load("yeastRokasSelacUNRESTgamma.Rdata")
#selac.wg <- result
#load("yeastRokasGTRgamma.Rdata")
#gtr.g <- result
#gtr.g$gamma.type = "median"

#pp <- GetAdequateManyReps(nreps=10, n.cores=1, model.to.reconstruct.under="selac", model.to.simulate.under="selac", selac.obj.to.reconstruct=selac.nog, selac.obj.to.simulate=selac.wg, aa.optim.input=NULL, fasta.rows.to.keep=1:7, taxon.to.drop=6, partition.number=17, numcode=1)
#save(pp, file="adequacy.gtr.g.Rsave")
#pp <- GetAdequateManyReps(nreps=100, n.cores=4, selac.obj=selac.nog, selac.obj2=selac.wg, aa.optim.input=NULL, fasta.rows.to.keep=1:7, taxon.to.drop=6, partition.number=17, numcode=1, data.type="codon")
#save(pp, file="adequacy.selac.nog.Rsave")
#pp <- GetAdequateManyReps(nreps=100, n.cores=4, selac.obj=selac.wg, selac.obj2=selac.wg, aa.optim.input=NULL, fasta.rows.to.keep=1:7, taxon.to.drop=6, partition.number=17, numcode=1, data.type="codon")
#save(pp, file="adequacy.selac.wg.Rsave")
###############################

## Adequacy of selac + Gamma ##
#load("yeastRokasSelacUNRESTgamma.Rdata")
#selac.obj2 <- result
#load("yeastRokasGTRgamma.Rdata")
#selac.obj <- result
#selac.obj$gamma.type <- "median"
#pp <- GetAdequateManyReps(100, n.cores=4, selac.obj=result, selac.obj2=NULL, fasta.rows.to.keep=1:7, taxon.to.drop=2, partition.number=17, data.type="codon")
###############################
#load("adequacy.gtr.g.Rsave")
#adequate.obj.gtr <- pp
#load("adequacy.selac.nog.Rsave")
#adequate.obj.selac.nog <- pp
#load("adequacy.selac.wg.Rsave")
#adequate.obj.selac.wg <- pp
#PlotAdequacyResults(adequate.obj.gtr=adequate.obj.gtr, adequate.obj.selac.nog=adequate.obj.selac.nog, adequate.obj.selac.wg=adequate.obj.selac.wg, file.name="adequacy_selac_Scer.pdf")
