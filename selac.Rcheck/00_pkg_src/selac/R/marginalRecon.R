
######################################################################################################################################
######################################################################################################################################
### Marginal reconstruction across genes and nodes in a tree
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu


GetMarginalSingleSite <- function(charnum=1, codon.data, phy, Q, root.p=NULL, taxon.to.drop=NULL){
    
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    nl <- nrow(Q)
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
            if(nl > 4){
                liks[i,c(49, 51, 57)] <- 0
            }
        }
    }
    
    if(!is.null(taxon.to.drop)){
        liks[taxon.to.drop,] <- 1
        if(nl > 4){
            liks[taxon.to.drop,c(49, 51, 57)] <- 0
        }
    }
    
    TIPS <- 1:nb.tip
    anc <- unique(phy$edge[,1])
    
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
    
    #A temporary likelihood matrix so that the original does not get written over:
    liks.down <- liks
    #A transpose of Q for assessing probability of j to i, rather than i to j:
    tranQ <- t(Q)
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
                v <- v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks.down[desNodes[desIndex],]
            }
        }
        comp[focal] <- sum(v)
        liks.down[focal, ] <- v/comp[focal]
    }
    root <- nb.tip + 1L
    #Enter the root defined root probabilities if they are supplied by the user:
    equil.root <- NULL
    for(i in 1:ncol(Q)){
        posrows <- which(Q[,i] >= 0)
        rowsum <- sum(Q[posrows,i])
        poscols <- which(Q[i,] >= 0)
        colsum <- sum(Q[i,poscols])
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
                    v <- v*expm(Q * phy$edge.length[sisterRows[sisterIndex]], method=c("Ward77")) %*% liks.down[sisterNodes[sisterIndex],]
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
        focalRows<-which(phy$edge[,2]==focal)
        #Now you are assessing the change along the branch subtending the focal by multiplying the probability of
        #everything at and above focal by the probability of the mother and all the sisters given time t:
        v <- liks.down[focal,]*expm(tranQ * phy$edge.length[focalRows], method=c("Ward77")) %*% liks.up[focal,]
        comp[focal] <- sum(v)
        liks.final[focal, ] <- v/comp[focal]
    }

    #Just add in the marginal at the root calculated on the original downpass or if supplied by the user:
    liks.final[root,] <- liks.down[root,] * root.p
    root.final <- liks.down[root,] * root.p
    comproot <- sum(root.final)
    liks.final[root,] <- root.final/comproot
    return(liks.final)
}


GetMarginalGene <- function(pars, codon.data, phy, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, aa.optim_array, numcode, nuc.model="UNREST", include.gamma=FALSE, gamma.type="quadrature", ncats=4, k.levels=0, diploid=TRUE, taxon.to.drop=NULL, aa.properties=NULL, data.type="codon"){
    
    if(data.type == "nucleotide"){
        nsites <- dim(codon.data)[2]-1
        marginal.recon_array <- array(1, dim=c(Ntip(phy)+Nnode(phy), 4, nsites))
        
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
            
            for (i in sequence(nsites)) {
                
                marginal.recon_ncats_array <- array(1, dim=c(Ntip(phy)+Nnode(phy), 4, ncats))
                for(k.cat in sequence(ncats)){
                    marginal.recon_ncats_array[,,k.cat] <- GetMarginalSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, Q=nuc.mutation.rates*rates.k[k.cat], root.p=root.p_array, taxon.to.drop=taxon.to.drop) * weights.k
                }
                marginal.recon_array[,,i] <- apply(marginal.recon_ncats_array, 2, rowSums)/ rowSums(apply(marginal.recon_ncats_array, 2, rowSums))
            }
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
            for (i in sequence(nsites)) {
                
                phy.sort <- reorder(phy, "pruningwise")
                
                marginal.recon_array[,,i] <- GetMarginalSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, Q=nuc.mutation.rates, root.p=root.p_array, taxon.to.drop=taxon.to.drop)
            }
        }
    }else{
        codon.index.matrix = CreateCodonMutationMatrixIndex()
        nsites <- dim(codon.data)[2]-1
        marginal.recon_array <- array(1, dim=c(Ntip(phy)+Nnode(phy), 64, nsites))
        
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
        gamma <-  0.0003990333
        
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
            for (i in sequence(nsites)) {
                marginal.recon_ncats_array <- array(1, dim=c(Ntip(phy)+Nnode(phy), 64, ncats))
                for(k.cat in sequence(ncats)){
                    if(k.levels > 0){
                        aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
                    }else{
                        aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
                    }
                    Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[k.cat], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999)
                    
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
                    phy <- reorder(phy, "pruningwise")
                    
                    #Generate matrix of root frequencies for each optimal AA:
                    root.p_array <- matrix(codon.freq.by.aa, nrow=dim(Q_codon_array)[2], ncol=21)
                    root.p_array <- t(root.p_array)
                    root.p_array <- root.p_array / rowSums(root.p_array)
                    rownames(root.p_array) <- .unique.aa
                    
                    phy.sort <- reorder(phy, "pruningwise")
                    
                    marginal.recon_ncats_array[,,k.cat] <- GetMarginalSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, Q=Q_codon_array[,,aa.optim_array[i]], root.p=root.p_array[aa.optim_array[i],], taxon.to.drop=taxon.to.drop) * weights.k
                }
                marginal.recon_array[,,i] <- apply(marginal.recon_ncats_array, 2, rowSums)/ rowSums(apply(marginal.recon_ncats_array, 2, rowSums))
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
            
            for (i in sequence(nsites)) {
                marginal.recon_array[,,i] <- GetMarginalSingleSite(charnum=i, codon.data=codon.data, phy=phy.sort, Q=Q_codon_array[,,aa.optim_array[i]], root.p=root.p_array[aa.optim_array[i],], taxon.to.drop=taxon.to.drop)
            }
        }
    }
    return(marginal.recon_array)
}


#' @title Get marginal reconstruction all genes
#'
#' @description
#' Calculates the marginal probability of each codon at all sites across all genes
#'
#' @param selac.obj An object of class SELAC.
#' @param aa.optim.input A list of optimal amino acids with each list element designating a character vector for each gene. The optimal amino acids be the MLE from a selac run (default) or a list of user defined optimal A.A.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.
#' @param taxon.to.drop A single taxon (defined by number in phy object) to be removed from the reconstruction.
#' @param partition.number If only a single gene is desired to be reconstructed, the input is the partition number in the selac object.
#'
#' @details
#' Provides marginal probabilities for all nodes across all genes. The function is fairly simple to use as it only requires as input the selac output object and the working directory that the original analysis took place.
GetMarginalAllGenes <- function(selac.obj, aa.optim.input=NULL, fasta.rows.to.keep=NULL, taxon.to.drop, partition.number=NULL) {
    
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
    data.type <- selac.obj$data.type
    parallel.type <- "by.gene"
    n.cores <- NULL
    if(!is.null(partition.number)){
        obj.final <- NULL
        for(partition.index in partition.number:partition.number){
            x <- c(selac.obj$mle.pars[partition.index,])
            gene.tmp <- read.dna(partitions[partition.index], format='fasta')
            if(!is.null(fasta.rows.to.keep)){
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
            }else{
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
            }
            
            if(data.type == "nucleotide"){
                codon.data <- DNAbinToNucleotideNumeric(gene.tmp)
                codon.data <- codon.data[phy$tip.label,]
                codon.freq.by.gene <- selac.obj$empirical.base.freqs[[partition.index]]
                codon.freq.by.aa=NULL
            }else{
                codon.data <- DNAbinToCodonNumeric(gene.tmp)
                codon.data <- codon.data[phy$tip.label,]
                codon.freq.by.aa <- selac.obj$codon.freq.by.aa[[partition.index]]
                codon.freq.by.gene <- selac.obj$codon.freq.by.gene[[partition.index]]
                if(is.null(aa.optim.input)){
                    aa.optim_array <- selac.obj$aa.optim[[partition.index]]
                }
            }
            
            tmp <- GetMarginalGene(pars=x, codon.data=codon.data, phy=phy, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, numcode=numcode, diploid=diploid, taxon.to.drop=taxon.to.drop, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, nuc.model=nuc.model, aa.properties=aa.properties, data.type=data.type)
            obj.gene <- tmp
            obj.final <- obj.gene
        }
    }else{
        obj.final <- as.list(1:length(partitions))
        for(partition.index in 1:length(partitions)){
            x <- c(selac.obj$mle.pars[partition.index,])
            gene.tmp <- read.dna(partitions[partition.index], format='fasta')
            if(!is.null(fasta.rows.to.keep)){
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
            }else{
                gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
            }
            codon.data <- DNAbinToCodonNumeric(gene.tmp)
            codon.data <- codon.data[phy$tip.label,]
            
            if(data.type == "nucleotide"){
                codon.freq.by.gene=selac.obj$empirical.base.freqs[[partition.index]]
                codon.freq.by.aa=NULL
            }else{
                codon.freq.by.aa=selac.obj$codon.freq.by.aa[[partition.index]]
                codon.freq.by.gene=selac.obj$codon.freq.by.gene[[partition.index]]
            }
            
            if(is.null(aa.optim.input)){
                aa.optim_array <- selac.obj$aa.optim[[partition.index]]
            }
            tmp = GetMarginalGene(pars=x, codon.data=codon.data, phy=phy, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, numcode=numcode, diploid=diploid, taxon.to.drop=taxon.to.drop, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, nuc.model=nuc.model, aa.properties=aa.properties, data.type=data.type)
            obj.gene <- tmp
            obj.final[[partition.index]] <- obj.gene
        }
    }
    return(obj.final)
}
