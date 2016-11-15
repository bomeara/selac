
library(selac)
library(seqinr)
library(parallel)
library(expm)

.codon.name <- c("aaa" ,"aac" ,"aag" ,"aat" ,"aca" ,"acc" ,"acg" ,"act" ,"aga" ,"agc" ,"agg" ,"agt" ,"ata" ,"atc" ,"atg" ,"att" ,"caa" ,"cac" ,"cag" ,"cat", "cca" ,"ccc",
"ccg" ,"cct" ,"cga" ,"cgc" ,"cgg" ,"cgt" ,"cta" ,"ctc" ,"ctg" ,"ctt" ,"gaa" ,"gac" ,"gag" ,"gat" ,"gca" ,"gcc" ,"gcg" ,"gct" ,"gga" ,"ggc", "ggg" ,"ggt",
"gta" ,"gtc" ,"gtg" ,"gtt" ,"taa" ,"tac" ,"tag" ,"tat" ,"tca" ,"tcc" ,"tcg" ,"tct" ,"tga" ,"tgc" ,"tgg" ,"tgt" ,"tta" ,"ttc" ,"ttg" ,"ttt")


.unique.aa <- c("K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E", "D", "A", "G", "V", "*", "Y", "C", "W", "F")


TranslateCodon <- function(codon.string, numcode) {
    return(translate(s2c(codon.string), numcode=numcode))
}


.aa.translation <- list(sapply(.codon.name, TranslateCodon, numcode=1),
sapply(.codon.name, TranslateCodon, numcode=2),
sapply(.codon.name, TranslateCodon, numcode=3),
sapply(.codon.name, TranslateCodon, numcode=4),
sapply(.codon.name, TranslateCodon, numcode=5),
sapply(.codon.name, TranslateCodon, numcode=6),
sapply(.codon.name, TranslateCodon, numcode=9),
sapply(.codon.name, TranslateCodon, numcode=10),
sapply(.codon.name, TranslateCodon, numcode=11),
sapply(.codon.name, TranslateCodon, numcode=12),
sapply(.codon.name, TranslateCodon, numcode=13),
sapply(.codon.name, TranslateCodon, numcode=14),
sapply(.codon.name, TranslateCodon, numcode=15),
sapply(.codon.name, TranslateCodon, numcode=16),
sapply(.codon.name, TranslateCodon, numcode=21),
sapply(.codon.name, TranslateCodon, numcode=22),
sapply(.codon.name, TranslateCodon, numcode=23))


GetTipMarginalSingleSite <- function(charnum=1, codon.data, phy, Q, root.p=NULL, taxon.to.drop=NULL, state=1){
    
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
        liks[taxon.to.drop,] <- 0
        liks[taxon.to.drop, state] <- 1
        if(nl > 4){
            liks[taxon.to.drop,c(49, 51, 57)] <- 0
        }
    }
    
    anc <- unique(phy$edge[,1])
    
    TIPS <- 1:nb.tip
    comp <- numeric(nb.tip + nb.node)
    
    if(any(root.p < 0) | any(is.na(root.p))){
        return(1000000)
    }
    
    for (i  in seq(from = 1, length.out = nb.node)) {
        #the ancestral node at row i is called focal
        focal <- anc[i]
        #Get descendant information of focal
        desRows<-which(phy$edge[,1]==focal)
        desNodes<-phy$edge[desRows,2]
        v <- 1
        
        for (desIndex in sequence(length(desRows))){
            v<-v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
        }
        comp[focal] <- sum(v)
        liks[focal, ] <- v/comp[focal]
    }
    
    
    #Specifies the root:
    root <- nb.tip + 1L
    #If any of the logs have NAs restart search:
    if(is.nan(sum(log(comp[-TIPS]))) || is.na(sum(log(comp[-TIPS])))){
        return(1000000)
    }
    else{
        loglik<- (sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
        if(is.infinite(loglik)){return(1000000)}
    }
    return(loglik)
}


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
    
    #Now get the states for the tips (will do, not available for general use):
    likelihoods.vector <- c()
    
    if(nl > 4){
        for (k in 1:64){
            if(k == 49 | k == 51 | k==57){
                likelihoods.vector <- c(likelihoods.vector, -1000000)
            }else{
                #the ancestral node at row i is called focal
                tmp <- GetTipMarginalSingleSite(charnum=charnum, codon.data=codon.data, phy=phy, Q=Q, root.p=root.p, taxon.to.drop=taxon.to.drop, state=k)
                likelihoods.vector <- c(likelihoods.vector, tmp)
            }
        }
    }else{
        for (k in 1:4){
            #the ancestral node at row i is called focal
            tmp <- GetTipMarginalSingleSite(charnum=charnum, codon.data=codon.data, phy=phy, Q=Q, root.p=root.p, taxon.to.drop=taxon.to.drop, state=k)
            likelihoods.vector <- c(likelihoods.vector, tmp)
        }
    }
    
    best.probs <- max(likelihoods.vector)
    marginal.probs.tmp <- likelihoods.vector - best.probs
    marginal.probs <- exp(marginal.probs.tmp) / sum(exp(marginal.probs.tmp))
    liks.final[taxon.to.drop,] <- marginal.probs
    
    #Just add in the marginal at the root calculated on the original downpass or if supplied by the user:
    liks.final[root,] <- liks.down[root,] * root.p
    root.final <- liks.down[root,] * root.p
    comproot <- sum(root.final)
    liks.final[root,] <- root.final/comproot
    return(liks.final)
}


GetMarginalGene <- function(pars, codon.data, phy, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, aa.optim_array, numcode, diploid, taxon.to.drop=NULL, include.gamma=FALSE, gamma.type="quadrature", ncats=4, k.levels=0, nuc.model="UNREST", aa.properties=NULL, data.type="codon"){
    
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
        nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(transition.rates, model=nuc.model)
        if(include.gamma==TRUE){
            if(gamma.type == "median"){
                rates.k <- selac:::DiscreteGamma(shape=shape, ncats=ncats)
                weights.k <- rep(1/ncats, ncats)
            }
            if(gamma.type == "quadrature"){
                rates.and.weights <- selac:::LaguerreQuad(shape=shape, ncats=ncats)
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
        codon.index.matrix = selac:::CreateCodonMutationMatrixIndex()
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
                nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
            }
            if(nuc.model == "GTR") {
                base.freqs=c(pars[4:6], 1-sum(pars[4:6]))
                nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(x[9:length(pars)], model=nuc.model, base.freqs=base.freqs)
            }
            if(nuc.model == "UNREST") {
                nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(x[6:length(pars)], model=nuc.model)
            }
        }else{
            if(nuc.model == "JC") {
                base.freqs=c(pars[4:6], 1-sum(pars[4:6]))
                nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
            }
            if(nuc.model == "GTR") {
                base.freqs=c(pars[4:6], 1-sum(pars[4:6]))
                nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(pars[7:length(pars)], model=nuc.model, base.freqs=base.freqs)
            }
            if(nuc.model == "UNREST") {
                nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(pars[4:length(pars)], model=nuc.model)
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
                rates.k <- selac:::DiscreteGamma(shape=shape, ncats=ncats)
                weights.k <- rep(1/ncats, ncats)
            }
            if(gamma.type == "quadrature"){
                rates.and.weights <- selac:::LaguerreQuad(shape=shape, ncats=ncats)
                rates.k <- rates.and.weights[1:ncats]
                weights.k <- rates.and.weights[(ncats+1):(ncats*2)]
            }
            for (i in sequence(nsites)) {
                marginal.recon_ncats_array <- array(1, dim=c(Ntip(phy)+Nnode(phy), 64, ncats))
                for(k.cat in sequence(ncats)){
                    if(k.levels > 0){
                        aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
                    }else{
                        aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
                    }
                    Q_codon_array <- selac:::FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[k.cat], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999)
                    
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
                aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
            }else{
                aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
            }
            Q_codon_array <- selac:::FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999)
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


GetMarginalAllGenes <- function(selac.obj, aa.optim.input=NULL, fasta.rows.to.keep=NULL, taxon.to.drop, partition.number, data.type="codon") {
    
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
                codon.data <- selac:::DNAbinToNucleotideNumeric(gene.tmp)
                codon.data <- codon.data[phy$tip.label,]
                codon.freq.by.gene=selac.obj$empirical.base.freqs[[partition.index]]
                codon.freq.by.aa=NULL
            }else{
                codon.data <- selac:::DNAbinToCodonNumeric(gene.tmp)
                codon.data <- codon.data[phy$tip.label,]
                codon.freq.by.aa=selac.obj$codon.freq.by.aa[[partition.index]]
                codon.freq.by.gene=selac.obj$codon.freq.by.gene[[partition.index]]
                if(is.null(aa.optim.input)){
                    aa.optim_array <- selac.obj$aa.optim[[partition.index]]
                }
            }
            
            tmp = GetMarginalGene(pars=x, codon.data=codon.data, phy=phy, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, numcode=numcode, diploid=diploid, taxon.to.drop=taxon.to.drop, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, nuc.model=nuc.model, aa.properties=aa.properties, data.type=data.type)
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
            codon.data <- selac:::DNAbinToCodonNumeric(gene.tmp)
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


GetFunctionalityMod <- function(gene.length, aa.data, optimal.aa, alpha, beta, gamma, aa.properties=NULL){
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
    aa.distances <- c()
    #Note using only the second row, because we are comparing empirical S. cervisae rates:
    for(site.index in 1:gene.length){
        if(aa.data[site.index]=="*"){
            gene.length = gene.length - 1
        }else{
            if(aa.data[site.index]!="NA"){
                aa.distances <- c(aa.distances, (1+((alpha*(aa.properties[aa.data[site.index],1] - aa.properties[optimal.aa[site.index],1])^2 + beta*(aa.properties[aa.data[site.index],2]-aa.properties[optimal.aa[site.index],2])^2+gamma*(aa.properties[aa.data[site.index],3]-aa.properties[optimal.aa[site.index],3])^2)^(1/2))))
            }else{
                aa.distances <- c(aa.distances, 0)
            }
        }
    }
    functionality = 1/((1/gene.length) * sum(aa.distances))
    return(functionality)
}


GetAdequateSelac <- function(selac.obj, aa.optim.input=NULL, fasta.rows.to.keep=NULL, taxon.to.drop=4, partition.number=55, numcode=1){
    tot.interval <- selac.obj$phy$edge.length[selac.obj$phy$edge[,2]==taxon.to.drop]
    prop.interval <- seq(0,1 , by=0.05)
    time.interval <- tot.interval * prop.interval
    functionality.dropped.taxon.intervals <- c()
    for(interval.index in 1:length(time.interval)){
        selac.obj$phy$edge.length[selac.obj$phy$edge[,2]==taxon.to.drop] <- time.interval[interval.index]
        marginals.across.tree <- GetMarginalAllGenes(selac.obj=selac.obj, aa.optim.input=NULL, fasta.rows.to.keep=NULL, taxon.to.drop=taxon.to.drop, partition.number=partition.number)
        reconstructed.sequence <- c()
        for(site.index in 1:selac.obj$nsites[partition.number]){
            reconstructed.sequence <- c(reconstructed.sequence, .aa.translation[[numcode]][sample(1:64, 1, prob=marginals.across.tree[taxon.to.drop,,site.index])])
        }
        reconstructed.sequence <- unname(reconstructed.sequence)
        functionality.dropped.taxon <- GetFunctionalityMod(gene.length=selac.obj$nsites[partition.number], aa.data=reconstructed.sequence, optimal.aa=selac.obj$aa.optim[[partition.number]], alpha=selac.obj$mle.pars[1,2], beta=selac.obj$mle.pars[1,3], gamma=selac.obj$volume.fixed.value, aa.properties=selac.obj$aa.properties)
        functionality.dropped.taxon.intervals <- c(functionality.dropped.taxon.intervals, functionality.dropped.taxon)
    }
    final.res <- cbind(time.interval, functionality.dropped.taxon.intervals)
    return(final.res)
}


GetAdequateGTR <- function(selac.obj, selac.obj2, fasta.rows.to.keep=NULL, taxon.to.drop=4, partition.number=55, numcode=1){
    tot.interval <- selac.obj$phy$edge.length[selac.obj$phy$edge[,2]==taxon.to.drop]
    prop.interval <- seq(0,1 , by=0.05)
    time.interval <- tot.interval * prop.interval
    functionality.dropped.taxon.intervals <- c()
    stop.codons.simulated <- c()
    for(interval.index in 1:length(time.interval)){
        selac.obj$phy$edge.length[selac.obj$phy$edge[,2]==taxon.to.drop] <- time.interval[interval.index]
        marginals.across.tree <- GetMarginalAllGenes(selac.obj=selac.obj, aa.optim.input=NULL, fasta.rows.to.keep=NULL, taxon.to.drop=taxon.to.drop, partition.number=partition.number, data.type="nucleotide")
        reconstructed.sequence <- c()
        for(site.index in 1:selac.obj$nsites[partition.number]){
            reconstructed.sequence <- c(reconstructed.sequence, n2s(0:3)[sample(1:4, 1, prob=marginals.across.tree[taxon.to.drop,,site.index])])
        }
        reconstructed.sequence.codon <- sapply(splitseq(seq = reconstructed.sequence, frame = frame, word = 3), selac:::CodonStringToNumeric)
        reconstructed.sequence <- c()
        for(site.index in 1:length(reconstructed.sequence.codon)){
            reconstructed.sequence <- c(reconstructed.sequence, .aa.translation[[numcode]][reconstructed.sequence.codon[site.index]])
        }
        functionality.dropped.taxon <- GetFunctionalityMod(gene.length=selac.obj2$nsites[partition.number], aa.data=reconstructed.sequence, optimal.aa=selac.obj2$aa.optim[[partition.number]], alpha=selac.obj2$mle.pars[1,2], beta=selac.obj2$mle.pars[1,3], gamma=selac.obj2$volume.fixed.value, aa.properties=selac.obj2$aa.properties)
        functionality.dropped.taxon.intervals <- c(functionality.dropped.taxon.intervals, functionality.dropped.taxon)
        stop.codons.simulated <- c(stop.codons.simulated, length(which(reconstructed.sequence=="*")))
    }
    final.res <- cbind(time.interval, functionality.dropped.taxon.intervals, stop.codons.simulated)
    return(final.res)
}


GetAdequateManyReps <- function(nreps, n.cores, selac.obj, selac.obj2=NULL, aa.optim.input=NULL, fasta.rows.to.keep=NULL, taxon.to.drop=2, partition.number=17, numcode=1, data.type="codon"){
    if(data.type == "nucleotide"){
        RunMany <- function(nrep){
            tmp <- GetAdequateGTR(selac.obj=selac.obj, selac.obj2=selac.obj2, fasta.rows.to.keep=fasta.rows.to.keep, taxon.to.drop=taxon.to.drop, partition.number=partition.number, numcode=numcode)
            return(tmp)
        }
        pp <- mclapply(1:nreps, RunMany, mc.cores=n.cores)
    }else{
        RunMany <- function(nrep){
            tmp <- GetAdequateSelac(selac.obj=selac.obj, fasta.rows.to.keep=fasta.rows.to.keep, taxon.to.drop=taxon.to.drop, partition.number=partition.number, numcode=numcode)
            return(tmp)
        }
        pp <- mclapply(1:nreps, RunMany, mc.cores=n.cores)
    }
    return(pp)
}


## Adequacy of selac + Gamma ##
#load("yeastRokasSelacUNRESTgamma.Rdata")
#pp <- GetAdequateManyReps(100, n.cores=24, selac.obj=result, selac.obj2=NULL, fasta.rows.to.keep=1:7, taxon.to.drop=2, partition.number=17, data.type="codon")
#save(pp, file="selac.selacG.adequacy.Rsave")
###############################

## Adequacy of selac + Gamma ##
#load("yeastRokasSelacUNRESTgamma.Rdata")
#selac.obj2 <- result
#load("yeastRokasGTRgamma.Rdata")
#selac.obj <- result
#selac.obj$gamma.type <- "median"
#pp <- GetAdequateManyReps(100, n.cores=24, selac.obj=selac.obj, selac.obj2=selac.obj2, fasta.rows.to.keep=1:7, taxon.to.drop=2, partition.number=17)
###############################




