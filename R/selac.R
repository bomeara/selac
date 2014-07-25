#use seqinr coding of nucleotides: see ?n2s: 0 -> "a", 1 -> "c", 2 -> "g", 3 -> "t"


#use corHMM rayDISC to calculate likelihoods

#' Create a distance of physiochemical distances between pairs of amino acids
#' 
#' This is based on the work of Grantham (1974), but allows user specification
#' of the weights
#'
#' @param alpha is the weight on composition (c), the atomic weight ratio of non carbon elements in end groups or rings to carbons in the side chain
#' @param beta is the weight on polarity (p)
#' @param gamma is the weight on molecular volume (v)
#' @param aa.properties is matrix of composition, polarity, and volume parameters for amino acids (three letter codes as rownames); Grantham's is used if not user-supplied
#' @param normalize determines whether the distance matrix is normalized so that the mean distance is 1
#' @return aa.distances, the symmetric matrix of physiochemical distances between amino acids
CreateAADistanceMatrix <- function(alpha=1.829272, beta=0.101799, gamma=0.0003990333, aa.properties=NULL, normalize=TRUE) {
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
  n.states <- dim(aa.properties)[1]
  if(n.states != 20) {
  	warning(paste("aa.properties given use", n.states, "states, normally there are 20 amino acids"))
  }
  aa.distances <- matrix(0,nrow=n.states,ncol=n.states)
  for (i in sequence(n.states)) {
    for (j in sequence(n.states)) {
      aa.distances[i, j] <- (alpha*(aa.properties[i,1] - aa.properties[j,1])^2 + beta*(aa.properties[i,2]-aa.properties[j,2])^2+gamma*(aa.properties[i,3]-aa.properties[j,3])^2)^(1/2)
    }
  }
  if(normalize) {
  	aa.distances <- aa.distances / (sum(aa.distances) / (n.states*n.states - n.states)) #normalize so mean is 1 across the non-diagonal entries	
  }	
  rownames(aa.distances) <- rownames(aa.properties)
  colnames(aa.distances) <- rownames(aa.properties)
  return(aa.distances)
}


#'Create a nucleotide instantaneous mutation matrix, rows=from, cols=to, nuc order as in seqinr coding
#'
#' @param rates are the rates for the given model
#' @param model is the model. Choices are JC, HKY, and GTR
#' @return nuc.mutation.rates matrix of instantaneous mutation rates for nucleotides
CreateNucleotideMutationMatrix <- function(rates=c(1), model="JC") {
  if(model == "JC") {
    nuc.mutation.rates <- matrix(data=rates[1], nrow=4, ncol=4)
    rownames(nuc.mutation.rates) <- n2s(0:3)
    colnames(nuc.mutation.rates) <- n2s(0:3)
    diag(nuc.mutation.rates) <- 0
    diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
    return(nuc.mutation.rates)
  } else {
  	stop("JC is all that is present so far")	
  }
}

#'Create a codon mutation model, rows=from, cols=to, nuc order as in seqinr coding
#'
#' Note it does not assume a symmetric nucleotide matrix nor a symmetric codon mutation rate matrix
#' Stop codons are included
#' @param nuc.mutation.rates are the rates from (row) to (col) based on the nucleotide model
#' @return codon.mutation.rates matrix of instantaneous mutation rates for codons
CreateCodonMutationMatrix <- function(nuc.mutation.rates) {
  codon.sets <- expand.grid(0:3, 0:3, 0:3)
  codon.sets <- data.frame(first=codon.sets[,3], second=codon.sets[,2], third=codon.sets[,1]) #reordering to group similar codons
  n.codons <- dim(codon.sets)[1]
  codon.mutation.rates <- matrix(data=0, nrow=n.codons, ncol=n.codons)
  rownames(codon.mutation.rates) <- rep("",n.codons)
  colnames(codon.mutation.rates) <- rep("",n.codons)
  for (i in sequence(n.codons)) {
    for (j in sequence(n.codons)) {
      if(sum(codon.sets[i,] == codon.sets[j,])==2) { #means that two of the bases match
      	mismatch.position <- which(codon.sets[i,] != codon.sets[j,])
      	codon.mutation.rates[i,j] <- nuc.mutation.rates[1+codon.sets[i,mismatch.position], 1+codon.sets[j, mismatch.position]] #nucs numbered from 0:3, rows are 1:4, thus the add 1
      }
    }
    codon.name <- paste(n2s(as.numeric(codon.sets[i,])), collapse="")
    rownames(codon.mutation.rates)[i] <- codon.name
    colnames(codon.mutation.rates)[i] <- codon.name

  }
  diag(codon.mutation.rates) <- -rowSums(codon.mutation.rates)
  return(codon.mutation.rates)
}

CodonNumericToString <- function(x) { #remember that codon numbers start at 1
	return(n2s(x, levels=words(length=3), base4=FALSE))
}

CodonStringToNumeric <- function(x) { #remember that codon numbers start at 1
	return(which(words(length=3)==x))
}


#Problem with this is that it requires an assumption about frequencies. Flat
#  amino acid frequencies != flat codon frequences != frequencies of codons at
#  equilibrium given nucleotide model for many models. So, let's just go
#  directly to the codon model
#'Create an amino acid mutation model, rows=from, cols=to
#'
#' Note it does not assume a symmetric nucleotide matrix nor a symmetric codon mutation rate matrix
#' @param codon.mutation.rates are the rates from (row) to (col) based on the nucleotide model
#' @param codon.freqs are the frequencies of each codon. By default, flat prior
#' @param numcode is the genetic code as in seqinr, set by default to standard
#' @return aa.mutation.rates matrix of instantaneous mutation rates for aa
#CreateAAMutationMatrix <- function(codon.mutation.rates, codon.freqs = rep(1, 64), numcode=1, include.stop.codon=FALSE) {
#	TranslateCodon <- function(codon, numcode=1) {
#		return(translate(s2c(codon), numcode=numcode))
#	}
#	aa.by.codon <- sapply(rownames(codon.mutation.rates), TranslateCodon, numcode=numcode)
#	codon.mutation.rates.scaled <- codon.mutation.rates*codon.freqs
#	aa.ordered <- unique(aa.by.codon)
#	n.aa <- length(aa.ordered)
#	n.codon <- dim(codon.mutation.rates)[1]
#	aa.mutation.rates <- matrix(data=0, nrow=n.aa, ncol=n.aa)
#	rownames(aa.mutation.rates) <- aa.ordered
#	colnames(aa.mutation.rates) <- aa.ordered
#	for (i in sequence(n.codon)) {
#     for (j in sequence(n.codon)) {
#      	if (i != j) {
#          from.aa <- aa.by.codon[i]
#          to.aa <- aa.by.codon[j]
#          from.location <- which(aa.ordered==from.aa)
#          to.location <- which(aa.ordered==to.aa)
#          aa.mutation.rates[from.location, to.location] <- aa.mutation.rates[from.location, to.location] + codon.mutation.rates[i, j]
#        }
#      }
#    }
#    if(!include.stop.codon) {
#    	aa.mutation.rates<-aa.mutation.rates[which(rownames(aa.mutation.rates)!="*"),]
#    	aa.mutation.rates<-aa.mutation.rates[,which(colnames(aa.mutation.rates)!="*")]
#    }
#    diag(aa.mutation.rates) <- 0
#    diag(aa.mutation.rates) <- -rowSums(aa.mutation.rates)
#    return(aa.mutation.rates)
#}

CompareVectors <- function(cd1,cd2){
  cmp = (cd1 != cd2)
  num = sum(cmp)
  pos = which(cmp==TRUE)
  return(list(num=num,pos=pos))
}


GetPairwiseProteinFixationProbabilityArbitraryLength <-function(protein1, protein2, protein_op, s, aa.distances, C=2,Phi=0.5,q=4e-7,Ne=5e6){
  d1 <- GetProteinProteinDistance(protein1,protein_op,aa.distances)
  d2 <- GetProteinProteinDistance(protein2,protein_op,aa.distances)
  if(length(d1)!=length(d2)) #throw error if length of proteins are not the same
    stop("error: 2 proteins are of different lengths!")
  if(length(d1)==1){ #only one amino acid
    return(GetPairwiseProteinFixationProbabilitySingleSite(d1,d2,s,C=C,Phi=Phi,q=q,Ne=Ne))
  }
  else{
    if((length(s)==1)&&(length(d1)!=1)) #if s is given as a scalar, then treat it to be the same across all sites
      s <- rep(s,length(d1))
    l = length(d1)
    cmp = CompareVectors(d1,d2)
    if(cmp$num > 1) return(0) #more than 1 position differ
    else if((cmp$num ==0)) return(1/(2*Ne)) #same fitness/functionality
    else{ #exactly 1 position differs
      pos = cmp$pos 
      return(GetPairwiseProteinFixationProbabilitySingleSite(d1[pos],d2[pos],s[pos],C=C,Phi=Phi,q=q,Ne=Ne))
    }
  }
}

GetPairwiseProteinFixationProbabilitySingleSite <- function(d1,d2,s,C=2,Phi=0.5,q=4e-7,Ne=5e6){
  if((d1==d2)||(s==0)) #When the fitnesses are the same, neutral case, pure drift
    return(1/(2*Ne))
  else{
      fit_ratio <- exp(-C*Phi*q*s*(d1-d2)) #f1/f2
     if(fit_ratio==Inf) #1 is much better than 2 (the mutant)
       return(0)
     else if(fit_ratio==1)
       return(1/(2*Ne))
     else
    return((1-fit_ratio)/(1-fit_ratio^(2*Ne)))
  }
}

TranslateCodon <- function(codon.string, numcode=1) {
  return(translate(s2c(codon.string), numcode=numcode))
}

GetProteinProteinDistance <- function(protein1, protein2,aa.distances){
  if(length(protein1)!=length(protein2)) #throw error if length of proteins are not the same
    stop("error: 2 proteins are of different lengths!")
  site_d <- function(k){
    if((is.na(protein1[k])) & (!is.na(protein2[k]))){
      return(mean(aa.distances[,protein2[k]]))
    }
    else if((is.na(protein2[k])) & (!is.na(protein1[k]))){
      return(mean(aa.distances[protein1[k],]))
    }
    else if((is.na(protein2[k])) & (!is.na(protein1[k]))){
      return(mean(aa.distances))
    }
    else
      return(aa.distances[protein1[k],protein2[k]])
  }
  d <- sapply(c(1:length(protein1)),site_d,simplify=TRUE) 
  return(d)
}

CreateCodonFixationRateMatrix <- function(aa_op, s, aa.distances, C=2, Phi=0.5,q=4e-7,Ne=5e6, include.stop.codon=FALSE, numcode=1){
  codon.sets <- expand.grid(0:3, 0:3, 0:3)
  codon.sets <- data.frame(first=codon.sets[,3], second=codon.sets[,2], third=codon.sets[,1]) #reordering to group similar codons
  n.codons <- dim(codon.sets)[1]
  codon.fixation.rates <- matrix(data=0, nrow=n.codons, ncol=n.codons)
  codon.names <- rep("", n.codons)
  for (i in sequence(n.codons)) {
  	codon.names[i] <- paste(n2s(as.numeric(codon.sets[i,])), collapse="")
  }
  rownames(codon.fixation.rates) <- codon.names
  colnames(codon.fixation.rates) <- codon.names
  codon.aa <- sapply(codon.names, TranslateCodon, numcode=numcode)
  
  for (i in sequence(n.codons)) {
    for (j in sequence(n.codons)) {
    	if(sum(codon.sets[i,] == codon.sets[j,])>=2) { #match at two or three sites of three
    		aa1 <- TranslateCodon(paste(n2s(as.numeric(codon.sets[i,])), collapse=""), numcode=numcode)
    		aa2 <- TranslateCodon(paste(n2s(as.numeric(codon.sets[j,])), collapse=""), numcode=numcode)
    		if(aa1!="*" && aa2!="*") { #says we cannot mutate to stop codons
    			d1 <- GetProteinProteinDistance(protein1=aa1, protein2=aa_op, aa.distances=aa.distances)
    			d2 <- GetProteinProteinDistance(protein1=aa2, protein2=aa_op, aa.distances=aa.distances)
    			codon.fixation.rates[i,j] <- GetPairwiseProteinFixationProbabilitySingleSite(d1, d2, s=s, C=C, Phi=Phi, q=q, Ne=Ne)
    		} else {
    			if(s==0) { #handles stop codon case where neutral, so could possibly go into and out of stop codons
    				codon.fixation.rates[i,j] <- 1/(2*Ne)
    			}
    		}
    	}
    }
  }
  if(!include.stop.codon) {
    codon.fixation.rates<-codon.fixation.rates[which(codon.aa!="*"),]
    codon.fixation.rates<-codon.fixation.rates[,which(codon.aa!="*")]
  }
  diag(codon.fixation.rates) <- 0 #b/c we don't want these included when calculating diag
  diag(codon.fixation.rates) <- -rowSums(codon.fixation.rates)
  return(codon.fixation.rates)
}

CreateAAFixationMatrix <- function(aa_op,s,aa.distances,C=2, Phi=0.5,q=4e-7,Ne=5e6){
  m = 20
  mat <- matrix(0,nrow=m,ncol=m)#set diagonal entries to be 0 at first
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      mat[i,j] <- GetPairwiseProteinFixationProbabilityArbitraryLength(i,j,aa_op,s,aa.distances,C,Phi,q,Ne) #fixation prob -> transition rate
      mat[j,i] <- GetPairwiseProteinFixationProbabilityArbitraryLength(j,i,aa_op,s,aa.distances,C,Phi,q,Ne) #symmetric entry (not the same rate!)
    }#end for j
  }#end for i
  return(mat)
}


#'Get likelihood for a given AA site given tree and Q matrix, assuming the optimal AA is known
#' @param aa.data is data in corHMM format (one column species, one column state)
#' @param phy is a phylo object
#' @param Q_aa is the transition matrix Q for a given optimal aa at this site
#' @param charnum is which character to use in the aa.data matrix
#' @param root.p is the root frequency: NULL if equilibrium, maddfitz, or a vector
#' @param return.all, if TRUE, allows return of everything from rayDISC rather than just the log likelihood
GetLikelihoodSAC_AAForSingleCharGivenOptimum <- function(aa.data, phy, Q_aa, charnum=1, root.p=NULL, return.all=FALSE) {
	result <- rayDISC(phy=phy, data=aa.data, ntraits=1, charnum=charnum, p=Q_aa, root.p=root.p)
	ifelse(return.all, return(result), return(result$loglik))
}

#'Get likelihood for a given codon site given tree and Q matrix, assuming the optimal AA is known
#' @param codon.data is data in corHMM format (one column species, one column state, states from CodonStringToNumeric)
#' @param phy is a phylo object
#' @param Q_codon is the transition matrix Q for a given optimal aa at this site
#' @param charnum is which character to use in the aa.data matrix
#' @param root.p is the root frequency: NULL if equilibrium, maddfitz, or a vector
#' @param return.all, if TRUE, allows return of everything from rayDISC rather than just the log likelihood
GetLikelihoodSAC_CodonForSingleCharGivenOptimum <- function(codon.data, phy, Q_codon, charnum=1, root.p=NULL, return.all=FALSE) {
	result <- rayDISC(phy=phy, data=codon.data, ntraits=1, charnum=charnum, p=Q_codon, root.p=root.p)
	ifelse(return.all, return(result), return(result$loglik))
}

PlotBubbleMatrix <- function(x, main="", special=Inf, cex=1){
diag(x) <- 0
	x<-x/max(x)

  plot(x=range(.5,.5+dim(x)[2]),y=-range(.5, .5+dim(x)[1]), xlab="", ylab="", type="n", main=main,xaxt='n',yaxt='n', asp=1,bty="n")
  	axis(side=2, at=-sequence(dim(x)[1]), labels=rownames(x), las=2, cex=cex)
	axis(side=3, at=sequence(dim(x)[2]), labels=colnames(x), las=2, cex=cex)

  for (i in sequence(dim(x)[2])) {
  	for (j in sequence(dim(x)[1])) {
  		bg="black"
  		if(i%in% special || j %in% special) {
  			bg="red"
  		}
  		if(x[j,i]>0) {
  			symbols(x=i, y=-j, circles=sqrt(x[j,i])/(2.1*sqrt(max(x))), inches=FALSE, add=TRUE, fg=bg, bg=bg)
  		}
  	}	
  }
}

GetGainLossRatios <- function(x) {
  diag(x) <- 0
  x<-x/max(x)
  ratio.matrix <- x*0
  for (i in sequence(dim(x)[2])) {
  	for (j in sequence(dim(x)[1])) {
  		if(i>j) {
  		 	gain.rate <- x[j, i]
  			loss.rate <- x[i, j]
  			ratio<-gain.rate/(loss.rate+gain.rate)-0.5
  			if(is.na(ratio)) {
  				ratio <- 0
  			}
  			ratio.matrix[j,i]<-ratio
  		}
  	}
  }
  return(ratio.matrix)
}

#blue is proportional increase (largest circle means gain rate is positive, loss rate is zero)
#red is proportional decrease (circle of 0.5 means loss rate is twice that of gain rate)
PlotBubbleRatio <- function(x, main="", cex=1){
  ratio.matrix<-GetGainLossRatios(x)
    plot(x=range(.5,.5+dim(x)[2]),y=-range(.5, .5+dim(x)[1]), xlab="", ylab="", type="n", main=main,xaxt='n',yaxt='n', asp=1,bty="n")
  	axis(side=2, at=-sequence(dim(x)[1]), labels=rownames(x), las=2, cex=cex)
	axis(side=3, at=sequence(dim(x)[2]), labels=colnames(x), las=2, cex=cex)
  for (i in sequence(dim(x)[2])) {
  	for (j in sequence(dim(x)[1])) {
  		if(ratio.matrix[j,i]!=0) {
  			gain.rate <- x[j, i]
  			loss.rate <- x[i, j]
  			ratio<-ratio.matrix[j,i]
  			bg="blue"
  			if (ratio < 0) {
  				bg="red"
  			}
  			symbols(x=i, y=-j, circles=0.5*sqrt(abs(ratio))/sqrt(max(abs(ratio.matrix))), inches=FALSE, add=TRUE, fg=bg, bg=bg)
  		}
  	}	
  }
}

PlotTransitionNetwork <- function(x, main="") {
	require(igraph)
	diag(x) <- 0
	x<-x/max(x)
	g <- graph.adjacency(x, weighted=TRUE, mode="directed")
	g.layout <- layout.fruchterman.reingold(g)
	plot(g, layout=g.layout, edge.width=10*get.edge.attribute(g, "weight"), edge.curved=TRUE)
}

DNAbinToCodonNumeric <- function(x, frame=0, corHMM.format=TRUE) {
	bound.characters <- sapply(as.character(x), paste, collapse="")
	#following fn is derived from code for uco in seqinr
	SplitToCodons <- function(seq.string, frame) {
		seq.string<-strsplit(seq.string, split="")[[1]]
		if (any(seq.string %in% LETTERS)) {
    	    seq.string <- tolower(seq.string)
   		 }
   		 return(sapply(splitseq(seq = seq.string, frame = frame, word = 3), CodonStringToNumeric))
   	}
	split.characters <- t(sapply(bound.characters, SplitToCodons, frame=frame))
	colnames(split.characters) <- sequence(dim(split.characters)[2])
	if(corHMM.format) {
		split.characters<-cbind(data.frame(Taxa=rownames(split.characters)), data.frame(split.characters))
	}
	return(split.characters)
}

