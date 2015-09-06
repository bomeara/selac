

######################################################################################################################################
######################################################################################################################################
### SELAC -- SELection on Amino acids and/or Codons 
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu and Brian O

###LOAD REQUIRED PACKAGES -- eventually move to namespace:
library(expm)
library(nnet)
library(nloptr)
library(seqinr)
library(phangorn)
library(parallel)

# Use seqinr coding of nucleotides: see ?n2s: 0 -> "a", 1 -> "c", 2 -> "g", 3 -> "t"

#Numcodes:
# 1 standard
# 2 vertebrate.mitochondrial
# 3 yeast.mitochondrial
# 4 protozoan.mitochondrial+mycoplasma
# 5 invertebrate.mitochondrial
# 6 ciliate+dasycladaceal
# 9 echinoderm+flatworm.mitochondrial
# 10 euplotid
# 11 bacterial+plantplastid
# 12 alternativeyeast
# 13 ascidian.mitochondrial
# 14 alternativeflatworm.mitochondrial
# 15 blepharism
# 16 chlorophycean.mitochondrial
# 21 trematode.mitochondrial
# 22 scenedesmus.mitochondrial
# 23 hraustochytrium.mitochondria

######################################################################################################################################
######################################################################################################################################
### Variation functions used by main function:
######################################################################################################################################
######################################################################################################################################

# Create a distance of physiochemical distances between pairs of amino acids
# 
# This is based on the work of Grantham (1974), but allows user specification
# of the weights
#
# @param alpha is the weight on composition (c), the atomic weight ratio of non carbon elements in end groups or rings to carbons in the side chain
# @param beta is the weight on polarity (p)
# @param gamma is the weight on molecular volume (v)
# @param aa.properties is matrix of composition, polarity, and volume parameters for amino acids (three letter codes as rownames); Granthams is used if not user-supplied
# @param normalize determines whether the distance matrix is normalized so that the mean distance is 1
# @return aa.distances, the symmetric matrix of physiochemical distances between amino acids
CreateAADistanceMatrixOriginal <- function(alpha=1.829272, beta=0.101799, gamma=0.0003990333, aa.properties=NULL, normalize=TRUE) {
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
		warning(paste("aa.properties given", n.states, "states, normally there are 20 amino acids"))
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


CreateAADistanceMatrix <- function(alpha=1.829272, beta=0.101799, gamma=0.0003990333, aa.properties=NULL, normalize=FALSE, poly.params=NULL, k=0) {
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
		warning(paste("aa.properties given", n.states, "states, normally there are 20 amino acids"))
	}
	aa.distances <- matrix(0,nrow=n.states,ncol=n.states)
	if(k > 0){
		for (i in sequence(n.states)) {
			for (j in sequence(n.states)) {
				aa.distances[i, j] <- PolynomialTransform(x=(alpha*(aa.properties[i,1] - aa.properties[j,1])^2 + beta*(aa.properties[i,2]-aa.properties[j,2])^2+gamma*(aa.properties[i,3]-aa.properties[j,3])^2)^(1/2), xi=0, poly.params=poly.params, k=1)
			}
		}
	}else{
		for (i in sequence(n.states)) {
			for (j in sequence(n.states)) {
				aa.distances[i, j] <- PolynomialTransform(x=(alpha*(aa.properties[i,1] - aa.properties[j,1])^2 + beta*(aa.properties[i,2]-aa.properties[j,2])^2+gamma*(aa.properties[i,3]-aa.properties[j,3])^2)^(1/2), xi=0, poly.params=NULL, k=0)
			}
		}		
	}
	if(normalize) {
		aa.distances <- aa.distances / (sum(aa.distances) / (n.states*n.states - n.states)) #normalize so mean is 1 across the non-diagonal entries	
	}
	rownames(aa.distances) <- rownames(aa.properties)
	colnames(aa.distances) <- rownames(aa.properties)
	return(aa.distances)
}


###TEMPORARY### Used to test whether other distances are worth using.
GenerateAAProperties <- function(rows){
	mike.table <- read.delim("table.of.aa.attributes.from.Sharma.et.al.2013.csv", sep=",")
	aa.properties <- structure(unlist(t(mike.table[rows,-1])), .Dim = c(20L, 3L), .Dimnames = list(c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"), c("c", "p", "v")))
	return(aa.properties)
}


# Follows Grantham distance equations for generating starting values
GetAADistanceStartingParameters <- function(aa.properties){
	if(is.null(aa.properties)){
		aa.properties <- structure(c(0, 2.75, 1.38, 0.92, 0, 0.74, 0.58, 0, 0.33, 0, 0, 
									 1.33, 0.39, 0.89, 0.65, 1.42, 0.71, 0, 0.13, 0.2, 8.1, 5.5, 13, 
									 12.3, 5.2, 9, 10.4, 5.2, 11.3, 4.9, 5.7, 11.6, 8, 10.5, 10.5, 
									 9.2, 8.6, 5.9, 5.4, 6.2, 31, 55, 54, 83, 132, 3, 96, 111, 119, 
									 111, 105, 56, 32.5, 85, 124, 32, 61, 84, 170, 136), .Dim = c(20L, 3L), .Dimnames = list(c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"), c("c", "p", "v"))) #properties from Grantham paper		
	}
	average.chemical.distance <- c()
	for(i in 1:3){
		average.chemical.distance <- c(average.chemical.distance, mean(dist(aa.properties[,i])))	
	}
	weighting.factors <- (1 / average.chemical.distance)^2
	return(weighting.factors)
}


#pg. 24 "Eq. 3.2.7: m2k_1 = bk,0 + bk,1x + bk,2x^2 + bk,3x^3 + ... + ak,2kx^2k+1, where bk,0=xi, bk,j=ak,(j-1)/j, j = 1,2, ... , 2k+1."
PolynomialTransform <- function(x, xi, poly.params, k){
    if(k == 0){
        mk2_1 <- xi + x
    }
    if(k == 1){
		mk2_1 <- xi + x + poly.params[1]*(x^2) + poly.params[2]*(x^3)
	}
	
	##Not necessary yet:
	#if(k == 2){
	#	mk2_1 <- xi + (coef.vec[2-1,]/2)*x + (coef.vec[3-1,]/3)*(x^2) + (coef.vec[4-1]/4)*(x^3) + (coef.vec[5-1]/5)*(x^4) + (coef.vec[6-1]/6)*(x^5)
	#}
	#if(k == 3){
	#	mk2_1 <- xi + (coef.vec[2-1,]/2)*x + (coef.vec[3-1,]/3)*(x^2) + (coef.vec[4-1]/4)*(x^3) + (coef.vec[5-1]/5)*(x^4) + (coef.vec[6-1]/6)*(x^5) + (coef.vec[7-1]/7)*(x^6) + (coef.vec[8-1]/8)*(x^7)
	#}
	####################
    
	return(mk2_1)
}


#Using matrix algebra from pg. 28.
CalculatePolynomialCoefficients <- function(alpha.poly, beta.poly, k){
	#coef.mat is the matrix of coefficients:
	coef.vec <- 1
	#k.set is a sequence of k ending with the max k which is specified at the function call:
	k.set = 1:k
	#phi is a transformation of alpha and beta:
	phi.poly = alpha.poly^2 + beta.poly
	#Matrix multiplication order does not matter for answer, but order matters for efficiency:
	for(k.index in 1:k){
		coef.vec <- CreatePolynomialMatrix(alpha.poly[k.index], phi.poly[k.index], k.set[k.index]) %*% coef.vec
	}
	return(coef.vec)
}


#See pg. 25-28 for structure of matrices. Easy.
CreatePolynomialMatrix <- function(alpha, phi, k){
	T.mat.k <- matrix(0, 2*k+1, 2*k-1)
	for(column.index in 1:(2*k-1)){
		if(column.index == 1){
			T.mat.k[1, column.index] = 1
			T.mat.k[2, column.index] = -2 * alpha
			T.mat.k[3, column.index] = phi
		}else{
			T.mat.k[column.index, column.index] = 1
			T.mat.k[column.index+1, column.index] = -2 * alpha
			T.mat.k[column.index+2, column.index] = phi
		}
	}
	return(T.mat.k)	
}


#Create a nucleotide instantaneous mutation matrix, rows=from, cols=to, nuc order as in seqinr coding
#
# @param rates are the rates for the given model
# @param model is the model. Choices are JC, HKY, and GTR
# @return nuc.mutation.rates matrix of instantaneous mutation rates for nucleotides
CreateNucleotideMutationMatrix <- function(rates, model="JC", base.freqs=NULL) {
	if(model == "JC") {
		nuc.mutation.rates <- matrix(data=rates[1], nrow=4, ncol=4)
		rownames(nuc.mutation.rates) <- n2s(0:3)
		colnames(nuc.mutation.rates) <- n2s(0:3)
		diag(nuc.mutation.rates) <- 0
		diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
		if(!is.null(base.freqs)){
			diag(nuc.mutation.rates) = 0
			nuc.mutation.rates = t(nuc.mutation.rates * base.freqs)	
			diag(nuc.mutation.rates) = -rowSums(nuc.mutation.rates)
		}
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
		rownames(nuc.mutation.rates) <- n2s(0:3)
		colnames(nuc.mutation.rates) <- n2s(0:3)
		nuc.mutation.rates[4,3] = nuc.mutation.rates[3,4] = 1
		diag(nuc.mutation.rates) <- 0
		diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
		if(!is.null(base.freqs)){
			diag(nuc.mutation.rates) = 0
			nuc.mutation.rates = t(nuc.mutation.rates * base.freqs)	
			diag(nuc.mutation.rates) = -rowSums(nuc.mutation.rates)
		}
		return(nuc.mutation.rates)			
	} 
	if(model == "UNREST") {
		index <- matrix(NA, 4, 4)
		np <- 11
		index[col(index) != row(index)] <- 1:np
		nuc.mutation.rates <- matrix(0, nrow=4, ncol=4)
		nuc.mutation.rates<-matrix(rates[index], dim(index))
		rownames(nuc.mutation.rates) <- n2s(0:3)
		colnames(nuc.mutation.rates) <- n2s(0:3)
		nuc.mutation.rates[4,3] = 1
		diag(nuc.mutation.rates) <- 0
		diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
		return(nuc.mutation.rates)
	}
}


# Create a codon mutation model index, rows=from, cols=to, nuc order as in seqinr coding
#
# Note it does not assume a symmetric nucleotide matrix nor a symmetric codon mutation rate matrix
# Stop codons are included
# @param nuc.mutation.rates are the rates from (row) to (col) based on the nucleotide model
# @return codon.mutation.rates matrix of instantaneous mutation rates for codons
CreateCodonMutationMatrixIndex <- function() {
	nuc.rates.index = matrix(1:16, 4, 4)
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
				codon.mutation.rates[i,j] <- nuc.rates.index[1+codon.sets[i,mismatch.position], 1+codon.sets[j, mismatch.position]] #nucs numbered from 0:3, rows are 1:4, thus the add 1
			}
		}
		codon.name <- paste(n2s(as.numeric(codon.sets[i,])), collapse="")
		rownames(codon.mutation.rates)[i] <- codon.name
		colnames(codon.mutation.rates)[i] <- codon.name
	}
	codon.mutation.rates[codon.mutation.rates==0] = 17
	return(codon.mutation.rates)
}


# Create a codon mutation model, rows=from, cols=to, nuc order as in seqinr coding
#
# Note it does not assume a symmetric nucleotide matrix nor a symmetric codon mutation rate matrix
# Stop codons are included
# @param nuc.mutation.rates are the rates from (row) to (col) based on the nucleotide model
# @return codon.mutation.rates matrix of instantaneous mutation rates for codons
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
	diag(codon.mutation.rates) <- 0
	diag(codon.mutation.rates) <- -rowSums(codon.mutation.rates)
	return(codon.mutation.rates)
}


# Create a codon mutation model based on Goldman and Yang model, rows=from, cols=to, nuc order as in seqinr coding
#
# Note it does not assume a symmetric nucleotide matrix nor a symmetric codon mutation rate matrix
# Stop codons are included
# @param nuc.mutation.rates are the rates from (row) to (col) based on the nucleotide model
# @return codon.mutation.rates matrix of instantaneous mutation rates for codons
CreateCodonMutationMatrixGoldmanYang <- function(x, base.freqs) {
	kappa.par = x[1]
	omega.par = x[2]
	codon.sets <- expand.grid(0:3, 0:3, 0:3)
	codon.sets <- data.frame(first=codon.sets[,3], second=codon.sets[,2], third=codon.sets[,1]) #reordering to group similar codons
	n.codons <- dim(codon.sets)[1]
	codon.mutation.rates <- matrix(data=0, nrow=n.codons, ncol=n.codons)
	rownames(codon.mutation.rates) <- rep("",n.codons)
	colnames(codon.mutation.rates) <- rep("",n.codons)
	codon.set.translate <- apply(codon.sets, 2, n2s)
	codon.name <- apply(codon.set.translate, 1, paste, collapse="")
	aa.translation <- sapply(codon.name,TranslateCodon, numcode=numcode)
	for (i in sequence(n.codons)) {
		for (j in sequence(n.codons)) { #synonymous
			if(aa.translation[i] == aa.translation[j]){
				if(sum(codon.sets[i,] == codon.sets[j,])==2) { #means that two of the bases match
					mismatch.position <- which(codon.sets[i,] != codon.sets[j,])
					if(codon.sets[i,mismatch.position] == 0 & codon.sets[j,mismatch.position] == 2 | codon.sets[i,mismatch.position] == 2 & codon.sets[j,mismatch.position] == 0 | codon.sets[i,mismatch.position] == 1 & codon.sets[j,mismatch.position] == 3 | codon.sets[i,mismatch.position] == 3 & codon.sets[j,mismatch.position] == 1){#AtoGtransition or GtoAtransition or CtoTtransition or TtoCtransition
						codon.mutation.rates[i,j] = base.freqs[j] * kappa.par
					}else{ #transversion
						codon.mutation.rates[i,j] = base.freqs[j]
					}
				}				
			}else{ #nonsynonymous
				if(sum(codon.sets[i,] == codon.sets[j,])==2) { #means that two of the bases match
					mismatch.position <- which(codon.sets[i,] != codon.sets[j,])
					if(codon.sets[i,mismatch.position] == 0 & codon.sets[j,mismatch.position] == 2 | codon.sets[i,mismatch.position] == 2 & codon.sets[j,mismatch.position] == 0 | codon.sets[i,mismatch.position] == 1 & codon.sets[j,mismatch.position] == 3 | codon.sets[i,mismatch.position] == 3 & codon.sets[j,mismatch.position] == 1){#AtoGtransition or GtoAtransition or CtoTtransition or TtoCtransition
						codon.mutation.rates[i,j] = base.freqs[j] * kappa.par * omega.par
					}else{ #transversion
						codon.mutation.rates[i,j] = base.freqs[j] * omega.par
					}
				}
			}
		}
	}
	rownames(codon.mutation.rates) <- colnames(codon.mutation.rates) <- codon.name
	diag(codon.mutation.rates) <- 0
	diag(codon.mutation.rates) <- -rowSums(codon.mutation.rates)
	return(codon.mutation.rates)
}


CodonNumericToString <- function(x) { #remember that codon numbers start at 1
	return(n2s(x, levels=words(length=3), base4=FALSE))
}


CodonStringToNumeric <- function(x) { #remember that codon numbers start at 1
	triplet <- which(words(length=3)==x)
	if(length(triplet) == 0){
		triplet <- NA
	}
	return(triplet)
}

NucleotideStringToNumeric <- function(x) { #remember that codon numbers start at 1
	singlet <- which(words(length=1)==x)
	if(length(singlet) == 0){
		singlet <- NA
	}
	return(singlet)
}


ConvertCodonNumericDataToAAData <- function(codon.data, numcode=1) {
	aa.data <- codon.data
	for (row.index in sequence(dim(codon.data)[1])) {
		for (col.index in 2:(dim(codon.data)[2])) {
			#if(is.na(codon.data[row.index, col.index])){
			if(codon.data[row.index, col.index] == 65){
				aa.data[row.index, col.index] = "NA"
			}else{
				aa.data[row.index, col.index] <- TranslateCodon(CodonNumericToString(codon.data[row.index, col.index]), numcode=numcode)
			}
		}
	}
	return(aa.data)
}


#Problem with this is that it requires an assumption about frequencies. Flat
#  amino acid frequencies != flat codon frequences != frequencies of codons at
#  equilibrium given nucleotide model for many models. So, lets just go
#  directly to the codon model
#Create an amino acid mutation model, rows=from, cols=to
#
# Note it does not assume a symmetric nucleotide matrix nor a symmetric codon mutation rate matrix
# @param codon.mutation.rates are the rates from (row) to (col) based on the nucleotide model
# @param codon.freqs are the frequencies of each codon. By default, flat prior
# @param numcode is the genetic code as in seqinr, set by default to standard
# @return aa.mutation.rates matrix of instantaneous mutation rates for aa
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


GetPairwiseProteinFixationProbabilityArbitraryLength <- function(protein1, protein2, protein_op, s, aa.distances, nsites, C=2, Phi=0.5, q=4e-7, Ne=5e6){
	d1 <- GetProteinProteinDistance(protein1,protein_op,aa.distances)
	d2 <- GetProteinProteinDistance(protein2,protein_op,aa.distances)
	if(length(d1)!=length(d2)) #throw error if length of proteins are not the same
    stop("error: 2 proteins are of different lengths!")
	if(length(d1)==1){ #only one amino acid
		return(GetPairwiseProteinFixationProbabilitySingleSite(d1, d2, s, nsites=nsites, C=C,Phi=Phi,q=q,Ne=Ne))
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
			return(GetPairwiseProteinFixationProbabilitySingleSite(d1[pos], d2[pos], s[pos], nsites=nsites, C=C,Phi=Phi,q=q,Ne=Ne))
		}
	}
}


GetPairwiseProteinFixationProbabilitySingleSite <- function(d1, d2, nsites, C=2, Phi=0.5, q=4e-7, Ne=5e6){
#	if((d1==d2)||(s==0)) #When the fitnesses are the same, neutral case, pure drift
#   return(1/(2*Ne))
#	else{
	fit_ratio <- exp(-(C+(C/nsites))*Phi*q*(d1-d2)) #f1/f2
	if(fit_ratio==Inf) #1 is much better than 2 (the mutant)
	return(0)
	else if(fit_ratio==1)
	return(1/(2*Ne))
	else
	return((1-fit_ratio)/(1-fit_ratio^(2*Ne)))
#	}
}


TranslateCodon <- function(codon.string, numcode) {
	return(translate(s2c(codon.string), numcode=numcode))
}


GetProteinProteinDistance <- function(protein1, protein2, aa.distances){
	if(length(protein1)!=length(protein2)) #throw error if length of proteins are not the same
    stop("error: 2 proteins are of different lengths!")
	site_d <- function(k){
		if(protein1[k]=="*" || protein2[k]=="*") {
			warning("You have a stop codon in your sequence. This was treated as having a very large difference from other amino acids, but you probably want to exclude such sites. It may also be that your numcode is not appropriate for your data, and perhaps you want one that works for invertebrate mitochondria, chloroplasts, etc.")
			return(100*max(aa.distances))
		}
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
		return(aa.distances[protein1[k], protein2[k]])
	}
	d <- sapply(c(1:length(protein1)),site_d,simplify=TRUE) 
	return(d)
}


FastCreateAllCodonFixationProbabilityMatrices <- function(aa.distances=CreateAADistanceMatrix(), nsites, C=2, Phi=0.5, q=4e-7, Ne=5e6, include.stop.codon=TRUE, numcode=1, flee.stop.codon.rate=0.9999999) {
	codon.sets <- expand.grid(0:3, 0:3, 0:3)
	codon.sets <- data.frame(first=codon.sets[,3], second=codon.sets[,2], third=codon.sets[,1]) #reordering to group similar codons
	n.codons <- dim(codon.sets)[1]
	codon.names <- rep("", n.codons)
	for (i in sequence(n.codons)) {
		codon.names[i] <- paste(n2s(as.numeric(codon.sets[i,])), collapse="")
	}
	codon.aa <- sapply(codon.names, TranslateCodon, numcode=numcode)
	unique.aa <- unique(codon.aa)
	codon.fixation.probs <- array(data=0, dim=c(n.codons, n.codons, length(unique.aa)), dimnames=list(codon.names, codon.names, unique.aa))
	for (i in sequence(n.codons)) {
		for (j in sequence(n.codons)) {
			if(sum(codon.sets[i,] == codon.sets[j,])>=2) { #match at two or three sites of three
				for (k in sequence(length(unique.aa))) {
					aa1 <- codon.aa[i]
					aa2 <- codon.aa[j]
					if(aa1!="*" && aa2!="*" && unique.aa[k]!="*") { #says we cannot mutate to stop codons and stop codons can never be optimal
						d1 <- GetProteinProteinDistance(protein1=aa1, protein2=unique.aa[k], aa.distances=aa.distances)
						d2 <- GetProteinProteinDistance(protein1=aa2, protein2=unique.aa[k], aa.distances=aa.distances)
						codon.fixation.probs[i,j, k] <- GetPairwiseProteinFixationProbabilitySingleSite(d1, d2, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne)
					}else {
#						We have dropped s from the model as it is now explained through grantham like distances:
#						if(s==0) { #handles stop codon case where neutral, so could possibly go into and out of stop codons
#							codon.fixation.probs[i,j, k] <- 0
#						}else {
							if(aa2!="*" && unique.aa[k]!="*") {
								codon.fixation.probs[i,j, k] <- 0 #Old = if we are somehow in a stop codon, have a very high rate of moving away from this; New = make is zero because in theory our model should use selection to kill these but infinite selection is rather harsh.
#							}
						}
					}
				}
			}
		}
	}
	codon.fixation.probs[,,"*"] = 0 
	return(codon.fixation.probs)
}


DiagArray <- function (dim){
    n <- dim[2]
    d <- seq(1, n*n, by=n+1)
    as.vector(outer(d, seq(0, by=n*n, length=prod(dim[-1:-2])), "+"))
}


FastCreateAllCodonFixationProbabilityMatricesSetToOne <- function(numcode=1) {
	codon.sets <- expand.grid(0:3, 0:3, 0:3)
	codon.sets <- data.frame(first=codon.sets[,3], second=codon.sets[,2], third=codon.sets[,1]) #reordering to group similar codons
	n.codons <- dim(codon.sets)[1]
	codon.names <- rep("", n.codons)
	for (i in sequence(n.codons)) {
		codon.names[i] <- paste(n2s(as.numeric(codon.sets[i,])), collapse="")
	}
	codon.aa <- sapply(codon.names, TranslateCodon, numcode=numcode)
	unique.aa <- unique(codon.aa)
	codon.fixation.rates <- array(data=1, dim=c(n.codons, n.codons, length(unique.aa)), dimnames=list(codon.names, codon.names, unique.aa))
	return(codon.fixation.rates)
}


CreateCodonFixationProbabilityMatrix <- function(aa_op, s, aa.distances, nsites, C=2, Phi=0.5,q=4e-7,Ne=5e6, include.stop.codon=TRUE, numcode=1){
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
					codon.fixation.rates[i,j] <- GetPairwiseProteinFixationProbabilitySingleSite(d1, d2, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne)
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
	diag(codon.fixation.rates) <- 0 #b/c we dont want these included when calculating diag
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


# Get likelihood for a given AA site given tree and Q matrix, assuming the optimal AA is known
# @param aa.data is data in corHMM format (one column species, one column state)
# @param phy is a phylo object
# @param Q_aa is the transition matrix Q for a given optimal aa at this site
# @param charnum is which character to use in the aa.data matrix
# @param root.p is the root frequency: NULL if equilibrium, maddfitz, or a vector
# @param return.all, if TRUE, allows return of everything from rayDISC rather than just the log likelihood
GetLikelihoodSAC_AAForSingleCharGivenOptimum <- function(aa.data, phy, Q_aa, charnum=1, root.p=NULL, return.all=FALSE) {
	#	result <- rayDISC(phy=phy, data=aa.data, ntraits=1, charnum=charnum, p=Q_aa, root.p=root.p, node.states="marginal")
	#Makes dev.rayDISC available:
	#dev.raydisc <- corHMM:::dev.raydisc
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode
	nl <- nrow(Q_aa)
	#Now we need to build the matrix of likelihoods to pass to dev.raydisc:
	liks <- matrix(0, nb.tip + nb.node, nl)
	#Now loop through the tips.
	for(i in 1:nb.tip){
		#The codon at a site for a species is not NA, then just put a 1 in the appropriate column.
		#Note: We add charnum+1, because the first column in the data is the species labels:
		if(!is.na(aa.data[i,charnum+1])){
			liks[i,aa.data[i,charnum+1]] <- 1
		}else{
			#If here, then the site has no data, so we treat it as ambiguous for all possible codons. Likely things might be more complicated, but this can modified later:
			liks[i,] <- 1
		}
	}	
	diag(Q_codon) = 0
	diag(Q_codon) = -rowSums(Q_codon)
	#The result here is just the likelihood: 
	result <- -dev.raydisc(p=NULL, phy=phy, liks=liks, Q=Q_aa, rate=NULL, root.p=root.p)
	ifelse(return.all, stop("return all not currently implemented"), return(result))
}


# Get likelihood for a given codon site given tree and Q matrix, assuming the optimal AA is known
# @param codon.data is data in corHMM format (one column species, one column state, states from CodonStringToNumeric)
# @param phy is a phylo object
# @param Q_codon is the transition matrix Q for a given optimal aa at this site
# @param charnum is which character to use in the aa.data matrix
# @param root.p is the root frequency: NULL if equilibrium, maddfitz, or a vector
# @param return.all, if TRUE, allows return of everything from rayDISC rather than just the log likelihood
GetLikelihoodSAC_CodonForSingleCharGivenOptimum <- function(charnum=1, codon.data, phy, Q_codon, root.p=NULL, scale.factor, return.all=FALSE) {
	#result <- rayDISC(phy=phy, data=codon.data, ntraits=1, charnum=charnum, p=Q_codon, root.p=root.p)
	#Makes dev.rayDISC available:
	#dev.raydisc <- corHMM:::dev.raydisc
	#Check for any missing data -- all taxa missing data for a site are simply removed and the likelihood is inferred from the reduced dataset. Consider, for example, that a taxon is missing the first codon position base. Many different AA could fit that and so it drives likelihood optimum towards low s.
	
	### Old way of dealing with missing data -- we would just remove the taxa. Now we set it as ambiguous and ignore it during our tree traversal of the final likelihood.
	#	if(any(is.na(codon.data[,charnum+1]))){
	#		print(charnum)
	#		bad.taxa = which(is.na(codon.data[,charnum+1]))
	#		phy = drop.tip(phy, phy$tip.label[bad.taxa])
	#		codon.data = codon.data[-bad.taxa,]
	#	}
	###
	
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode
	nl <- nrow(Q_codon[[1]])
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
		}
	}
	
	#The result here is just the likelihood: 
	result <- -FinishLikelihoodCalculation(phy=phy, liks=liks, Q=Q_codon, root.p=root.p)
	ifelse(return.all, stop("return all not currently implemented"), return(result))
}


GetLikelihoodSAC_CodonForManyCharGivenFixedOptimumAndQAndRoot <- function(codon.data, phy, Q_codon, root.p=NULL, return.all=FALSE) {
	return(sum(sapply(seq(from=1, to=dim(codon.data)[2]-1, by=1), GetLikelihoodSAC_CodonForSingleCharGivenOptimum, codon.data=codon.data, phy=phy, Q_codon=Q_codon, root.p=root.p, return.all=return.all)))
}


GetLikelihoodSAC_CodonForManyCharVaryingBySite <- function(codon.data, phy, Q_codon_array, root.p_array=NULL, aa.optim_array, codon_mutation_matrix, Ne, rates, numcode) {
	
	nsites <- dim(codon.data$unique.site.patterns)[2]-1
	final.likelihood.vector <- rep(NA, nsites)
	equilibrium.codon.freq <- root.p_array / sum(root.p_array)	
	unique.aa <- GetMatrixAANames(numcode)
	#print(Q_codon_array)
	#Q_array codon mutation matrix multiplication here -- need to do this here because we need our scaling factor:
	for(k in 1:21){ 
		Q_codon_array[,,unique.aa[k]] = 2 * Ne * (codon_mutation_matrix * Q_codon_array[,,unique.aa[k]])  
		diag(Q_codon_array[,,unique.aa[k]]) = 0
		diag(Q_codon_array[,,unique.aa[k]]) = -rowSums(Q_codon_array[,,unique.aa[k]])
	}
	#Put the na.rm=TRUE bit here just in case -- when the amino acid is a stop codon, there is a bunch of NaNs. Should be fixed now.
	scale.factor <- -sum(Q_codon_array[DiagArray(dim(Q_codon_array))] * equilibrium.codon.freq, na.rm=TRUE)
	
	## This is obviously not very elegant, but not sure how else to code it to store this stuff in this way -- WORK IN PROGRESS:
	expQt <- NULL
	expQt$K <- GetExpQt(phy=phy, Q=Q_codon_array[,,"K"], scale.factor=scale.factor, rates=rates) 
	expQt$N <- GetExpQt(phy=phy, Q=Q_codon_array[,,"N"], scale.factor=scale.factor, rates=rates)
	expQt$T <- GetExpQt(phy=phy, Q=Q_codon_array[,,"T"], scale.factor=scale.factor, rates=rates)
	expQt$R <- GetExpQt(phy=phy, Q=Q_codon_array[,,"R"], scale.factor=scale.factor, rates=rates)
	expQt$S <- GetExpQt(phy=phy, Q=Q_codon_array[,,"S"], scale.factor=scale.factor, rates=rates)
	expQt$I <- GetExpQt(phy=phy, Q=Q_codon_array[,,"I"], scale.factor=scale.factor, rates=rates)
	expQt$M <- GetExpQt(phy=phy, Q=Q_codon_array[,,"M"], scale.factor=scale.factor, rates=rates)
	expQt$Q <- GetExpQt(phy=phy, Q=Q_codon_array[,,"Q"], scale.factor=scale.factor, rates=rates)
	expQt$H <- GetExpQt(phy=phy, Q=Q_codon_array[,,"H"], scale.factor=scale.factor, rates=rates)
	expQt$P <- GetExpQt(phy=phy, Q=Q_codon_array[,,"P"], scale.factor=scale.factor, rates=rates)
	expQt$L <- GetExpQt(phy=phy, Q=Q_codon_array[,,"L"], scale.factor=scale.factor, rates=rates)
	expQt$E <- GetExpQt(phy=phy, Q=Q_codon_array[,,"E"], scale.factor=scale.factor, rates=rates)
	expQt$D <- GetExpQt(phy=phy, Q=Q_codon_array[,,"D"], scale.factor=scale.factor, rates=rates)
	expQt$A <- GetExpQt(phy=phy, Q=Q_codon_array[,,"A"], scale.factor=scale.factor, rates=rates)
	expQt$G <- GetExpQt(phy=phy, Q=Q_codon_array[,,"G"], scale.factor=scale.factor, rates=rates)
	expQt$V <- GetExpQt(phy=phy, Q=Q_codon_array[,,"V"], scale.factor=scale.factor, rates=rates)
	expQt$Y <- GetExpQt(phy=phy, Q=Q_codon_array[,,"Y"], scale.factor=scale.factor, rates=rates)
	expQt$C <- GetExpQt(phy=phy, Q=Q_codon_array[,,"C"], scale.factor=scale.factor, rates=rates)
	expQt$W <- GetExpQt(phy=phy, Q=Q_codon_array[,,"W"], scale.factor=scale.factor, rates=rates)
	expQt$F <- GetExpQt(phy=phy, Q=Q_codon_array[,,"F"], scale.factor=scale.factor, rates=rates)

	#Generate matrix of root frequencies for each optimal AA:
	root.p_array <- matrix(root.p_array, nrow=dim(Q_codon_array)[2], ncol=21)
	root.p_array <- t(root.p_array)
	root.p_array <- root.p_array / rowSums(root.p_array)
	rownames(root.p_array) <- unique.aa
	
	for (i in sequence(nsites)) {
		final.likelihood.vector[i] <- GetLikelihoodSAC_CodonForSingleCharGivenOptimum(charnum=i, codon.data=codon.data$unique.site.patterns, phy=phy, Q_codon=expQt[[aa.optim_array[i]]], root.p=root.p_array[aa.optim_array[i],], scale.factor=scale.factor, return.all=FALSE)
	}
	return(final.likelihood.vector)
}


GetLikelihoodGoldYang_CodonForManyCharVaryingBySite <- function(codon.data, phy, root.p_array=NULL, codon_mutation_matrix, numcode) {
	nsites <- dim(codon.data$unique.site.patterns)[2] - 1
	final.likelihood.vector <- rep(NA, nsites)
	if(is.null(root.p_array)) {
		#Generate matrix of frequencies for each site -- NULL represents the empirical distribution at the tips:
		root.p_array <- matrix(0, nrow=nsites, ncol=dim(codon_mutation_matrix)[2])
		for(i in 1:nsites){
			codon.freqs <- table(codon.data$unique.site.patterns[,i+1])
			tmp.vector <- numeric(dim(codon_mutation_matrix)[2])
			tmp.vector[as.numeric(names(codon.freqs))] = codon.freqs
			root.p_array[i,] <- tmp.vector/sum(tmp.vector)
		}
	}
	for (i in sequence(nsites)) {
		final.likelihood.vector[i] <- GetLikelihoodSAC_CodonForSingleCharGivenOptimum(charnum=i, codon.data=codon.data$unique.site.patterns, phy=phy, Q_codon=codon_mutation_matrix, return.all=FALSE)
	}
	return(sum(final.likelihood.vector * codon.data$site.pattern.counts))
}


GetLikelihoodNucleotideForManyCharVaryingBySite <- function(nuc.data, phy, nuc.mutation.rates, include.gamma=FALSE, rates.k=NULL, ncats=NULL, root.p_array=NULL) {
	nsites <- dim(nuc.data$unique.site.patterns)[2]-1
	final.likelihood.vector <- rep(NA, nsites)
	if(is.null(root.p_array)) {
	#Generate matrix of equal frequencies for each site:
		root.p_array <- rep(0.25, 4)
	}
	#Rescaling Q matrix in order to have a 1 nucleotide change per site if the branch length was 1:
	diag(nuc.mutation.rates) = 0
	nuc.mutation.rates = t(nuc.mutation.rates * root.p_array)	
	diag(nuc.mutation.rates) = -rowSums(nuc.mutation.rates)
	scale.factor <- -sum(diag(nuc.mutation.rates) * root.p_array)
	
	expQt <- GetExpQt(phy=phy, Q=nuc.mutation.rates, scale.factor=scale.factor, rates=rates.k) 
	for(i in sequence(nsites)) {
		final.likelihood.vector[i] <- GetLikelihoodSAC_CodonForSingleCharGivenOptimum(charnum=i, codon.data=nuc.data$unique.site.patterns, phy=phy, Q_codon=expQt, root.p=root.p_array, scale.factor=scale.factor, return.all=FALSE)
	}
	return(final.likelihood.vector)
}


GetLikelihoodSAC_CodonForManyCharGivenAllParams <- function(x, codon.data, phy, aa.optim_array=NULL, root.p_array=NULL, numcode=1, aa.properties=NULL, nuc.model, codon.index.matrix, include.gamma, ncats, k.levels=0, logspace=FALSE, verbose=TRUE, neglnl=FALSE) {
	if(logspace) {
		x = exp(x)
	}
	if(include.gamma == TRUE){
		shape = x[length(x)]
		x = x[-length(x)]
	}
	C=2
	Phi <- 0.5
	q <- 4e-7
	alpha <- x[1]
	beta <- x[2]
	gamma = x[3]
	Ne <- x[4]
	
	if(k.levels > 0){
		if(nuc.model == "JC") {
			base.freqs=c(x[5:7], 1-sum(x[5:7]))
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
		}
		if(nuc.model == "GTR") {
			base.freqs=c(x[5:7], 1-sum(x[5:7]))
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[10:length(x)], model=nuc.model, base.freqs=base.freqs)
		}
		if(nuc.model == "UNREST") {
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[10:length(x)], model=nuc.model)
		}		
	}else{
		if(nuc.model == "JC") {
			base.freqs=c(x[5:7], 1-sum(x[5:7]))
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
		}
		if(nuc.model == "GTR") {
			base.freqs=c(x[5:7], 1-sum(x[5:7]))
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[8:length(x)], model=nuc.model, base.freqs=base.freqs)
		}
		if(nuc.model == "UNREST") {
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[8:length(x)], model=nuc.model)
		}		
	}
	#codon_mutation_matrix = CreateCodonMutationMatrix(nuc.mutation.rates) #We now make an index matrix first then just place the nucleotide rates into it:
	codon_mutation_matrix = c(as.vector(nuc.mutation.rates), 0)[codon.index.matrix]
	nsites <- dim(codon.data$unique.site.patterns)[2]-1
	
	if(include.gamma==TRUE){
		rates.k <- DiscreteGamma(shape, ncats)
		final.likelihood.mat = matrix(0, nrow=ncats, ncol=nsites)		
		for(k in sequence(ncats)){
			if(k.levels > 0){
				aa.distances <- CreateAADistanceMatrix(alpha=alpha*rates.k[k], beta=beta*rates.k[k], gamma=gamma*rates.k[k], aa.properties=aa.properties, normalize=FALSE, poly.params=x[8:9], k=k.levels)
			}else{
				aa.distances <- CreateAADistanceMatrix(alpha=alpha*rates.k[k], beta=beta*rates.k[k], gamma=gamma*rates.k[k], aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
			}
			Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, flee.stop.codon.rate=0.9999999)
			final.likelihood.mat[k,] = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data, phy, Q_codon_array, root.p_array=root.p_array, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode)
		}
		likelihood <- sum(log(colMeans(exp(final.likelihood.mat))) * codon.data$site.pattern.counts)
	}else{
		if(k.levels > 0){
			aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=x[8:9], k=k.levels)
		}else{
			aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
		}		
		Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, flee.stop.codon.rate=0.9999999)
		final.likelihood = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data, phy, Q_codon_array, root.p_array=root.p_array, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode)
		likelihood <- sum(final.likelihood * codon.data$site.pattern.counts)
	}
	
	if(neglnl) {
		likelihood <- -1 * likelihood
	}
	if(verbose) {
		results.vector <- c(likelihood, C, Phi, q, alpha, beta, gamma, Ne)
		names(results.vector) <- c("likelihood", "C", "Phi", "q", "alpha", "beta", "gamma", "Ne")
		print(results.vector)
	}
	if(is.na(likelihood) || is.nan(likelihood)){
		return(10000000000)
	}else{
		return(likelihood)
	}
}


GetLikelihoodGoldYang_CodonForManyCharGivenAllParams <- function(x, codon.data, phy, root.p_array=NULL, numcode=1, logspace=FALSE, verbose=TRUE, neglnl=FALSE) {
	if(logspace) {
		x = exp(x)
	}
	if(root.p_array==NULL){
		codon.freq <- rep(1/64, 64)
	}else{
		codon.freq <- as.matrix(codon.data[,-1])
		codon.freq <- table(codon.freq)/sum(table(codon.freq))
	}
	codon_mutation_matrix = CreateCodonMutationMatrixGoldYang(x, base.freqs=empirical.base.freq)
	likelihood <- GetLikelihoodGoldYang_CodonForManyCharVaryingBySite(codon.data, phy, root.p_array=empirical.base.freq, codon_mutation_matrix=codon_mutation_matrix)
	if(neglnl) {
		likelihood <- -1 * likelihood
	}
	if(verbose) {
		results.vector <- c(likelihood, x)
		names(results.vector) <- c("likelihood")
		print(results.vector)
	}
	return(likelihood)
}


GetLikelihoodNucleotideForManyCharGivenAllParams <- function(x, nuc.data, phy, root.p_array=NULL, numcode=1, nuc.model, include.gamma=FALSE, rates.k=NULL, ncats=NULL, logspace=FALSE, verbose=TRUE, neglnl=FALSE) {
	if(logspace) {
		x = exp(x)
	}
	if(include.gamma == TRUE){
		shape = x[1]
		x = x[-1]
	}

	if(length(x)==0){
		transition.rates <- 1	
	}else{
		transition.rates <- x[1:length(x)]
	}
	nsites <- dim(nuc.data$unique.site.patterns)[2]-1
	nuc.mutation.rates <- CreateNucleotideMutationMatrix(transition.rates, model=nuc.model)

	if(include.gamma==TRUE){
		rates.k <- DiscreteGamma(shape, ncats)
		final.likelihood.mat = matrix(0, nrow=ncats, ncol=nsites)		
		for(k in sequence(ncats)){
			final.likelihood.mat[k,] = GetLikelihoodNucleotideForManyCharVaryingBySite(nuc.data=nuc.data, phy=phy, nuc.mutation.rates=nuc.mutation.rates, rates.k=rates.k[k], root.p_array=root.p_array)
		}
		likelihood <- sum(log(colMeans(exp(final.likelihood.mat))) * nuc.data$site.pattern.counts)
	}else{
		final.likelihood = GetLikelihoodNucleotideForManyCharVaryingBySite(nuc.data=nuc.data, phy=phy, nuc.mutation.rates=nuc.mutation.rates, rates.k=NULL, root.p_array=root.p_array)
		likelihood <- sum(final.likelihood * nuc.data$site.pattern.counts)
	}
	
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


# Get likelihood for a given codon site given tree and Q matrix and determine the best likelihood across all possible amino acids
# @param codon.data is data in corHMM format (one column species, one column state, states from CodonStringToNumeric)
# @param phy is a phylo object
# @param Q_codon is the transition matrix Q for a given optimal aa at this site
# @param charnum is which character to use in the aa.data matrix
# @param root.p is the root frequency: NULL if equilibrium, maddfitz, or a vector
# @param return.all, if TRUE, allows return of everything from rayDISC rather than just the log likelihood
GetOptimalAAPerSite <- function(x, codon.data, phy, aa.optim_array=NULL, root.p_array=NULL, numcode=1, aa.properties=NULL, nuc.model, codon.index.matrix, include.gamma=FALSE, ncats=4, k.levels=0, logspace=FALSE, verbose=TRUE, neglnl=FALSE) {
	if(logspace) {
		x = exp(x)
	}
	if(include.gamma == TRUE){
		shape = x[length(x)]
		x = x[-length(x)]
	}
	C=2
	Phi <- 0.5
	q <- 4e-7
	alpha <- x[1]
	beta <- x[2]
	gamma = x[3]
	Ne <- x[4]
	
	if(k.levels > 0){
		if(nuc.model == "JC") {
			base.freqs=c(x[5:7], 1-sum(x[5:7]))
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
		}
		if(nuc.model == "GTR") {
			base.freqs=c(x[5:7], 1-sum(x[5:7]))
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[10:length(x)], model=nuc.model, base.freqs=base.freqs)
		}
		if(nuc.model == "UNREST") {
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[10:length(x)], model=nuc.model)
		}		
	}else{
		if(nuc.model == "JC") {
			base.freqs=c(x[5:7], 1-sum(x[5:7]))
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
		}
		if(nuc.model == "GTR") {
			base.freqs=c(x[5:7], 1-sum(x[5:7]))
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[8:length(x)], model=nuc.model, base.freqs=base.freqs)
		}
		if(nuc.model == "UNREST") {
			nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[8:length(x)], model=nuc.model)
		}		
	}
	
	if(!is.null(codon.data$unique.site.patterns)){
		codon.data.list <- codon.data
		nsites <- dim(codon.data$unique.site.patterns)[2]-1
	}else{
		nsites <- dim(codon.data)[2]-1
		codon.data.list <- NULL
		codon.data.list$unique.site.patterns <- codon.data
		codon.data.list$site.pattern.counts <- rep(1, nsites)
	}
	
	codon_mutation_matrix = c(as.vector(nuc.mutation.rates), 0)[codon.index.matrix]
	optimal.vector.by.site <- rep(NA, nsites)
	unique.aa <- GetMatrixAANames(numcode)	
	optimal.aa.likelihood.mat <- matrix(0, nrow=length(unique.aa), ncol=nsites)
	
	for(i in 1:length(unique.aa)){
		if(unique.aa[i]=="*"){
			optimal.aa.likelihood.mat[i,] <- rep(-10000000000, nsites)
		}else{
			aa.optim_array = rep(unique.aa[i], nsites)
			if(include.gamma==TRUE){
				rates.k <- DiscreteGamma(shape, ncats)
				final.likelihood.mat = matrix(0, nrow=ncats, ncol=nsites)		
				for(k in sequence(ncats)){
					if(k.levels > 0){
						aa.distances <- CreateAADistanceMatrix(alpha=alpha*rates.k[k], beta=beta*rates.k[k], gamma=gamma*rates.k[k], aa.properties=aa.properties, normalize=FALSE, poly.params=x[8:9], k=k.levels)
					}else{
						aa.distances <- CreateAADistanceMatrix(alpha=alpha*rates.k[k], beta=beta*rates.k[k], gamma=gamma*rates.k[k], aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
					}							
					Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, flee.stop.codon.rate=1000) 
					tmp = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data.list, phy, Q_codon_array, root.p_array=root.p_array, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode)
					tmp[is.na(tmp)] = -10000000000
					final.likelihood.mat[k,] = tmp
				}
				optimal.aa.likelihood.mat[i,] <- log(colMeans(exp(final.likelihood.mat)))
			}else{
				if(k.levels > 0){
					aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=x[8:9], k=k.levels)
				}else{
					aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
				}						
				Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, flee.stop.codon.rate=1000) 
				tmp = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data.list, phy, Q_codon_array, root.p_array=root.p_array, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode)
				tmp[is.na(tmp)] = -10000000000
				final.likelihood = tmp
				optimal.aa.likelihood.mat[i,] <- final.likelihood
			}
		}
	}
	for(j in 1:nsites){
		optimal.vector.by.site[j] <- unique.aa[which.is.max(optimal.aa.likelihood.mat[,j])]
	}
	return(optimal.vector.by.site)
}


OptimizeEdgeLengthsGlobal <- function(x, codon.site.data, codon.site.counts, n.partitions, index.matrix, phy, aa.optim_array=NULL, root.p_array=NULL, numcode=1, aa.properties=NULL, nuc.model, codon.index.matrix=NULL, include.gamma=FALSE, ncats, k.levels, logspace=FALSE, verbose=TRUE, n.cores=NULL, neglnl=FALSE) {
	if(logspace) {
		x <- exp(x)
	}
	par.mat <- index.matrix
	par.mat[] <- c(x, 0)[index.matrix]
	#Puts edge lengths on tree:
	if(is.null(aa.optim_array)){
		if(nuc.model == "JC"){
			max.par = 0
		}			
		if(nuc.model == "GTR"){
			max.par = 5
		}
		if(nuc.model == "UNREST"){
			max.par = 11
		}
		if(include.gamma == TRUE){
			max.par = max.par + 1
		}		
		if(is.null(n.cores)){
			likelihood.vector <- c()
			for(partition.index in sequence(n.partitions)){
				phy$edge.length = par.mat[partition.index,(max.par+1):ncol(par.mat)]
				nuc.data = NULL
				nuc.data$unique.site.patterns = codon.site.data[[partition.index]]
				nuc.data$site.pattern.counts = codon.site.counts[[partition.index]]
				likelihood.vector = c(likelihood.vector, GetLikelihoodNucleotideForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, root.p_array=root.p_array[[partition.index]], numcode=numcode, nuc.model=nuc.model, include.gamma=include.gamma, ncats=ncats, logspace=logspace, verbose=verbose, neglnl=neglnl))
			}
			likelihood = sum(likelihood.vector)
		}else{
			MultiCoreLikelihood <- function(partition.index){
				phy$edge.length = par.mat[partition.index,(max.par+1):ncol(par.mat)]
				nuc.data = NULL
				nuc.data$unique.site.patterns = codon.site.data[[partition.index]]
				nuc.data$site.pattern.counts = codon.site.counts[[partition.index]]
				likelihood.tmp = GetLikelihoodNucleotideForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, root.p_array=root.p_array[[partition.index]], numcode=numcode, nuc.model=nuc.model, include.gamma=include.gamma, ncats=ncats, logspace=logspace, verbose=verbose, neglnl=neglnl)
				return(likelihood.tmp)
			}
			likelihood <- sum(unlist(mclapply(1:n.partitions, MultiCoreLikelihood, mc.cores=n.cores)))
		}
	}else{
		if(nuc.model == "JC"){
			max.par = 7
		}
		if(nuc.model == "GTR"){
			max.par = 7 + 5
		}
		if(nuc.model == "UNREST"){
			max.par = 7 + 11
		}
		if(include.gamma == TRUE){
			max.par = max.par + 1
		}		
		if(k.levels > 0){
			max.par = max.par + 2
		}
		if(is.null(n.cores)){
			likelihood.vector <- c()
			for(partition.index in sequence(n.partitions)){
				phy$edge.length = par.mat[partition.index,(max.par+1):ncol(par.mat)]
				codon.data = NULL
				codon.data$unique.site.patterns = codon.site.data[[partition.index]]
				codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
				likelihood.vector = c(likelihood.vector, GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], root.p_array=root.p_array[[partition.index]], numcode=numcode, aa.properties=aa.properties, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl))
			}
			likelihood = sum(likelihood.vector)
		}else{
			MultiCoreLikelihood <- function(partition.index){
				phy$edge.length = par.mat[partition.index,(max.par+1):ncol(par.mat)]
				codon.data = NULL
				codon.data$unique.site.patterns = codon.site.data[[partition.index]]
				codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
				likelihood.tmp = GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], root.p_array=root.p_array[[partition.index]], numcode=numcode, aa.properties=aa.properties, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl)
				return(likelihood.tmp)
			}
			likelihood <- sum(unlist(mclapply(1:n.partitions, MultiCoreLikelihood, mc.cores=n.cores)))
		}			
	}
	return(likelihood)
}


ComputeStartingBranchLengths <- function(phy, chars=NULL) {
	phy <- compute.brlen(phy, method="Grafen")
	if(!is.null(chars)) {
		dna.distances <- dist(chars[,-1], diag=TRUE, upper=TRUE)/length(chars[,-1])
		max.root.tip.distance <- max(dna.distances)/2 
		#Note: /2 b/c distance is from tip to root to tip; we only want from the tip to the root. This might not be the true max length (if tree isnt ultrametric) but it gets us in right ballpark to start 
		max.branching.times <- max(branching.times(phy))
		phy$edge.length <- phy$edge.length * max.root.tip.distance / max.branching.times
	} else {
		rescaling=phy$edge.length[1]
		phy$edge.length <- phy$edge.length / rescaling 
		#so that the first brlen is fixed at 1
	}
	return(phy)
}


DiscreteGamma <- function (shape, ncats){
    quantiles <- qgamma((1:(ncats - 1))/ncats, shape = shape, rate = shape)
    return(diff(c(0, pgamma(quantiles * shape, shape + 1), 1)) * ncats)
}


PlotBubbleMatrix <- function(x, main="", special=Inf, cex=1){
	diag(x) <- 0
	x<-x/max(x)
	plot(x=range(.5,.5+dim(x)[2]),y=-range(.5, .5+dim(x)[1]), xlab="", ylab="", type="n", sub=main,xaxt='n',yaxt='n', asp=1,bty="n")
	axis(side=2, at=-sequence(dim(x)[1]), labels=rownames(x), las=2, cex.axis=cex)
	axis(side=3, at=sequence(dim(x)[2]), labels=colnames(x), las=2, cex.axis=cex)
	
	#abline(h=-1:(-dim(x)[2]), v=1:(dim(x)[1]), col="gray", lty=3)
	abline(h=-range(special)[1],v=range(special)[1], lty=2)
	abline(h=-range(special)[2],v=range(special)[2], lty=2)
	for (i in sequence(dim(x)[2])) {
		for (j in sequence(dim(x)[1])) {
			bg="gray"
			if(i %in% special || j %in% special) {
				if(i %in% special){
					bg="green4"
				}else{
					bg="magenta3"	
				}
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
	split.characters[is.na(split.characters)] = 65
	return(split.characters)
}


DNAbinToNucleotideNumeric <- function(x, frame=0, corHMM.format=TRUE) {
	bound.characters <- sapply(as.character(x), paste, collapse="")
	#following fn is derived from code for uco in seqinr
	SplitToCodons <- function(seq.string, frame) {
		seq.string<-strsplit(seq.string, split="")[[1]]
		if (any(seq.string %in% LETTERS)) {
			seq.string <- tolower(seq.string)
		}
		return(sapply(splitseq(seq = seq.string, frame = frame, word = 1), NucleotideStringToNumeric))
	}
	split.characters <- t(sapply(bound.characters, SplitToCodons, frame=frame))
	colnames(split.characters) <- sequence(dim(split.characters)[2])
	if(corHMM.format) {
		split.characters<-cbind(data.frame(Taxa=rownames(split.characters)), data.frame(split.characters))
	}
	split.characters[is.na(split.characters)] = 65
	return(split.characters)
}


SitePattern <- function(codon.data, corHMM.format=TRUE, includes.optimal.aa=FALSE){
	if(includes.optimal.aa == TRUE){
		char.strings <- sapply(codon.data[,-1], paste, collapse="_")	
		tabled.strings <- table(char.strings)
		site.pattern.totals <- c()
		reduced.codon.data <- c()
		reduced.optimal.aa <- c()
		for(i in 1:length(tabled.strings)){
			site.pattern.totals <- c(site.pattern.totals, tabled.strings[i])
			reduced.codon.data.tmp <- unlist(strsplit(names(tabled.strings)[i], "_"))
			reduced.codon.data <- cbind(reduced.codon.data, as.numeric(reduced.codon.data.tmp[1:(length(reduced.codon.data.tmp)-1)]))
			reduced.optimal.aa <- c(reduced.optimal.aa, reduced.codon.data.tmp[length(reduced.codon.data.tmp)])
		}
		if(corHMM.format) {
			names.for.rows <- rownames(codon.data)[-dim(codon.data)[1]]
			split.characters<-cbind(data.frame(Taxa=names.for.rows), data.frame(reduced.codon.data))
			rownames(split.characters) <- names.for.rows
		}	
		names(site.pattern.totals) <- NULL
		obj <- NULL
		obj$unique.site.patterns <- split.characters
		obj$site.pattern.counts <- site.pattern.totals
		obj$optimal.aa <- reduced.optimal.aa
	}else{
		char.strings <- sapply(codon.data[,-1], paste, collapse="_")	
		tabled.strings <- table(char.strings)
		site.pattern.totals <- c()
		reduced.codon.data <- c()
		for(i in 1:length(tabled.strings)){
			site.pattern.totals <- c(site.pattern.totals, tabled.strings[i])
			reduced.codon.data <- cbind(reduced.codon.data, as.numeric(unlist(strsplit(names(tabled.strings)[i], "_"))))
		}
		if(corHMM.format) {
			split.characters<-cbind(data.frame(Taxa=rownames(codon.data)), data.frame(reduced.codon.data))
			rownames(split.characters) <- rownames(codon.data)
		}	
		names(site.pattern.totals) <- NULL
		obj <- NULL
		obj$unique.site.patterns <- split.characters
		obj$site.pattern.counts <- site.pattern.totals
	}
	return(obj)
}


GetMatrixAANames <-function(numcode){
	codon.sets <- expand.grid(0:3, 0:3, 0:3)
	codon.sets <- data.frame(first=codon.sets[,3], second=codon.sets[,2], third=codon.sets[,1]) #reordering to group similar codons
	codon.set.translate <- apply(codon.sets, 2, n2s)
	codon.name <- apply(codon.set.translate, 1, paste, collapse="")
	codon.aa <- sapply(codon.name,TranslateCodon, numcode=numcode)
	names(codon.aa ) = NULL
	unique.aa <- unique(codon.aa)
	return(unique.aa)
}


CodonEquilibriumFrequencies <- function(codon.data, aa.opt.vector, numcode){
	codon.sets <- expand.grid(0:3, 0:3, 0:3)
	codon.sets <- data.frame(first=codon.sets[,3], second=codon.sets[,2], third=codon.sets[,1]) #reordering to group similar codons
	codon.set.translate <- apply(codon.sets, 2, n2s)
	codon.name <- apply(codon.set.translate, 1, paste, collapse="")
	aa.translation <- sapply(codon.name,TranslateCodon, numcode=numcode)
	names(aa.translation) = NULL
	unique.aa <- unique(aa.translation)
	eq.freqs <- c()
	for(aa.id.index in sequence(21)) {
		cols <- which(aa.opt.vector == unique.aa[aa.id.index])
		eq.freqs.tmp <- rep(0, 64)
		for(col.index in sequence(length(cols))) {
			for(row.index in sequence(dim(codon.data)[1])) {
				if(codon.data[row.index, cols[col.index]]<65){
					eq.freqs.tmp[codon.data[row.index, cols[col.index]]] <- eq.freqs.tmp[codon.data[row.index, cols[col.index]]] + 1
				}
			}
		}
		eq.freqs <- c(eq.freqs, eq.freqs.tmp)
	}
	return(eq.freqs)
}


GetMaxName <- function(x) {
	return(names(table(x))[(which.is.max(table(x)))]) #note that this breaks ties at random
}



######################################################################################################################################
######################################################################################################################################
### Likelihood calculator -- Two step process
######################################################################################################################################
######################################################################################################################################

#Step 1: We perform exponentiation as few times as possible; in the case of selac, for example, we do it for each of the possible optimal amino acids and store the matrices.
GetExpQt <- function(phy, Q, scale.factor, rates=NULL){
	Q.scaled = Q * (1/scale.factor)
	if(!is.null(rates)){
		Q.scaled = Q.scaled * rates	
	}
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	expQt <- as.list(numeric(nb.tip + nb.node))
	TIPS <- 1:nb.tip
	comp <- numeric(nb.tip + nb.node)
	phy <- reorder(phy, "pruningwise")
	#Obtain an object of all the unique ancestors
	anc <- unique(phy$edge[,1])
	for (i  in seq(from = 1, length.out = nb.node)) {
		#the ancestral node at row i is called focal
		focal <- anc[i]
		#Get descendant information of focal
		desRows<-which(phy$edge[,1]==focal)
		desNodes<-phy$edge[desRows,2]
		for (desIndex in sequence(length(desRows))){
			expQt[[desNodes[desIndex]]] <- expm(Q.scaled * phy$edge.length[desRows[desIndex]], method=c("Ward77"))
		}
	}
	return(expQt)
}


#Step 2: Finish likelihood by taking our already exponentiated Q down the tee and simply re-traverse the tree and multiply by the observed likelihood. 
FinishLikelihoodCalculation <- function(phy, liks, Q, root.p){	
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	TIPS <- 1:nb.tip
	comp <- numeric(nb.tip + nb.node)
	phy <- reorder(phy, "pruningwise")
	#Obtain an object of all the unique ancestors
	anc <- unique(phy$edge[,1])
	for (i  in seq(from = 1, length.out = nb.node)) {
		#the ancestral node at row i is called focal
		focal <- anc[i]
		#Get descendant information of focal
		desRows<-which(phy$edge[,1]==focal)
		desNodes<-phy$edge[desRows,2]
		v <- 1
		for (desIndex in sequence(length(desRows))){
			if(desNodes[desIndex] <= nb.tip){
				if(sum(liks[desNodes[desIndex],]) < 2){
					v <- v * Q[[desNodes[desIndex]]] %*% liks[desNodes[desIndex],]
				}
			}else{
				v <- v * Q[[desNodes[desIndex]]] %*% liks[desNodes[desIndex],]	
			}
		}
		comp[focal] <- sum(v)
		liks[focal, ] <- v/comp[focal]
	}
	#Specifies the root:
	root <- nb.tip + 1L
	#If any of the logs have NAs restart search:
	if(is.nan(sum(log(comp[-TIPS]))) || is.na(sum(log(comp[-TIPS])))){
		return(10000000000)
	}
	else{
		loglik<- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
		if(is.infinite(loglik)){return(10000000000)}
	}
	loglik
}



######################################################################################################################################
######################################################################################################################################
### Main function -- SINGLE OR MANY GENE PARTITIONS IMPLEMENTATION -- Keeps alpha, beta, gamma, and Ne constant across genes
######################################################################################################################################
######################################################################################################################################

EstimateParametersCodonGlobal <- function(codon.data.path, n.partitions=NULL, phy, edge.length="optimize", edge.linked=TRUE, optimal.aa="optimize", nuc.model="GTR", gold.yang=FALSE, include.gamma=FALSE, ncats=4, numcode=1, k.levels=0, aa.properties=NULL, verbose=FALSE, n.cores=NULL, max.tol=.Machine$double.eps^0.25, fasta.rows.to.keep=NULL) {
	
	cat("Initializing data and model parameters...", "\n")
	
	partitions <- system(paste("ls -1 ", codon.data.path, "*.fasta", sep=""), intern=TRUE)

	if(is.null(n.partitions)){
		n.partitions <- length(partitions)
	}else{
		n.partitions = n.partitions
	}
	site.pattern.data.list <- as.list(numeric(n.partitions))
	site.pattern.count.list <- as.list(numeric(n.partitions))
	nsites.vector <- c()
	if(optimal.aa == "none"){
		empirical.base.freq.list <- as.list(numeric(n.partitions))
		for (partition.index in sequence(n.partitions)) {
			gene.tmp <- read.dna(partitions[partition.index], format='fasta')
			if(!is.null(fasta.rows.to.keep)){
				gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
			}else{
				gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
			}
			nucleotide.data <- DNAbinToNucleotideNumeric(gene.tmp)
			nucleotide.data <- nucleotide.data[phy$tip.label,]
			nsites.vector = c(nsites.vector, dim(nucleotide.data)[2] - 1)
			empirical.base.freq <- as.matrix(nucleotide.data[,-1])
			empirical.base.freq <- table(empirical.base.freq, deparse.level = 0)/sum(table(empirical.base.freq, deparse.level = 0))
			empirical.base.freq.list[[partition.index]] <- as.vector(empirical.base.freq[1:4])
			nucleotide.data <- SitePattern(nucleotide.data, includes.optimal.aa=FALSE)
			site.pattern.data.list[[partition.index]] = nucleotide.data$unique.site.patterns
			site.pattern.count.list[[partition.index]] = nucleotide.data$site.pattern.counts
		}
	}else{
		empirical.codon.freq.list <- as.list(numeric(n.partitions))
		aa.optim.list <- as.list(numeric(n.partitions))
		aa.optim.full.list <- as.list(numeric(n.partitions))
		for (partition.index in sequence(n.partitions)) {
			gene.tmp <- read.dna(partitions[partition.index], format='fasta')
			if(!is.null(fasta.rows.to.keep)){
				gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
			}else{
				gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
			}
			codon.data <- DNAbinToCodonNumeric(gene.tmp)
			codon.data <- codon.data[phy$tip.label,]
			nsites.vector = c(nsites.vector, dim(codon.data)[2] - 1)
			aa.data <- ConvertCodonNumericDataToAAData(codon.data, numcode=numcode)
			aa.optim <- apply(aa.data[, -1], 2, GetMaxName) #starting values for all, final values for majrule
			aa.optim.full.list[[partition.index]] <- aa.optim
			empirical.codon.freq.list[[partition.index]] <- CodonEquilibriumFrequencies(codon.data[,-1], aa.optim, numcode=numcode)
			aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
			colnames(aa.optim.frame.to.add) <- colnames(codon.data)
			codon.data <- rbind(codon.data, aa.optim.frame.to.add)
			codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
			site.pattern.data.list[[partition.index]] = codon.data$unique.site.patterns
			site.pattern.count.list[[partition.index]] = codon.data$site.pattern.counts
			aa.optim.list[[partition.index]] = codon.data$optimal.aa		
		}
	}
	opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = "100000", "ftol_rel" = max.tol)
	results.final <- c()
	if(nuc.model == "JC"){
		nuc.ip = NULL
		max.par.model.count = 0
	}
	if(nuc.model == "GTR"){
		nuc.ip = rep(1, 5)
		max.par.model.count = 5
	}
	if(nuc.model == "UNREST"){
		nuc.ip = rep(1, 11)
		max.par.model.count = 11
	}	
	if(optimal.aa=="none") {
		codon.index.matrix = NA
		if(include.gamma == TRUE){
			ip = c(1,nuc.ip)
			max.par.model.count = max.par.model.count + 1
		}else{
			ip = nuc.ip
		}
		upper = rep(21, length(ip))
		lower = rep(-21, length(ip))
		if(edge.length == "optimize"){
			index.matrix = matrix(0, n.partitions, length(ip)+length(phy$edge.length))
			index.matrix[1,] = 1:ncol(index.matrix)
			phy <- ComputeStartingBranchLengths(phy, chars=site.pattern.data.list[[partition.index]])
			ip.vector = c(ip, phy$edge.length)
			upper.vector = c(upper, rep(log(5), length(phy$edge.length)))
			lower.vector = c(lower, rep(-21, length(phy$edge.length)))
			if(edge.linked == TRUE){
				for(partition.index in 2:n.partitions){
					ip.vector = c(ip.vector, ip)
					upper.vector = c(upper.vector, upper)
					lower.vector = c(lower.vector, lower)
					index.matrix.tmp = numeric(length(ip))
					index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
					index.matrix.tmp = c(index.matrix.tmp,index.matrix[1,(max.par.model.count+1):ncol(index.matrix)])
					index.matrix[partition.index,] <- index.matrix.tmp
				}				
			}else{
				for(partition.index in 2:n.partitions){
					phy <- ComputeStartingBranchLengths(phy, chars=site.pattern.data.list[[partition.index]])
					ip.vector = c(ip.vector, c(ip, phy$edge.length))
					upper.vector = c(upper.vector, c(upper, rep(log(5), length(phy$edge.length))))
					lower.vector = c(lower.vector, c(lower, rep(-21, length(phy$edge.length))))
					index.matrix.tmp = numeric(length(ip)+length(phy$edge.length))
					index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
					index.matrix[partition.index,] <- index.matrix.tmp
				}
			}
			cat("Finished. Optimizing model parameters...", "\n")
			results.final <- nloptr(x0=log(ip.vector), eval_f = OptimizeEdgeLengthsGlobal, ub=upper.vector, lb=lower.vector, opts=opts, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, n.partitions=n.partitions, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=empirical.base.freq.list, numcode=numcode, aa.properties=aa.properties, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
			cat("Finished. Summarizing results...", "\n")			
		}
		mle.pars.mat <- index.matrix
		mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
		phy$edge.length <- apply(mle.pars.mat[,7:ncol(mle.pars.mat)], 2, weighted.mean, w=nsites.vector)
		loglik <- -(results.final$objective) #to go from neglnl to lnl
		np = max(index.matrix)
		s <- 0
		obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, AICc = NULL, mle.pars=mle.pars.mat, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, aa.optim=NULL, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, aa.properties=aa.properties, empirical.base.freqs=empirical.base.freq.list, max.tol=max.tol) 
		class(obj) = "selac"				
	}
	if(optimal.aa=="majrule" | optimal.aa=="optimize") {
		codon.index.matrix = CreateCodonMutationMatrixIndex()		
		cpv.starting.parameters <- GetAADistanceStartingParameters(aa.properties=aa.properties)
		if(include.gamma == TRUE){
			if(nuc.model == "JC"){
				if(k.levels == 0){
					ip = c(cpv.starting.parameters, 5e6, 0.25, 0.25, 0.25, 1)
					upper = c(rep(21, 4), 0, 0, 0, 21, 21)
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 0 + 1
				}else{
					ip = c(cpv.starting.parameters, 5e6, 0.25, 0.25, 0.25, 1, 1, 1)
					upper = c(rep(21, 4), 0, 0, 0, 21, 21, 21)
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 0 + 1 + 1
				}
			}
			if(nuc.model == "GTR") {
				if(k.levels == 0){
					ip = c(cpv.starting.parameters, 5e6, 0.25, 0.25, 0.25, nuc.ip, 1)
					upper = c(rep(21, 4), 0, 0, 0, rep(21, length(nuc.ip)), 21)
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 5 + 1
				}else{
					ip = c(cpv.starting.parameters, 5e6, 0.25, 0.25, 0.25, 1, 1, nuc.ip, 1)
					upper = c(rep(21, 4), 0, 0, 0, 21, 21, rep(21, length(nuc.ip)), 21)
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 5 + 1	+ 1					
				}
			}
			if(nuc.model == "UNREST") {
				if(k.levels == 0){
					ip = c(cpv.starting.parameters, 5e6, 0.25, 0.25, 0.25, nuc.ip, 1)
					upper = c(rep(21, 4), 0, 0, 0, rep(21, length(nuc.ip)), 21)
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 11 + 1
				}else{
					ip = c(cpv.starting.parameters, 5e6, 0.25, 0.25, 0.25, 1, 1, nuc.ip, 1)
					upper = c(rep(21, 4), 0, 0, 0, 21, 21, rep(21, length(nuc.ip)), 21)
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 11 + 1 + 1	
				}
			}
			if(edge.length == "optimize"){
				index.matrix = matrix(0, n.partitions, max.par.model.count+length(phy$edge.length))
				index.matrix[1,] = 1:ncol(index.matrix)
				phy <- ComputeStartingBranchLengths(phy, chars=site.pattern.data.list[[partition.index]])
				ip.vector = c(ip, phy$edge.length)
				upper.vector = c(upper, rep(log(5), length(phy$edge.length)))
				lower.vector = c(lower, rep(-21, length(phy$edge.length)))
				if(edge.linked == TRUE){
					for(partition.index in 2:n.partitions){
						if(nuc.model == "JC"){
							ip.vector = c(ip.vector, ip[5], ip[6], ip[7])
							upper.vector = c(upper.vector, c(upper[5], upper[6], upper[7]))
							lower.vector = c(lower.vector, c(lower[5], lower[6], lower[7]))
						}else{
							ip.vector = c(ip.vector, ip[5], ip[6], ip[7], nuc.ip)
							upper.vector = c(upper.vector, c(upper[5], upper[6], upper[7], rep(21, length(nuc.ip))))
							lower.vector = c(lower.vector, c(lower[5], lower[6], lower[7], rep(-21, length(nuc.ip))))								
						}
						index.matrix.tmp = numeric(max.par.model.count + length(phy$edge.length))
						#This fixes alpha, beta, gamma, Ne, and the edge lengths across all partitions:
						if(k.levels == 0){
							index.matrix.tmp[c(1:4,max.par.model.count)] = c(1:4,max.par.model.count)
						}else{
							index.matrix.tmp[c(1:4,8:9,max.par.model.count)] = c(1:4,8:9,max.par.model.count)
						}
						index.matrix.tmp[(max.par.model.count+1):ncol(index.matrix)] = index.matrix[1,(max.par.model.count+1):ncol(index.matrix)]
						index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
						index.matrix[partition.index,] <- index.matrix.tmp							
					}					
				}else{
					for(partition.index in 2:n.partitions){
						phy <- ComputeStartingBranchLengths(phy, chars=site.pattern.data.list[[partition.index]])
						if(nuc.model == "JC"){
							ip.vector = c(ip.vector, ip[5], ip[6], ip[7])
							upper.vector = c(upper.vector, c(upper[5], upper[6], upper[7]))
							lower.vector = c(lower.vector, c(lower[5], lower[6], lower[7]))
						}else{
							ip.vector = c(ip.vector, ip[5], ip[6], ip[7], nuc.ip)
							upper.vector = c(upper.vector, c(upper[5], upper[6], upper[7], rep(21, length(nuc.ip)), rep(log(5), length(phy$edge.length))))
							lower.vector = c(lower.vector, c(lower[5], lower[6], lower[7], rep(-21, length(nuc.ip)), rep(-21, length(phy$edge.length))))								
						}
						index.matrix.tmp = numeric(max.par.model.count + length(phy$edge.length))
						#This fixes alpha, beta, gamma, and Ne across all partitions:
						if(k.levels == 0){
							index.matrix.tmp[c(1:4,max.par.model.count)] = c(1:4,max.par.model.count)
						}else{
							index.matrix.tmp[c(1:4,8:9,max.par.model.count)] = c(1:4,8:9,max.par.model.count)
						}
						index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
						index.matrix[partition.index,] <- index.matrix.tmp							
					}
				}
				if(optimal.aa == "optimize"){
					cat("Finished. Optimizing model parameters...", "\n")
					results.final <- nloptr(x0=log(ip.vector), eval_f = OptimizeEdgeLengthsGlobal, ub=upper.vector, lb=lower.vector, opts=opts, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, n.partitions=n.partitions, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=empirical.codon.freq.list, numcode=numcode, aa.properties=aa.properties, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
					cat("Finished. Finding optimal aa...", "\n")
					mle.pars.mat <- index.matrix
					mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
					aa.optim.list <- as.list(numeric(n.partitions))
					aa.optim.full.list <- as.list(numeric(n.partitions))
					for(partition.index in sequence(n.partitions)) {
						phy$edge.length = mle.pars.mat[partition.index, (max.par.model.count+1):ncol(mle.pars.mat)]
						gene.tmp <- read.dna(partitions[partition.index], format='fasta')
						if(!is.null(fasta.rows.to.keep)){
							gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
						}else{
							gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
						}
						codon.data <- DNAbinToCodonNumeric(gene.tmp)
						codon.data <- codon.data[phy$tip.label,]
						aa.optim.full.list[[partition.index]] = GetOptimalAAPerSite(x=log(mle.pars.mat[partition.index,1:max.par.model.count]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim.list[[partition.index]], root.p_array=empirical.codon.freq.list[[partition.index]], numcode=numcode, aa.properties=aa.properties, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, neglnl=TRUE)
						empirical.codon.freq.list[[partition.index]] <- CodonEquilibriumFrequencies(codon.data[,-1], aa.optim.full.list[[partition.index]], numcode=numcode)
						aa.optim.frame.to.add <- matrix(c("optimal", aa.optim.full.list[[partition.index]]), 1, dim(codon.data)[2])
						colnames(aa.optim.frame.to.add) <- colnames(codon.data)
						codon.data <- rbind(codon.data, aa.optim.frame.to.add)
						codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
						site.pattern.data.list[[partition.index]] = codon.data$unique.site.patterns
						site.pattern.count.list[[partition.index]] = codon.data$site.pattern.counts
						aa.optim.list[[partition.index]] = codon.data$optimal.aa		
					}
					cat("Finished. Final search given optimal aa...", "\n")
					results.final <- nloptr(x0=results.final$solution, eval_f = OptimizeEdgeLengthsGlobal, ub=upper.vector, lb=lower.vector, opts=opts, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, n.partitions=n.partitions, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=empirical.codon.freq.list, numcode=numcode, aa.properties=aa.properties, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
					cat("Finished. Summarizing results...", "\n")			
				}else{
					cat("Finished. Optimizing model parameters...", "\n")
					results.final <- nloptr(x0=log(ip.vector), eval_f = OptimizeEdgeLengthsGlobal, ub=upper.vector, lb=lower.vector, opts=opts, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, n.partitions=n.partitions, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=empirical.codon.freq.list, numcode=numcode, aa.properties=aa.properties, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
					cat("Finished. Summarizing results...", "\n")
				}
			}
			mle.pars.mat <- index.matrix
			mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
			phy$edge.length <- apply(mle.pars.mat[,(max.par.model.count+1):ncol(mle.pars.mat)], 2, weighted.mean, w=nsites.vector)
		}else{
			if(nuc.model == "JC"){
				if(k.levels == 0){
					ip = c(cpv.starting.parameters, 5e6, 0.25, 0.25, 0.25)
					upper = c(rep(21, 4), 0, 0, 0)
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 0 + 0
				}else{
					ip = c(cpv.starting.parameters, 5e6, 0.25, 0.25, 0.25, 0.001, 0.001)
					upper = c(rep(21, 4), 0, 0, 0, 21, 21)
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 0 + 0 + 2					
				}
			}
			if(nuc.model == "GTR") {
				if(k.levels == 0){
					ip = c(cpv.starting.parameters, 5e6, 0.25, 0.25, 0.25, nuc.ip)
					upper = c(rep(21, 4), 0, 0, 0, rep(21, length(nuc.ip)))
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 5 + 0
				}else{
					ip = c(cpv.starting.parameters, 5e6, 0.25, 0.25, 0.25, 0.001, 0.001, nuc.ip)
					upper = c(rep(21, 4), 0, 0, 0, 21, 21, rep(21, length(nuc.ip)))
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 5 + 0 + 2					
				}
			}
			if(nuc.model == "UNREST") {
				if(k.levels == 0){
					ip = c(cpv.starting.parameters, 5e6, nuc.ip)
					upper = c(rep(21, 4), 0, 0, 0, rep(21, length(nuc.ip)))
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 11 + 0
				}else{
					ip = c(cpv.starting.parameters, 5e6, 0.001, 0.001, nuc.ip)
					upper = c(rep(21, 4), 0, 0, 0, 21, 21, rep(21, length(nuc.ip)))
					lower = rep(-21, length(ip))
					max.par.model.count = 7 + 11 + 0 + 2					
				}
			}
			if(edge.length == "optimize"){
				index.matrix = matrix(0, n.partitions, max.par.model.count+length(phy$edge.length))
				index.matrix[1,] = 1:ncol(index.matrix)
				phy <- ComputeStartingBranchLengths(phy, chars=site.pattern.data.list[[partition.index]])
				ip.vector = c(ip, phy$edge.length)
				upper.vector = c(upper, rep(log(5), length(phy$edge.length)))
				lower.vector = c(lower, rep(-21, length(phy$edge.length)))
				if(edge.linked == TRUE){
					for(partition.index in 2:n.partitions){
						if(nuc.model == "JC"){
							ip.vector = c(ip.vector, ip[5], ip[6], ip[7])
							upper.vector = c(upper.vector, c(upper[5], upper[6], upper[7]))
							lower.vector = c(lower.vector, c(lower[5], lower[6], lower[7]))
						}else{
							ip.vector = c(ip.vector, ip[5], ip[6], ip[7], nuc.ip)
							upper.vector = c(upper.vector, c(upper[5], upper[6], upper[7], rep(21, length(nuc.ip))))
							lower.vector = c(lower.vector, c(lower[5], lower[6], lower[7], rep(-21, length(nuc.ip))))								
						}
						index.matrix.tmp = numeric(max.par.model.count + length(phy$edge.length))
						#This fixes s, alpha, beta, gamma, Ne, and the edge lengths across all partitions:
						if(k.levels == 0){
							index.matrix.tmp[c(1:4)] = c(1:4)
						}else{
							index.matrix.tmp[c(1:4,8:9)] = c(1:4,8:9)
						}
						index.matrix.tmp[(max.par.model.count+1):ncol(index.matrix)] = index.matrix[1,(max.par.model.count+1):ncol(index.matrix)]
						index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
						index.matrix[partition.index,] <- index.matrix.tmp														
					}					
				}else{
					for(partition.index in 2:n.partitions){
						phy <- ComputeStartingBranchLengths(phy, chars=site.pattern.data.list[[partition.index]])
						if(nuc.model == "JC"){
							ip.vector = c(ip.vector, ip[5], ip[6], ip[7])
							upper.vector = c(upper.vector, c(upper[5], upper[6], upper[7]))
							lower.vector = c(lower.vector, c(lower[5], lower[6], lower[7]))
						}else{
							ip.vector = c(ip.vector, ip[5], ip[6], ip[7], nuc.ip, phy$edge.length)
							upper.vector = c(upper.vector, c(upper[5], upper[6], upper[7], rep(21, length(nuc.ip)), rep(log(5), length(phy$edge.length))))
							lower.vector = c(lower.vector, c(lower[5], lower[6], lower[7], rep(-21, length(nuc.ip)), rep(-21, length(phy$edge.length))))								
						}
						index.matrix.tmp = numeric(max.par.model.count+length(phy$edge.length))
						#This fixes s, alpha, beta, gamma, and Ne across all partitions:
						if(k.levels == 0){
							index.matrix.tmp[c(1:4)] = c(1:4)
						}else{
							index.matrix.tmp[c(1:4,8:9)] = c(1:4,8:9)
						}
						index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
						index.matrix[partition.index,] <- index.matrix.tmp							
					}
				}
				if(optimal.aa == "optimize"){
					cat("Finished. Optimizing model parameters...", "\n")
					results.final <- nloptr(x0=log(ip.vector), eval_f = OptimizeEdgeLengthsGlobal, ub=upper.vector, lb=lower.vector, opts=opts, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, n.partitions=n.partitions, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=empirical.codon.freq.list, numcode=numcode, aa.properties=aa.properties, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
					print(results.final)
					cat("Finished. Finding optimal aa...", "\n")
					mle.pars.mat <- index.matrix
					mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
					aa.optim.list <- as.list(numeric(n.partitions))
					aa.optim.full.list <- as.list(numeric(n.partitions))
					for(partition.index in sequence(n.partitions)) {
						phy$edge.length = mle.pars.mat[partition.index, (max.par.model.count+1):ncol(mle.pars.mat)]
						gene.tmp <- read.dna(partitions[partition.index], format='fasta')
						if(!is.null(fasta.rows.to.keep)){
							gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
						}else{
							gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
						}
						codon.data <- DNAbinToCodonNumeric(gene.tmp)
						codon.data <- codon.data[phy$tip.label,]
						aa.optim.full.list[[partition.index]] = GetOptimalAAPerSite(x=log(mle.pars.mat[partition.index,1:max.par.model.count]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim.list[[partition.index]], root.p_array=empirical.codon.freq.list[[partition.index]], numcode=numcode, aa.properties=aa.properties, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, neglnl=TRUE)
						empirical.codon.freq.list[[partition.index]] <- CodonEquilibriumFrequencies(codon.data[,-1], aa.optim.full.list[[partition.index]], numcode=numcode)
						aa.optim.frame.to.add <- matrix(c("optimal", aa.optim.full.list[[partition.index]]), 1, dim(codon.data)[2])
						colnames(aa.optim.frame.to.add) <- colnames(codon.data)
						codon.data <- rbind(codon.data, aa.optim.frame.to.add)
						codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
						site.pattern.data.list[[partition.index]] = codon.data$unique.site.patterns
						site.pattern.count.list[[partition.index]] = codon.data$site.pattern.counts
						aa.optim.list[[partition.index]] = codon.data$optimal.aa		
					}
					cat("Finished. Final search given optimal aa...", "\n")
					results.final <- nloptr(x0=results.final$solution, eval_f = OptimizeEdgeLengthsGlobal, ub=upper.vector, lb=lower.vector, opts=opts, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, n.partitions=n.partitions, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=empirical.codon.freq.list, numcode=numcode, aa.properties=aa.properties, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
					cat("Finished. Summarizing results...", "\n")			
					
				}else{
					cat("Finished. Optimizing model parameters...", "\n")
					results.final <- nloptr(x0=log(ip.vector), eval_f = OptimizeEdgeLengthsGlobal, ub=upper.vector, lb=lower.vector, opts=opts, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, n.partitions=n.partitions, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=empirical.codon.freq.list, numcode=numcode, aa.properties=aa.properties, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores=n.cores, neglnl=TRUE)
					cat("Finished. Summarizing results...", "\n")
				}
			}
			mle.pars.mat <- index.matrix
			mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
			phy$edge.length <- apply(mle.pars.mat[,(max.par.model.count+1):ncol(mle.pars.mat)], 2, weighted.mean, w=nsites.vector)
		}
		loglik <- -(results.final$objective) #to go from neglnl to lnl
		#Counting parameters: Do we count the nsites too? Yup.
		np = max(index.matrix) + sum(nsites.vector)
		#Transforming parameters to obtain estimates of s across all partitions:
		#C.Phi.q.s <- mle.pars.mat[,1]
		C=2
		#Phi.q.s <- C.Phi.q.s / C
		Phi <- 0.5
		#q.s <- Phi.q.s / Phi
		q <- 4e-7
		#s <- q.s / q
		###########################
		obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, AICc = NULL, mle.pars=mle.pars.mat, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, aa.optim=aa.optim.full.list, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, ncats=ncats, k.levels=k.levels, aa.properties=aa.properties, empirical.codon.freqs=empirical.codon.freq.list, max.tol=max.tol) 
		class(obj) = "selac"		
	}
	return(obj)
}



######################################################################################################################################
######################################################################################################################################
### Print function for the selac class:
######################################################################################################################################
######################################################################################################################################

####This needs serious work####
print.selac <- function(x,...){
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,ntips,sum(x$nsites), x$k.levels, row.names="")
	names(output)<-c("-lnL","AIC", "ntax", "nsites", "k.levels")
	cat("\nFit\n")
	print(output)
	cpv.starting.parameters <- GetAADistanceStartingParameters(aa.properties=x$aa.properties)	
	if(x$aa.optim.type=="majrule" | x$aa.optim.type=="optimize"){
		if(x$nuc.model == "JC"){
			cat("\n")
			cat("\nSELAC Parameters\n")
			if(x$include.gamma==TRUE){
				output<-data.frame(x$mle.pars[1,1],x$mle.pars[1,2], x$mle.pars[1,3], x$mle.pars[1,4], x$mle.pars[1,9], row.names="")
				names(output)<-c("c","p","v","Ne","gamma")					
			}else{
				output<-data.frame(x$mle.pars[1,1],x$mle.pars[1,2], x$mle.pars[1,3], x$mle.pars[1,4], row.names="")
				names(output)<-c("c","p","v","Ne")				
			}
			print(output)
			cat("\n")
			cat("\nBase frequencies per partition\n")
			tmp <- cbind(x$mle.pars[,5:7], 1-colSums(t(x$mle.pars[,5:7])))
			base.freqs <- data.frame(t(tmp), row.names=c("A","C","G","T"))
			print(base.freqs)
		}else{
			cat("\n")
			cat("\nSELAC parameters\n")
			if(x$include.gamma==TRUE){
				output<-data.frame(x$mle.pars[1,1],x$mle.pars[1,2], x$mle.pars[1,3], x$mle.pars[1,4], x$mle.pars[1,8], row.names="")
				names(output)<-c("c","p","v","Ne","gamma")					
			}else{
				output<-data.frame(x$mle.pars[1,1],x$mle.pars[1,2], x$mle.pars[1,3], x$mle.pars[1,4], row.names="")
				names(output)<-c("c","p","v","Ne")				
			}
			print(output)
			cat("\n")
			tmp <- cbind(x$mle.pars[,5:7], 1-colSums(t(x$mle.pars[,5:7])))
			base.freqs <- data.frame(t(tmp), row.names=c("A","C","G","T"))
			print(base.freqs)
		}
		cat("\n")
	}
}


######################################################################################################################################
######################################################################################################################################
### Development NOTES -- JMB
######################################################################################################################################
######################################################################################################################################
