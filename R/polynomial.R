
######################################################################################################################################
######################################################################################################################################
### Simulation run
######################################################################################################################################
######################################################################################################################################

#Step 1: Simulate trait values.
#Step 2: Optimize model parameters.

res<-c()
for(i in 1:100){
    tmp.trait <- SimulateMonoPolyValue()
    res<-rbind(res, OptimizePolynomialVariables(tmp.trait, k=1)$k.1)
}
colMeans(res)


######################################################################################################################################
######################################################################################################################################
### If all else fails, let us simulate....
######################################################################################################################################
######################################################################################################################################

#To simulate probabilities assuming k=1, I believe we simply need to do the following:
#1/(1+e^-mi(theta), where mi(theta) = bk,0 + bk,1x + bk,2x^2 + bk,3x^3 + ... + ak,2kx^2k+1, where bk,0=xi, bk,j=ak,(j-1)/j, j = 1,2, ... , 2k+1.
#where you simply randomly choose values from a normal (or any other distribution). I ordered from smallest to largest -- looks right, but weird.

SimulateMonoPolyProbs <- function(nsize=100, xi=0, alpha=-0.5, beta =.1, k=1){
    #Note lambda is assumed to equal 1 automatically.
    theta = rnorm(nsize)
    phi = (alpha)^2 + beta
    coef.vec = CalculatePolynomialCoefficients(alpha.vec=alpha, phi.vec=phi, k=k)
    mi <- xi + (coef.vec[1,1]/2)*theta + (coef.vec[2,1]/3)*(theta^2) + (coef.vec[3,1]/4)*(theta^3)
    p.theta = 1 / (1+exp(-mi))
    return(p.theta[order(p.theta)])
}


#To simulate raw values assuming k=1, I believe we simply need to do the following:
#mi(theta) = bk,0 + bk,1x + bk,2x^2 + bk,3x^3 + ... + ak,2kx^2k+1, where bk,0=xi, bk,j=ak,(j-1)/j, j = 1,2, ... , 2k+1.
#where you simply randomly choose values from a normal (or any other distribution). I ordered from smallest to largest -- looks right, but weird.
SimulateMonoPolyValue <- function(nsize=100, xi=0, alpha=-0.5, beta=.1, k=1){
    #Note lambda is assumed to equal 1 automatically.
    theta = rnorm(nsize)
    phi = (alpha)^2 + beta
    coef.vec = CalculatePolynomialCoefficients(alpha.vec=alpha, phi.vec=phi, k=k)
    mi <- xi + (coef.vec[1,1]/2)*theta + (coef.vec[2,1]/3)*(theta^2) + (coef.vec[3,1]/4)*(theta^3)
    return(mi[order(mi)])
}

######################################################################################################################################
######################################################################################################################################
### Functions for calculating the polynomial 
######################################################################################################################################
######################################################################################################################################

#pg 30: "The optimal estimates for the monotonic polynomial will be searched for at each stage of k. Usually we start k=0, where the positive polynomial is a real positive constant lambda. JMB NOTE in our case k=0, lambda=1 -- makes things easier. Also note that beta, which is parameter 3, must be >0."
OptimizePolynomialVariables <- function(x, k){
	k.pars.list <- NULL
	k.sequence <- 1:k
    opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = "100000", "ftol_rel" = .Machine$double.eps^0.25)
    
    #Step 1: Optimize for k=1, assumes lambda=1. Parameter order: xik, alphak, phik:
    get.k1.pars <- nloptr(x0=c(0, -1/2, 1), eval_f = PolynomialLikelihoodFunction, ub=c(100, 100, 100), lb=c(-100, -100, 1e-100), opts=opts, x=x, k=k.sequence[1], k.minus.1.p=NULL, k.minus.2.p=NULL)
	k.pars.list$k.1 = get.k1.pars$solution
    k.pars.list$k.1.lik = -get.k1.pars$objective

    ###Not necessary yet:
    #Step 2: Optimize for k=2, if necessary:
    #if(!is.na(k.sequence[2] == 2)){
    #	get.k2.pars <- optim(c(1, 1, 2), fn=PolynomialLikelihoodFunction, method="L-BFGS-B", upper=c(Inf, Inf, Inf), lower=c(-Inf, -Inf, 1e-10), x=x, k=k.sequence[2], k.minus.1.p=k.pars.list$k.1, k.minus.2.p=NULL)
    #	k.pars.list$k.2 <- get.k2.pars$par
    #}
    #Step 3: Optimize for k=3, if necessary:
    #if(!is.na(k.sequence[3] == 3)){
    #	get.k3.pars <- optim(c(1, 1, 2), fn=PolynomialLikelihoodFunction, method="L-BFGS-B", upper=c(Inf, Inf, Inf), lower=c(-Inf, -Inf, 1e-10), x=x, k=k.sequence[3], k.minus.1.p=k.pars.list$k.1, k.minus.2.p=k.pars.list$k.2)
    #	k.pars.list$k.3 <- get.k3.pars$par
    #}
    #####################
    
	return(k.pars.list)
}


#pg 29: "Let fk(x) be the density function derived from Fk(x), i.e., fk(x) = h(m_2k_1(x)) * m'k(x). The ML estimates will be obtained by minimizing the negative log likelihood function F=-sum(log(fk(xj)))."
PolynomialLikelihoodFunction <- function(p, x, k, k.minus.1.p, k.minus.2.p){
	if(k == 1){
		coef.vec <- CalculatePolynomialCoefficients(alpha.vec=c(p[2]), phi.vec=c(p[3]), k)
	}
    
    ###Not necessary yet:
    #if(k == 2){
    #   coef.vec <- CalculatePolynomialCoefficients(alpha.vec=c(k.minus.1.p[2], p[2]), beta.vec=c(k.minus.1.p[3], p[3]), k)
    #}
    #if(k == 3){
    #	coef.vec <- CalculatePolynomialCoefficients(alpha.vec=c(k.minus.1.p[2], k.minus.2.p[2], p[2]), beta.vec=c(k.minus.1.p[3], k.minus.2.p[3], p[3]), k)
    #}
    ####################
    
    res = c()
	for(j in 1:length(x)){
        res <- c(res, MonotonicPolynomialTransform(x[j], p[1], coef.vec, k) * PositivePolynomialTransform(x[j], coef.vec, k))
    }
    
	if(any(res <= 0 )){
		return(10000000000)
	}else{
		return(sum(log(res)))
	}
}


#pg. 24 "Eq. 3.2.6: m'k = ak,0 + ak,1x + ak,2x^2 + ... + ak,2kx^2k."
PositivePolynomialTransform <- function(x, coef.vec, k){
	if(k == 0){
		mk <- 1
	}
	if(k == 1){
		mk <- coef.vec[1,1] + coef.vec[2,1]*x + coef.vec[3,1]*(x^2)
	}
    
    ###Not necessary yet:
    #if(k == 2){
    #	mk <- coef.vec[1,] + coef.vec[2,]*x + coef.vec[3,]*(x^2) + coef.vec[4,]*(x^3) + coef.vec[5,]*(x^4)
    #}
    #if(k == 3){
    #	mk <- coef.vec[1,] + coef.vec[2,]*x + coef.vec[3,]*(x^2) + coef.vec[4,]*(x^3) + coef.vec[5,]*(x^4) + coef.vec[6,]*(x^5) + coef.vec[7,]*(x^6)
    #}
    #####################
    
	return(mk)
}


#pg. 24 "Eq. 3.2.7: m2k_1 = bk,0 + bk,1x + bk,2x^2 + bk,3x^3 + ... + ak,2kx^2k+1, where bk,0=xi, bk,j=ak,(j-1)/j, j = 1,2, ... , 2k+1."
MonotonicPolynomialTransform <- function(x, xi, coef.vec, k){
    if(k == 0){
        mk2_1 <- xi + 1
    }
    if(k == 1){
		mk2_1 <- xi + (coef.vec[1,1]/2)*x + (coef.vec[2,1]/3)*(x^2) + (coef.vec[3,1]/4)*(x^3)
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
CalculatePolynomialCoefficients <- function(alpha.vec, phi.vec, k){
    #coef.mat is the matrix of coefficients:
	coef.vec <- 1
    #k.set is a sequence of k ending with the max k which is specified at the function call:
	k.set = 1:k
    #Matrix multiplication order does not matter for answer, but order matters for efficiency:
	for(k.index in 1:k){
		coef.vec <- CreatePolynomialMatrix(alpha.vec[k.index], phi.vec[k.index], k.set[k.index]) %*% coef.vec
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


######################################################################################################################################
######################################################################################################################################
### Relevant functions for incorporating this into selac
######################################################################################################################################
######################################################################################################################################

GetPairwiseProteinFixationProbabilitySingleSite <- function(d1, d2, s, C=2, Phi=0.5, q=4e-7, Ne=5e6){
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


FastCreateAllCodonFixationProbabilityMatrices <- function(s=0.01, aa.distances=CreateAADistanceMatrix(), C=2, Phi=0.5, q=4e-7, Ne=5e6, include.stop.codon=TRUE, numcode=1, flee.stop.codon.rate=0.9999999) {
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
						codon.fixation.probs[i,j, k] <- GetPairwiseProteinFixationProbabilitySingleSite(d1, d2, s=s, C=C, Phi=Phi, q=q, Ne=Ne)
					}else {
						if(s==0) { #handles stop codon case where neutral, so could possibly go into and out of stop codons
							codon.fixation.probs[i,j, k] <- 0
						}else {
							if(aa2!="*" && unique.aa[k]!="*") {
								codon.fixation.probs[i,j, k] <- 0 #JMB: Old = if we are somehow in a stop codon, have a very high rate of moving away from this; New = make is zero because in theory our model should use selection to kill these but infinite selection is rather harsh.
							}
						}
					}
				}
			}
		}
	}
	codon.fixation.probs[,,"*"] = 0 
	return(codon.fixation.probs)
}

