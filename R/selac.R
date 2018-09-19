######################################################################################################################################
######################################################################################################################################
### SELAC -- SELection on Amino acids and/or Codons
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu and Brian O

###LOAD REQUIRED PACKAGES -- eventually move to namespace:
## only set to TRUE when testing. Set to FALSE when committing changes
# if(FALSE){
#     library(ape)
#     library(expm)
#     library(nnet)
#     library(nloptr)
#     library(seqinr)
#     library(phangorn)
#     library(MASS)
#     library(parallel)
#     library(Rcpp)
#     library(RcppArmadillo)
#     library(inline)
#     library(deSolve)
#     ##load compiled library independent of working directory
#     wd <- getwd();
#     ##get last part of wd that ends in 'selac'
#     selac.dir <- regmatches(wd, regexpr(".*/selac", wd))
#     if(length(selac.dir) == 0) selac.dir <- "./selac"
#     so.locale <- paste(selac.dir, "/src/selacHMM.so",sep="")
#     dyn.load(so.locale)
#     rm(selac.dir, so.locale)
# }

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
### A collection of constants used by various functions
######################################################################################################################################
######################################################################################################################################


.codon.sets <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                        0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,
                        0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3), ncol=3)


.codon.set.translate <- matrix(c("a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "g", "g", "g", "g", "g", "g", "g", "g", "g", "g", "g", "g", "g", "g", "g", "g", "t", "t", "t", "t", "t", "t", "t", "t", "t", "t", "t", "t", "t", "t", "t", "t",
                                 "a", "a", "a", "a", "c", "c", "c", "c", "g", "g", "g", "g", "t", "t", "t", "t", "a", "a", "a", "a", "c", "c", "c", "c", "g", "g", "g", "g", "t", "t", "t", "t", "a", "a", "a", "a", "c", "c", "c", "c", "g", "g", "g", "g", "t", "t", "t", "t", "a", "a", "a", "a", "c", "c", "c", "c", "g", "g", "g", "g", "t", "t", "t", "t",
                                 "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t", "a", "c", "g", "t"), ncol=3)


.codon.name <- c("aaa" ,"aac" ,"aag" ,"aat" ,"aca" ,"acc" ,"acg" ,"act" ,"aga" ,"agc" ,"agg" ,"agt" ,"ata" ,"atc" ,"atg" ,"att" ,"caa" ,"cac" ,"cag" ,"cat", "cca" ,"ccc",
                 "ccg" ,"cct" ,"cga" ,"cgc" ,"cgg" ,"cgt" ,"cta" ,"ctc" ,"ctg" ,"ctt" ,"gaa" ,"gac" ,"gag" ,"gat" ,"gca" ,"gcc" ,"gcg" ,"gct" ,"gga" ,"ggc", "ggg" ,"ggt",
                 "gta" ,"gtc" ,"gtg" ,"gtt" ,"taa" ,"tac" ,"tag" ,"tat" ,"tca" ,"tcc" ,"tcg" ,"tct" ,"tga" ,"tgc" ,"tgg" ,"tgt" ,"tta" ,"ttc" ,"ttg" ,"ttt")


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


.unique.aa <- c("K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E", "D", "A", "G", "V", "*", "Y", "C", "W", "F")

.numcode.translation.idx <- c(1:6, 1, 1, 7:14, rep(1,5), 15:17)


######################################################################################################################################
######################################################################################################################################
### Various functions used by main function:
######################################################################################################################################
######################################################################################################################################

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
  if(model == "HKY") {
      index <- matrix(NA, 4, 4)
      rates <- c(1,rates)
      sel <- col(index) < row(index)
      index[sel] <- c(1,2,1,1,2,1)
      index <- t(index)
      index[sel] <- c(1,2,1,1,2,1)
      nuc.mutation.rates <- matrix(0, nrow=4, ncol=4)
      nuc.mutation.rates<-matrix(rates[index], dim(index))
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
    nuc.mutation.rates[4,3] <- nuc.mutation.rates[3,4] <- 1
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
    np <- 12
    index[col(index) != row(index)] <- 1:np
    nuc.mutation.rates <- matrix(0, nrow=4, ncol=4)
    nuc.mutation.rates<-matrix(rates[index], dim(index))
    rownames(nuc.mutation.rates) <- n2s(0:3)
    colnames(nuc.mutation.rates) <- n2s(0:3)
    nuc.mutation.rates[3,4] = 1
    diag(nuc.mutation.rates) <- 0
    diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
    #Next we take our rates and find the homogeneous solution to Q*pi=0 to determine the base freqs:
    base.freqs <- Null(nuc.mutation.rates)
    #Rescale base.freqs so that they sum to 1:
    base.freqs.scaled <- c(base.freqs/sum(base.freqs))
    base.freqs.scaled.matrix <- rbind(base.freqs.scaled, base.freqs.scaled, base.freqs.scaled, base.freqs.scaled)
    diag(nuc.mutation.rates) <- 0
    #Rescale Q to account for base.freqs:
    nuc.mutation.rates <- nuc.mutation.rates * base.freqs.scaled.matrix
    diag(nuc.mutation.rates) <- -rowSums(nuc.mutation.rates)
    return(nuc.mutation.rates)
  }
}


CreateCodonMutationMatrixIndex <- function() {
  nuc.rates.index = matrix(1:16, 4, 4)
  #codon.sets <- CreateCodonSets()
  n.codons <- dim(.codon.sets)[1]
  codon.mutation.rates <- matrix(data=0, nrow=n.codons, ncol=n.codons)
  rownames(codon.mutation.rates) <- rep("",n.codons)
  colnames(codon.mutation.rates) <- rep("",n.codons)
  for (i in sequence(n.codons)) {
    for (j in sequence(n.codons)) {
      if(sum(.codon.sets[i,] == .codon.sets[j,])==2) { #means that two of the bases match
        mismatch.position <- which(.codon.sets[i,] != .codon.sets[j,])
        codon.mutation.rates[i,j] <- nuc.rates.index[1+.codon.sets[i,mismatch.position], 1+.codon.sets[j, mismatch.position]] #nucs numbered from 0:3, rows are 1:4, thus the add 1
      }
    }
    #codon.name <- paste(n2s(as.numeric(.codon.sets[i,])), collapse="")
    rownames(codon.mutation.rates)[i] <- .codon.name[i]
    colnames(codon.mutation.rates)[i] <- .codon.name[i]
  }
  codon.mutation.rates[codon.mutation.rates==0] = NA
  return(codon.mutation.rates)
}


CreateCodonMutationMatrixIndexEvolveAA <- function() {
  codon.mut.index.mat <- CreateCodonMutationMatrixIndex()

  n.codons <- dim(.codon.sets)[1]
  mat.dim <- 21*n.codons
  evolv.codon.mut.index.mat <- matrix(data=0, nrow=mat.dim, ncol=mat.dim)

  for(i in 1:21)
  {
    index.vec.diag <- (1+(i-1)*64):(64+(i-1)*64)
    # fill in a complete coloumn of matrices, has to be done before the matrix on the diagonal
    evolv.codon.mut.index.mat[index.vec.diag, -index.vec.diag] <- diag(17, 64, 64)#do.call(cbind, replicate(20, diag(18, 64, 64), simplify = FALSE))
    evolv.codon.mut.index.mat[index.vec.diag, index.vec.diag] <- codon.mut.index.mat # fill in matrix on the diagonal
  }

  evolv.codon.mut.index.mat[evolv.codon.mut.index.mat == 0] = NA
  return(evolv.codon.mut.index.mat)
}


CreateCodonMutationMatrix <- function(nuc.mutation.rates) {
  #codon.sets <- CreateCodonSets()
  n.codons <- dim(.codon.sets)[1]
  codon.mutation.rates <- matrix(data=0, nrow=n.codons, ncol=n.codons)
  rownames(codon.mutation.rates) <- rep("",n.codons)
  colnames(codon.mutation.rates) <- rep("",n.codons)
  for (i in sequence(n.codons)) {
    for (j in sequence(n.codons)) {
      if(sum(.codon.sets[i,] == .codon.sets[j,])==2) { #means that two of the bases match
        mismatch.position <- which(.codon.sets[i,] != .codon.sets[j,])
        codon.mutation.rates[i,j] <- nuc.mutation.rates[1+.codon.sets[i,mismatch.position], 1+.codon.sets[j, mismatch.position]] #nucs numbered from 0:3, rows are 1:4, thus the add 1
      }
    }
    #codon.name <- paste(n2s(as.numeric(.codon.sets[i,])), collapse="")
    rownames(codon.mutation.rates)[i] <- .codon.name
    colnames(codon.mutation.rates)[i] <- .codon.name

  }
  diag(codon.mutation.rates) <- 0
  diag(codon.mutation.rates) <- -rowSums(codon.mutation.rates)
  rownames(codon.mutation.rates) <- colnames(codon.mutation.rates) <- .codon.name
  return(codon.mutation.rates)
}


CreateCodonMutationMatrixMutSel <- function(omega.par, fitness.pars, nuc.mutation.rates, numcode) {
    #codon.sets <- CreateCodonSets()
    n.codons <- dim(.codon.sets)[1]
    codon.mutation.rates <- matrix(data=0, nrow=n.codons, ncol=n.codons)
    rownames(codon.mutation.rates) <- rep("",n.codons)
    colnames(codon.mutation.rates) <- rep("",n.codons)
    #codon.set.translate <- apply(.codon.sets, 2, n2s)
    #codon.name <- apply(.codon.set.translate, 1, paste, collapse="")
    aa.translations <- .aa.translation[[numcode]][.codon.name]
    for (i in sequence(n.codons)) {
        for (j in sequence(n.codons)) {
            if(aa.translations[i] == aa.translations[j]){ #synonymous
                if(sum(.codon.sets[i,] == .codon.sets[j,])==2) { #means that two of the bases match
                    mismatch.position <- which(.codon.sets[i,] != .codon.sets[j,])
                    matched.position <- which(.codon.sets[i,] == .codon.sets[j,])
                    if((fitness.pars[j]-fitness.pars[i]) == 0){
                        codon.mutation.rates[i,j] = nuc.mutation.rates[1+.codon.sets[i,mismatch.position], 1+.codon.sets[j, mismatch.position]]
                    }else{
                        codon.mutation.rates[i,j] <- nuc.mutation.rates[1+.codon.sets[i,mismatch.position], 1+.codon.sets[j, mismatch.position]] * ((fitness.pars[j] - fitness.pars[i]) / (1-exp(fitness.pars[i] - fitness.pars[j])))
                    }
                }
            }else{ #nonsynonymous
                if(sum(.codon.sets[i,] == .codon.sets[j,])==2) { #means that two of the bases match
                    mismatch.position <- which(.codon.sets[i,] != .codon.sets[j,])
                    matched.position <- which(.codon.sets[i,] == .codon.sets[j,])
                    if((fitness.pars[j]-fitness.pars[i]) == 0){
                        codon.mutation.rates[i,j] = omega.par * nuc.mutation.rates[1+.codon.sets[i,mismatch.position], 1+.codon.sets[j, mismatch.position]]
                    }else{
                        codon.mutation.rates[i,j] <- omega.par * nuc.mutation.rates[1+.codon.sets[i,mismatch.position], 1+.codon.sets[j, mismatch.position]] * ((fitness.pars[j] - fitness.pars[i]) / (1-exp(fitness.pars[i] - fitness.pars[j])))
                    }
                }
            }
        }
    }
    #Remove stop codon rates -- they should be removed already, but just in case...
    codon.mutation.rates[which(aa.translations == "*"),] = codon.mutation.rates[,which(aa.translations == "*")] = 0
    #Now let us finish up the matrix:
    rownames(codon.mutation.rates) <- colnames(codon.mutation.rates) <- .codon.name
    diag(codon.mutation.rates) <- 0
    diag(codon.mutation.rates) <- -rowSums(codon.mutation.rates)
    return(codon.mutation.rates)
}


CreateCodonMutationMatrixYN98 <- function(x, codon.freqs, numcode) {
    omega.par = x[1]
    kappa.par <- x[2]
    #The last value is arbitrarily set to 0 per Yang and Nielsen (2008):
    #codon.sets <- CreateCodonSets()
    n.codons <- dim(.codon.sets)[1]
    codon.mutation.rates <- matrix(data=0, nrow=n.codons, ncol=n.codons)
    rownames(codon.mutation.rates) <- rep("",n.codons)
    colnames(codon.mutation.rates) <- rep("",n.codons)
    #codon.set.translate <- apply(.codon.sets, 2, n2s)
    #codon.name <- apply(.codon.set.translate, 1, paste, collapse="")
    #We add this in because the stop codons are not included in Grantham's distance calculation:
    aa.translations <- .aa.translation[[numcode]][.codon.name]
    for (i in sequence(n.codons)) {
        for (j in sequence(n.codons)) {
            if(aa.translations[i] == aa.translations[j]){ #synonymous -- set distance to zero.
                if(sum(.codon.sets[i,] == .codon.sets[j,])==2) { #means that two of the bases match
                    mismatch.position <- which(.codon.sets[i,] != .codon.sets[j,])
                    matched.position <- which(.codon.sets[i,] == .codon.sets[j,])
                    if(.codon.sets[i, mismatch.position] == 0 & .codon.sets[j,mismatch.position] == 2 | .codon.sets[i, mismatch.position] == 2 & .codon.sets[j,mismatch.position] == 0 | .codon.sets[i, mismatch.position] == 1 & .codon.sets[j,mismatch.position] == 3 | .codon.sets[i, mismatch.position] == 3 & .codon.sets[j,mismatch.position] == 1){
                        codon.mutation.rates[i,j] <- kappa.par * codon.freqs[j]
                    }else{
                        codon.mutation.rates[i,j] <- codon.freqs[j]
                    }
                }
            }else{ #nonsynonymous -- so we need to know Grantham's distance.
                if(sum(.codon.sets[i,] == .codon.sets[j,])==2) { #means that two of the bases match
                    mismatch.position <- which(.codon.sets[i,] != .codon.sets[j,])
                    matched.position <- which(.codon.sets[i,] == .codon.sets[j,])
                    if(.codon.sets[i, mismatch.position] == 0 & .codon.sets[j,mismatch.position] == 2 | .codon.sets[i, mismatch.position] == 2 & .codon.sets[j,mismatch.position] == 0 | .codon.sets[i, mismatch.position] == 1 & .codon.sets[j,mismatch.position] == 3 | .codon.sets[i, mismatch.position] == 3 & .codon.sets[j,mismatch.position] == 1){
                        codon.mutation.rates[i,j] <- kappa.par * codon.freqs[j] * omega.par
                    }else{
                        codon.mutation.rates[i,j] <- codon.freqs[j] * omega.par
                    }
                }
            }
        }
    }
    #Remove stop codon rates -- they should be removed already, but just in case...
    codon.mutation.rates[which(aa.translations == "*"),] = codon.mutation.rates[,which(aa.translations == "*")] = 0
    #Now let us finish up the matrix:
    rownames(codon.mutation.rates) <- colnames(codon.mutation.rates) <- .codon.name
    diag(codon.mutation.rates) <- 0
    diag(codon.mutation.rates) <- -rowSums(codon.mutation.rates)
    return(codon.mutation.rates)
}



CreateCodonMutationMatrixGY94 <- function(x, aa.distances, codon.freqs, numcode) {
  v.par <- x[1]
  kappa.par = x[2]
  #The last value is arbitrarily set to 0 per Yang and Nielsen (2008):
  #codon.sets <- CreateCodonSets()
  n.codons <- dim(.codon.sets)[1]
  codon.mutation.rates <- matrix(data=0, nrow=n.codons, ncol=n.codons)
  rownames(codon.mutation.rates) <- rep("",n.codons)
  colnames(codon.mutation.rates) <- rep("",n.codons)
  #codon.set.translate <- apply(.codon.sets, 2, n2s)
  #codon.name <- apply(.codon.set.translate, 1, paste, collapse="")
  #We add this in because the stop codons are not included in Grantham's distance calculation:
  aa.distances <- rbind(aa.distances, "*"=0, deparse.level=2)
  aa.distances <- cbind(aa.distances, "*"=0, deparse.level=2)
  aa.translations <- .aa.translation[[numcode]][.codon.name]
  for (i in sequence(n.codons)) {
    for (j in sequence(n.codons)) {
      if(aa.translations[i] == aa.translations[j]){ #synonymous -- set distance to zero.
        if(sum(.codon.sets[i,] == .codon.sets[j,])==2) { #means that two of the bases match
          mismatch.position <- which(.codon.sets[i,] != .codon.sets[j,])
          matched.position <- which(.codon.sets[i,] == .codon.sets[j,])
          if(.codon.sets[i, mismatch.position] == 0 & .codon.sets[j,mismatch.position] == 2 | .codon.sets[i, mismatch.position] == 2 & .codon.sets[j,mismatch.position] == 0 | .codon.sets[i, mismatch.position] == 1 & .codon.sets[j,mismatch.position] == 3 | .codon.sets[i, mismatch.position] == 3 & .codon.sets[j,mismatch.position] == 1){
            codon.mutation.rates[i,j] <- kappa.par * codon.freqs[j] * exp(-0/v.par)
          }else{
            codon.mutation.rates[i,j] <- codon.freqs[j] * exp(-0/v.par)
          }
        }
      }else{ #nonsynonymous -- so we need to know Grantham's distance.
        if(sum(.codon.sets[i,] == .codon.sets[j,])==2) { #means that two of the bases match
          mismatch.position <- which(.codon.sets[i,] != .codon.sets[j,])
          matched.position <- which(.codon.sets[i,] == .codon.sets[j,])
          if(.codon.sets[i, mismatch.position] == 0 & .codon.sets[j,mismatch.position] == 2 | .codon.sets[i, mismatch.position] == 2 & .codon.sets[j,mismatch.position] == 0 | .codon.sets[i, mismatch.position] == 1 & .codon.sets[j,mismatch.position] == 3 | .codon.sets[i, mismatch.position] == 3 & .codon.sets[j,mismatch.position] == 1){
            codon.mutation.rates[i,j] <- kappa.par * codon.freqs[j] * exp(-aa.distances[aa.translations[i], aa.translations[j]]/v.par)
          }else{
            codon.mutation.rates[i,j] <- codon.freqs[j] * exp(-aa.distances[aa.translations[i], aa.translations[j]]/v.par)
          }
        }
      }
    }
  }
  #Remove stop codon rates -- they should be removed already, but just in case...
  codon.mutation.rates[which(aa.translations == "*"),] = codon.mutation.rates[,which(aa.translations == "*")] = 0
  #Now let us finish up the matrix:
  rownames(codon.mutation.rates) <- colnames(codon.mutation.rates) <- .codon.name
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


CodonStringToCharacter <- function(x) { #remember that codon numbers start at 1
  triplet <- x
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


NucleotideStringToCharacter <- function(x) { #remember that codon numbers start at 1
  singlet <- x
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


GetPairwiseProteinFixationProbabilityArbitraryLength <- function(protein1, protein2, protein_op, aa.distances, nsites, C=4, Phi=0.5, q=4e-7, Ne=5e6){
  d1 <- GetProteinProteinDistance(protein1,protein_op,aa.distances)
  d2 <- GetProteinProteinDistance(protein2,protein_op,aa.distances)
  if(length(d1)!=length(d2)) #throw error if length of proteins are not the same
    stop("error: 2 proteins are of different lengths!")
  if(length(d1)==1){ #only one amino acid
    return(GetPairwiseProteinFixationProbabilitySingleSite(d1, d2, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne))
  }
  else{
    if((length(d1)!=1)) #if s is given as a scalar, then treat it to be the same across all sites
      l = length(d1)
    cmp = CompareVectors(d1,d2)
    if(cmp$num > 1) return(0) #more than 1 position differ
    else if((cmp$num ==0)) return(1/(2*Ne)) #same fitness/functionality
    else{ #exactly 1 position differs
      pos = cmp$pos
      return(GetPairwiseProteinFixationProbabilitySingleSite(d1[pos], d2[pos], nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne))
    }
  }
}


GetFitness <- function(focal.protein, optimal.protein, aa.distances, nsites, C=4, Phi=0.5, q=4e-7) {
  focal.d <- GetProteinProteinDistance(focal.protein, optimal.protein, aa.distances)
  optimal.d <- GetProteinProteinDistance(optimal.protein, optimal.protein, aa.distances)
  return(exp(-(C+(C/nsites))*Phi*q*(focal.d-optimal.d)))
}


GetPairwiseProteinFixationProbabilitySingleSite <- function(d1, d2, nsites, C=4, Phi=0.5, q=4e-7, Ne=5e6, diploid=TRUE){
  if(diploid==TRUE){
    b = 1
  }else{
    b = 2
  }
  if(d1==d2){ #When the fitnesses are the same, neutral case, pure drift
    return(1/(2*Ne))
  }else{
    fit_ratio <- exp(-(C+(C/nsites))*Phi*q*(d1-d2)) #f1/f2
    if(fit_ratio==Inf){ #1 is much better than 2 (the mutant)
      return(0)
    }
    else{
      if(fit_ratio==1){
        return(1/(2*Ne))
      }
      else{
        return((1-fit_ratio^b)/(1-fit_ratio^(2*Ne)))
      }
    }
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
  #d <- sapply(c(1:length(protein1)),site_d,simplify=TRUE)
  if(length(protein1) == 1){
    d <- site_d(1)
  }else{
    d <- sapply(1:length(protein1), site_d, simplify=TRUE)
  }
  return(d)
}


#FastCreateAllCodonFixationProbabilityMatrices <- function(aa.distances=CreateAADistanceMatrix(), nsites, C=2, Phi=0.5, q=4e-7, Ne=5e6, include.stop.codon=TRUE, numcode=1, diploid=TRUE, flee.stop.codon.rate=0.9999999) {
##	#codon.sets <- CreateCodonSets()
#	.codon.sets <- expand.grid(0:3, 0:3, 0:3)
#	.codon.sets <- data.frame(first=.codon.sets[,3], second=.codon.sets[,2], third=.codon.sets[,1]) #reordering to group similar codons
#	n.codons <- dim(.codon.sets)[1]
#	codon.names <- rep("", n.codons)
#	for (i in sequence(n.codons)) {
#		codon.names[i] <- paste(n2s(as.numeric(.codon.sets[i,])), collapse="")
#	}
#	codon.aa <- sapply(codon.names, TranslateCodon, numcode=numcode)
#	unique.aa <- unique(codon.aa)
#	codon.fixation.probs <- array(data=0, dim=c(n.codons, n.codons, length(unique.aa)), dimnames=list(codon.names, codon.names, unique.aa))
#	for (i in sequence(n.codons)) {
#		for (j in sequence(n.codons)) {
#			if(sum(.codon.sets[i,] == .codon.sets[j,])>=2) { #match at two or three sites of three
#				for (k in sequence(length(unique.aa))) {
#					aa1 <- codon.aa[i]
#					aa2 <- codon.aa[j]
#					if(aa1!="*" && aa2!="*" && unique.aa[k]!="*") { #says we cannot mutate to stop codons and stop codons can never be optimal
#						d1 <- GetProteinProteinDistance(protein1=aa1, protein2=unique.aa[k], aa.distances=aa.distances)
#						d2 <- GetProteinProteinDistance(protein1=aa2, protein2=unique.aa[k], aa.distances=aa.distances)
#						codon.fixation.probs[i,j, k] <- GetPairwiseProteinFixationProbabilitySingleSite(d1, d2, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, diploid=diploid)
#					}else {
##						We have dropped s from the model as it is now explained through grantham like distances:
##						if(s==0) { #handles stop codon case where neutral, so could possibly go into and out of stop codons
##							codon.fixation.probs[i,j, k] <- 0
##						}else {
#							if(aa2!="*" && unique.aa[k]!="*") {
#								codon.fixation.probs[i,j, k] <- 0 #Old = if we are somehow in a stop codon, have a very high rate of moving away from this; New = make is zero because in theory our model should use selection to kill these but infinite selection is rather harsh.
##							}
#						}
#					}
#				}
#			}
#		}
#	}
#	codon.fixation.probs[,,"*"] = 0
#	return(codon.fixation.probs)
#}


CreateAAFixationMatrixForEverything <- function(aa.distances=CreateAADistanceMatrix(), nsites, C=4, Phi=0.5, q=4e-7, Ne=5e6, include.stop.codon=TRUE, numcode=1, diploid=TRUE) {
  states <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  fixation.array <- array(data=0, dim=rep(length(states)+include.stop.codon,3)) #adding the boolean to leave space for a stop codon if needed
  for (row.index in sequence(length(states))) {
    for (col.index in sequence(length(states))) {
      for (optimal.index in sequence(length(states))) {
        fixation.array[row.index, col.index, optimal.index] <- GetPairwiseProteinFixationProbabilitySingleSite(GetProteinProteinDistance(protein1=states[row.index], protein2=states[optimal.index], aa.distances=aa.distances),  GetProteinProteinDistance(protein1=states[col.index], protein2=states[optimal.index], aa.distances=aa.distances), nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, diploid=diploid)

      }
    }
  }
  if(include.stop.codon) {
    states <- c(states, "*")
  }
  dimnames(fixation.array) <- list(states, states, states)
  return(fixation.array)
}


FastCreateAllCodonFixationProbabilityMatrices <- function(aa.distances=CreateAADistanceMatrix(), nsites, C=4, Phi=0.5, q=4e-7, Ne=5e6, include.stop.codon=TRUE, numcode=1, diploid=TRUE, flee.stop.codon.rate=0.9999999) {
  #codon.sets <- CreateCodonSets()
  #codon.sets <- expand.grid(0:3, 0:3, 0:3)
  #codon.sets[,c(3,2,1)] <- .codon.sets[,c(1,2,3)] #re-ordering as in the original one
  colnames(.codon.sets) <- c("first", "second", "third")
  n.codons <- dim(.codon.sets)[1]
  codon.names <- rep("", n.codons)
  aa.fixation.probs <- CreateAAFixationMatrixForEverything(aa.distances=aa.distances, nsites, C, Phi, q, Ne, include.stop.codon, numcode, diploid)
  for (i in sequence(n.codons)) {
    codon.names[i] <- paste(n2s(as.numeric(.codon.sets[i,])), collapse="")
  }
  codon.aa <- sapply(codon.names, TranslateCodon, numcode=numcode)
  #unique.aa <- unique(codon.aa)

  codon.fixation.probs <- array(data=0, dim=c(n.codons, n.codons, length(.unique.aa)), dimnames=list(codon.names, codon.names, .unique.aa))
  for (i in sequence(n.codons)) {
    for (j in sequence(n.codons)) {
      if(sum(.codon.sets[i,] == .codon.sets[j,])>1 ) { #match at two or more sites
        for (k in sequence(length(.unique.aa))) {
          aa1 <- codon.aa[i]
          aa2 <- codon.aa[j]
          codon.fixation.probs[i,j, k] <- aa.fixation.probs[aa1, aa2, .unique.aa[k]]
        }
      }
    }
  }
  codon.fixation.probs[,,"*"] = 0
  return(codon.fixation.probs)
}


FastCreateOptAATransitionMatrices <- function(aa.distances=CreateAADistanceMatrix(), C, Phi, q, Ne, diploid, numcode=1, importance = 1) { #Cedric: added importance

  if(diploid == TRUE) {
    Ne <- 2*Ne
  } #Cedric: pay attention to diploid flag

  aa.dist.names <- colnames(aa.distances)
  aa.distances <- cbind(aa.distances, 0)
  aa.distances <- rbind(aa.distances, 0)
  colnames(aa.distances) <- c(aa.dist.names, "*")
  rownames(aa.distances) <- c(aa.dist.names, "*")

  numcode.idx <- .numcode.translation.idx[numcode]
  aa.names <- .aa.translation[[numcode.idx]]

  #aa.trans.mat <- (1.0/(aa.distances[.unique.aa, .unique.aa])^importance)/Ne #Cedric: adjusting for imporatance parameter and using 1/d instead of d
  aa.trans.mat <- (exp(-importance*aa.distances[.unique.aa, .unique.aa]))/Ne
  diag(aa.trans.mat) <- 0 # because R CAN divide by 0, some real Chuck Norris stuff here
  aa.trans.mat[,colnames(aa.trans.mat) == "*"] <- 0 # find better solution
  aa.trans.mat[colnames(aa.trans.mat) == "*",] <- 0

  # normalize distance
  norm.const <- rowSums(aa.trans.mat, na.rm = T)
  for(i in nrow(aa.trans.mat))
    aa.trans.mat[i,] <- aa.trans.mat[i,] / norm.const[i]


  aa.trans.matrices <- vector("list", 21)
  for(j in 1:21) {
    trans.matrix <- matrix(0, ncol=1344, nrow=64)
    for(i in 1:21) {
      index.vec.diag <- (1+(i-1)*64):(64+(i-1)*64)
      trans.matrix[, index.vec.diag] <- diag(aa.trans.mat[.unique.aa[j], .unique.aa[i]], ncol=64, nrow=64)
    }
    aa.trans.matrices[[j]] <- trans.matrix
  }
  names(aa.trans.matrices) <- .unique.aa
  return(aa.trans.matrices)
}


FastCreateEvolveAACodonFixationProbabilityMatrix <- function(aa.distances = CreateAADistanceMatrix(), nsites, C = 4, Phi = 0.5, q = 4e-7, Ne = 5e6, include.stop.codon = TRUE, numcode = 1, diploid = TRUE, flee.stop.codon.rate = 0.9999999, importance = 1) { #Cedric: Added the importance parameter
  codon.fixation.probs <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances, nsites, C, Phi, q, Ne, include.stop.codon, numcode, diploid, flee.stop.codon.rate)
  opt.aa.transition.rate <- FastCreateOptAATransitionMatrices(aa.distances=aa.distances, C=C, Phi=Phi, q=q, Ne=Ne, diploid=diploid, numcode=numcode, importance) #Cedric: passing importance through

  n.codons <- dim(.codon.sets)[1]
  mat.dim <- 21*n.codons
  evolv.codon.fixation.probs <- matrix(data=0, nrow=mat.dim, ncol=mat.dim)

  for(i in 1:21) {
    index.vec.diag <- (1+(i-1)*64):(64+(i-1)*64)
    evolv.codon.fixation.probs[index.vec.diag, ] <- opt.aa.transition.rate[[i]] # fill in a complete coloumn of matrices, has to be done before the matrix on the diagonal
    evolv.codon.fixation.probs[index.vec.diag, index.vec.diag] <- codon.fixation.probs[,,.unique.aa[i]] # fill in matrix on the diagonal
  }

  rownames(evolv.codon.fixation.probs) <- paste(rep(.codon.name, times=21), rep(.unique.aa, each=64), sep="")
  colnames(evolv.codon.fixation.probs) <- paste(rep(.codon.name, times=21), rep(.unique.aa, each=64), sep="")

  return(evolv.codon.fixation.probs)
}


## Work in progress ##
#cppFunction('NumericMatrix CreateCodonFixationProbabilityMatrixGivenOptimalAA(NumericMatrix codon_sets, StringVector aa, NumericMatrix aa_distances, int nsites, double C, double Phi, double q, double Ne, bool include_stop_codon,  int numcode, bool diploid, Function GetProteinProteinDistance, Function GetPairwiseProteinFixationProbabilitySingleSite, int optimalaa_offset0index) {
#	NumericMatrix codon_fixation_probs_aa(codon_sets.nrow(), codon_sets.nrow());
#//	Rcpp::Rcout << "ok";
#//	Rcout << aa;
#	for (int i=0; i<codon_sets.nrow(); ++i) {
#		for(int j=0; j<codon_sets.nrow(); ++j) {
#			int mismatches=0;
#			for(int pos=0; pos<3; ++pos) {
#				if(codon_sets(i, pos)==codon_sets(j, pos)) {
#					mismatches++;
#				}
#			}
#			if(mismatches<2) {
#				StringVector proteinI(1);
#				StringVector proteinJ(1);
#				StringVector proteinOptimal(1);
#				proteinI(0) = aa(i);
#				proteinJ(0) = aa(j);
#				proteinOptimal(0) = aa(optimalaa_offset0index);
#				codon_fixation_probs_aa(i,j) = as<double>(GetPairwiseProteinFixationProbabilitySingleSite(GetProteinProteinDistance(_["protein1"]=proteinI, _["protein2"]=proteinOptimal, _["aa.distances"]=aa_distances),GetProteinProteinDistance(_["protein1"]=proteinJ, _["protein2"]=proteinOptimal, _["aa.distances"]=aa_distances), _["nsites"]=nsites, _["C"]=C, _["Phi"]=Phi, _["q"]=q, _["Ne"]=Ne, _["diploid"]=diploid));
#				//codon_fixation_probs_aa(i,j) = as<double>(GetPairwiseProteinFixationProbabilitySingleSite(GetProteinProteinDistance(_["protein1"]=aa(i), _["protein2"]=aa(optimalaa_offset0index), _["aa.distances"]=aa_distances),  GetProteinProteinDistance(_["protein1"]=aa(j), _["protein2"]=aa(optimalaa_offset0index), _["aa.distances"]=aa_distances), _["nsites"]=nsites, _["C"]=C, _["Phi"]=Phi, _["q"]=q, _["Ne"]=Ne, _["diploid"]=diploid));
#			}
#		}
#	}
#	return(codon_fixation_probs_aa);
#}')
#result <- CreateCodonFixationProbabilityMatrixGivenOptimalAA(.codon.sets, unname(codon.aa),  aa.distances, nsites, C, Phi, q, Ne, include.stop.codon,  numcode, diploid, GetProteinProteinDistance, GetPairwiseProteinFixationProbabilitySingleSite, k-1)

## Work in progress ##
#CCQuestionablyFastCreateAllCodonFixationProbabilityMatrices <- function(aa.distances=CreateAADistanceMatrix(), nsites, C=4.0, Phi=0.5, q=4e-7, Ne=5e6, include.stop.codon=TRUE, numcode=1, diploid=TRUE, flee.stop.codon.rate=0.9999999) {
#	.codon.sets <- CreateCodonSets()
#	.codon.sets <- expand.grid(0:3, 0:3, 0:3)
#	.codon.sets[,c(3,2,1)] <- .codon.sets[,c(1,2,3)] #re-ordering as in the original one
#	colnames(.codon.sets) <- c("first", "second", "third")
#	n.codons <- dim(.codon.sets)[1]
#	codon.names <- rep("", n.codons)
#	aa.fixation.probs <- CreateAAFixationMatrixForEverything(aa.distances=aa.distances, nsites, C, Phi, q, Ne, include.stop.codon, numcode, diploid)
#	for (i in sequence(n.codons)) {
#		codon.names[i] <- paste(n2s(as.numeric(.codon.sets[i,])), collapse="")
#	}
#	codon.aa <- sapply(codon.names, TranslateCodon, numcode=numcode)
#	unique.aa <- unique(codon.aa)

#	codon.fixation.probs <- array(data=0, dim=c(n.codons, n.codons, length(unique.aa)), dimnames=list(codon.names, codon.names, unique.aa))
#	for (k in sequence(length(unique.aa))) {
#		codon.fixation.probs[,,k] <- CreateCodonFixationProbabilityMatrixGivenOptimalAA(.codon.sets, unname(codon.aa),  aa.distances, nsites, C, Phi, q, Ne, include.stop.codon,  numcode, diploid, GetProteinProteinDistance, GetPairwiseProteinFixationProbabilitySingleSite, k-1) #k-1 due to C++ counting from zero
#	}
#	return(codon.fixation.probs)
#}


DiagArray <- function (dim){
  n <- dim[2]
  d <- seq(1, n*n, by=n+1)
  as.vector(outer(d, seq(0, by=n*n, length=prod(dim[-1:-2])), "+"))
}


FastCreateAllCodonFixationProbabilityMatricesSetToOne <- function(numcode=1) {
  #   codon.sets <- CreateCodonSets()
  #	.codon.sets <- expand.grid(0:3, 0:3, 0:3)
  #	.codon.sets <- data.frame(first=.codon.sets[,3], second=.codon.sets[,2], third=.codon.sets[,1]) #reordering to group similar codons
  n.codons <- dim(.codon.sets)[1]
  codon.names <- rep("", n.codons)
  for (i in sequence(n.codons)) {
    codon.names[i] <- paste(n2s(as.numeric(.codon.sets[i,])), collapse="")
  }
  codon.aa <- sapply(codon.names, TranslateCodon, numcode=numcode)
  #unique.aa <- unique(codon.aa)
  codon.fixation.rates <- array(data=1, dim=c(n.codons, n.codons, length(.unique.aa)), dimnames=list(codon.names, codon.names, .unique.aa))
  return(codon.fixation.rates)
}


CreateCodonFixationProbabilityMatrix <- function(aa_op, s, aa.distances, nsites, C=4, Phi=0.5, q=4e-7, Ne=5e6, include.stop.codon=TRUE, numcode=1){
  #   codon.sets <- CreateCodonSets()
  #	.codon.sets <- expand.grid(0:3, 0:3, 0:3)
  #	.codon.sets <- data.frame(first=.codon.sets[,3], second=.codon.sets[,2], third=.codon.sets[,1]) #reordering to group similar codons
  n.codons <- dim(.codon.sets)[1]
  codon.fixation.rates <- matrix(data=0, nrow=n.codons, ncol=n.codons)
  codon.names <- rep("", n.codons)
  for (i in sequence(n.codons)) {
    codon.names[i] <- paste(n2s(as.numeric(.codon.sets[i,])), collapse="")
  }
  rownames(codon.fixation.rates) <- codon.names
  colnames(codon.fixation.rates) <- codon.names
  codon.aa <- sapply(codon.names, TranslateCodon, numcode=numcode)

  for (i in sequence(n.codons)) {
    for (j in sequence(n.codons)) {
      if(sum(.codon.sets[i,] == .codon.sets[j,])>=2) { #match at two or three sites of three
        aa1 <- TranslateCodon(paste(n2s(as.numeric(.codon.sets[i,])), collapse=""), numcode=numcode)
        aa2 <- TranslateCodon(paste(n2s(as.numeric(.codon.sets[j,])), collapse=""), numcode=numcode)
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


CreateAAFixationMatrix <- function(aa_op,s,aa.distances,C=4, Phi=0.5, q=4e-7, Ne=5e6){
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


CreateCodonSets <- function() {
  codon.sets <- expand.grid(0:3, 0:3, 0:3)
  codon.sets[,c(3,2,1)] <- codon.sets[,c(1,2,3)] #re-ordering as in the original one
  colnames(.codon.sets) <- c("first", "second", "third")
  codon.sets <- as.matrix(.codon.sets)
  return(.codon.sets)
}


GetLikelihoodSAC_AAForSingleCharGivenOptimum <- function(aa.data, phy, Q_aa, charnum=1, root.p=NULL, return.all=FALSE) {
  #result <- rayDISC(phy=phy, data=aa.data, ntraits=1, charnum=charnum, p=Q_aa, root.p=root.p, node.states="marginal")
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
  #result <- -dev.raydisc(p=NULL, phy=phy, liks=liks, Q=Q_aa, rate=NULL, root.p=root.p)
  result <- -Inf
  ifelse(return.all, stop("return all not currently implemented"), return(result))
}


GetLikelihoodSAC_CodonForSingleCharGivenOptimumHMMScoring <- function(charnum=1, codon.data, phy, Q_codon_array, expQt.codon, root.p=NULL, scale.factor, anc.indices, return.all=FALSE) {
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode

  nl <- 64
  # #Now we need to build the matrix of likelihoods to store node and tip state:
  # if(all(codon.data[,charnum+1] < 65)){
  #   #no need to subset
  #   liks <- Matrix::sparseMatrix(i=1:nb.tip,j=codon.data[,charnum+1],x=1,dims=c(nb.tip + nb.node, nl))
  # } else {
  #   key<-codon.data[,charnum+1] < 65
  #   liks <- Matrix::sparseMatrix(i=which(key),j=codon.data[which(key),charnum+1],x=1,dims=c(nb.tip + nb.node, nl))
  #   liks[which(!key),-c(49, 51, 57)] <- 1
  # }
  # ## Now HMM this matrix by pasting these together:
  # liks.HMM <- cbind(liks,liks,liks,liks,
  #                   liks,liks,liks,liks,
  #                   liks,liks,liks,liks,
  #                   liks,liks,liks,liks,
  #                   Matrix::Matrix(0, nb.tip + nb.node, nl),
  #                   liks,liks,liks,liks)
  liks.HMM <- Matrix::Matrix(0, nb.node, 21* nl)
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  
  if(any(root.p < 0) | any(is.na(root.p))){
    return(1000000)
  }
  #Obtain an object of all the unique ancestors
  for (focal in anc.indices) {
    #the ancestral node at row i is called focal
    #focal <- anc[i]
    #Get descendant information of focal
    desRows<-which(phy$edge[,1]==focal)
    # desNodes<-phy$edge[desRows,2]
    v <- 1
    for (rowIndex in desRows){
      desNode= phy$edge[rowIndex,2]
      if(desNode <= nb.tip){
        if(codon.data[ desNode,charnum+1] < 65){
          v <- v * expQt.codon[[codon.data[ desNode,charnum+1] ]][[desNode]]
          # v <- v * internal_expAtv(A=Q_codon_array ,t=phy$edge.length[rowIndex], v=liks.HMM[phy$edge[rowIndex,2],])
        }
      }else{
        v <- v * internal_expAtv(A=Q_codon_array ,t=phy$edge.length[rowIndex], v=liks.HMM[desNode-nb.tip,])
      }
    }
    comp[focal] <- sum(v)
    liks.HMM[focal-nb.tip,] <- v/comp[focal]
  }
  #Specifies the root:
  # root <- nb.tip + 1L
  #If any of the logs have NAs restart search:
  if(is.nan(sum(log(comp[-TIPS]))) || is.na(sum(log(comp[-TIPS])))){
    return(1000000)
  }else{
    loglik<- (sum(log(comp[-TIPS])) + log(sum(root.p * liks.HMM[1L,])))
    if(is.infinite(loglik)){return(1000000)}
  }
  return(loglik)
  
  
  #The result here is just the likelihood:
  # result <- -TreeTraversalODE(phy=phy, Q_codon_array_vectored=Q_codon_array_vectored, liks.HMM=liks.HMM, bad.likelihood=-100000, root.p=root.p)
  # ifelse(return.all, stop("return all not currently implemented"), return(result))
}


GetLikelihoodSAC_CodonForSingleCharGivenOptimum <- function(charnum=1, codon.data, phy, Q_codon, root.p=NULL, scale.factor, anc.indices, return.all=FALSE) {
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode

  nl <- nrow(Q_codon[[1]])
  #Now we need to build the matrix of likelihoods to pass to dev.raydisc:
  liks <- matrix(0, nb.tip + nb.node, nl)
  if(all(codon.data[,charnum+1] < 65)){
    #no need to subset
    liks[cbind(1:nb.tip,codon.data[,charnum+1])] <- 1
  } else {
    key<-codon.data[,charnum+1] < 65
    liks[cbind(which(key),codon.data[which(key),charnum+1])] <- 1
    liks[which(!key),] <- 1
    if(nl > 4){
      liks[which(!key),c(49, 51, 57)] <- 0
    }
  }
  
  #The result here is just the likelihood:
  result <- -FinishLikelihoodCalculation(phy=phy, liks=liks, Q=Q_codon, root.p=root.p, anc=anc.indices)
  if(return.all) stop("return all not currently implemented");
  return(result)
}


GetLikelihoodSAC_CodonForManyCharGivenFixedOptimumAndQAndRoot <- function(codon.data, phy, Q_codon, root.p=NULL, return.all=FALSE) {
  return(sum(sapply(seq(from=1, to=dim(codon.data)[2]-1, by=1), GetLikelihoodSAC_CodonForSingleCharGivenOptimum, codon.data=codon.data, phy=phy, Q_codon=Q_codon, root.p=root.p, return.all=return.all)))
}


GetLikelihoodSAC_CodonForManyCharVaryingBySiteEvolvingAA <- function(codon.data, phy, Q_codon_array, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, aa.optim_array, codon_mutation_matrix, Ne, rates, numcode, diploid, n.cores.by.gene.by.site=1, verbose=FALSE){

  nsites.unique <- dim(codon.data$unique.site.patterns)[2]-1
  final.likelihood.vector <- rep(NA, nsites.unique)

  #We rescale the codon matrix only:
  diag(codon_mutation_matrix) = 0
  diag(codon_mutation_matrix) = -rowSums(codon_mutation_matrix)
  scale.factor <- -sum(diag(codon_mutation_matrix) * codon.freq.by.gene, na.rm=TRUE)
  codon_mutation_matrix_scaled = codon_mutation_matrix * (1/scale.factor)
  #Finish the Q_array codon mutation matrix multiplication here:
  if(diploid == TRUE){
    Q_codon_array = (2 * Ne) * codon_mutation_matrix_scaled * Q_codon_array
  }else{
    Q_codon_array = Ne * codon_mutation_matrix_scaled * Q_codon_array
  }
  diag(Q_codon_array) = 0
  diag(Q_codon_array) = -rowSums(Q_codon_array)
  Q_codon_array <- Matrix::Matrix(Q_codon_array)
  #Put the na.rm=TRUE bit here just in case -- when the amino acid is a stop codon, there is a bunch of NaNs. Should be fixed now.
  #scale.factor <- -sum(Q_codon_array[DiagArray(dim(Q_codon_array))] * equilibrium.codon.freq, na.rm=TRUE)
  
  #Generate matrix of root frequencies for each optimal AA:
  root.p_array <- codon.freq.by.gene
  #root.p_array <- t(root.p_array)
  #root.p_array <- root.p_array / rowSums(root.p_array)
  #rownames(root.p_array) <- .unique.aa
  phy.sort <- reorder(phy, "pruningwise")
  # Q_codon_array_vectored <- c(t(Q_codon_array)) # has to be transposed
  # Q_codon_array_vectored <- Q_codon_array_vectored[.non_zero_pos]
  
  ## This evaluates exp(Qt)*liks for all of the tip edges for all 64 codon arrangements.
  ## Evaluating multiple tips at once is cheaper than separate evaluation.
  ## In the future, perhaps reduce the number of evaluations for unused codons
  nb.tip <- length(phy.sort$tip.label)
  edge.length.tip <- numeric(nb.tip)
  tip.edge=which(phy.sort$edge[,2]<=nb.tip)
  edge.length.tip[phy.sort$edge[tip.edge,2]] <- phy.sort$edge.length[tip.edge]
  
  #Future: find all tips which contain a codon
  tip.codons <- lapply(seq_len(nb.tip), function(j) unique(codon.data$unique.site.patterns[j,-1]) )
  codon.tips <- lapply(1:64,function(i) sapply(tip.codons,function(j) i%in% j ) )  
  if(!all(sapply(codon.tips,any))) {
    TipLikelihood_codon <- function(codon_number){
      # relevent lines:
      # rowIndex <- which(phy$edge[,1]==focal)
      # v <- v * internal_expAtv(A=Q_codon_array ,t=phy$edge.length[rowIndex], v=liks.HMM[phy$edge[rowIndex,2],])
      if(!any(codon.tips[[codon_number]])) return(list())
      exp_A_tvec_codon(A = Q_codon_array,codon = codon_number,tvec = edge.length.tip)
    }
  }else {
  TipLikelihood_codon <- function(codon_number){
    # relevent lines:
    # rowIndex <- which(phy$edge[,1]==focal)
    # v <- v * internal_expAtv(A=Q_codon_array ,t=phy$edge.length[rowIndex], v=liks.HMM[phy$edge[rowIndex,2],])
    exp_A_tvec_codon(A = Q_codon_array,codon = codon_number,tvec = edge.length.tip)
  }
  }
  # set this up so that:
  # internal_expAtv(A=Q_codon_array ,t=phy$edge.length[rowIndex], v=liks.HMM[phy$edge[rowIndex,2],])
  # can be replaced by:
  # expQt.codon[[ codon.data[ phy$edge[rowIndex,2], charnum+1 ] ]][[ phy$edge[rowIndex,2] ]]
  expQt.codon <- list(list())
  if(n.cores.by.gene.by.site == 1    ){
    expQt.codon <- lapply(1:64, TipLikelihood_codon)
  } else if(64 > n.cores.by.gene.by.site &&  sum(sapply(codon.tips,any)) %/% n.cores.by.gene.by.site < 10 ) {
    expQt.codon <- mclapply(1:64, TipLikelihood_codon, mc.cores=n.cores.by.gene.by.site, mc.preschedule=FALSE)
  } else {
    expQt.codon <- mclapply(1:64, TipLikelihood_codon, mc.cores=n.cores.by.gene.by.site, mc.preschedule=T)
  }
  anc.indices <- unique(phy.sort$edge[,1])
  if(verbose){ 
    MultiCoreLikelihoodBySite <- function(nsite.index){
      tmp <- GetLikelihoodSAC_CodonForSingleCharGivenOptimumHMMScoring(charnum=nsite.index, codon.data=codon.data$unique.site.patterns, 
                                                                       phy=phy.sort, Q_codon_array=Q_codon_array, 
                                                                       expQt.codon=expQt.codon,
                                                                       root.p=root.p_array, scale.factor=scale.factor, 
                                                                       anc.indices=anc.indices, return.all=FALSE)
      cat(".")
      return(tmp)
    }
    
  } else {
  MultiCoreLikelihoodBySite <- function(nsite.index){
    tmp <- GetLikelihoodSAC_CodonForSingleCharGivenOptimumHMMScoring(charnum=nsite.index, codon.data=codon.data$unique.site.patterns, phy=phy.sort, 
                                                                     Q_codon_array=Q_codon_array, 
                                                                     expQt.codon=expQt.codon,
                                                                     root.p=root.p_array, scale.factor=scale.factor, 
                                                                     anc.indices=anc.indices, return.all=FALSE)
    return(tmp)
  }
  }
  if(n.cores.by.gene.by.site == 1    ){
    final.likelihood.vector <- unlist(lapply(1:nsites.unique, MultiCoreLikelihoodBySite))
  } else if(nsites.unique > n.cores.by.gene.by.site &&  nsites.unique %/% n.cores.by.gene.by.site < 10 ) {
    final.likelihood.vector <- unlist(mclapply(1:nsites.unique, MultiCoreLikelihoodBySite, mc.cores=n.cores.by.gene.by.site, mc.preschedule=FALSE))
  } else {
    final.likelihood.vector <- unlist(mclapply(1:nsites.unique, MultiCoreLikelihoodBySite, mc.cores=n.cores.by.gene.by.site, mc.preschedule=T))
  }
  if(verbose) cat("|\n")
  return(final.likelihood.vector)
}


GetLikelihoodSAC_CodonForManyCharVaryingBySite <- function(codon.data, phy, Q_codon_array, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, aa.optim_array, codon_mutation_matrix, Ne, rates, numcode, diploid, n.cores.by.gene.by.site=1) {

  nsites.unique <- dim(codon.data$unique.site.patterns)[2]-1
  final.likelihood.vector <- rep(NA, nsites.unique)
  #unique.aa <- GetMatrixAANames(numcode)

  #We rescale the codon matrix only:
  diag(codon_mutation_matrix) = 0
  diag(codon_mutation_matrix) = -rowSums(codon_mutation_matrix)
  scale.factor <- -sum(diag(codon_mutation_matrix) * codon.freq.by.gene, na.rm=TRUE)
  codon_mutation_matrix_scaled = codon_mutation_matrix * (1/scale.factor)
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

  ## This is obviously not very elegant, but not sure how else to code it to store this stuff in this way -- WORK IN PROGRESS:
  tempGetAAExpQt <- local({
    p0=phy
    Qca=Q_codon_array
    r0=rates
    function(aa) {
      GetExpQt(phy=p0, Q=Qca[,,aa], scale.factor=NULL, rates=r0)
    } })
  expQt <- NULL
  expQt <- mclapply(c("K", "N", "T", "R", "S", 
                      "I", "M", "Q", "H", "P", 
                      "L", "E",  "D", "A", "G", 
                      "V", "Y", "C", "W", "F"),
                    FUN=tempGetAAExpQt, 
                    mc.cores=n.cores.by.gene.by.site)
  names(expQt) <- c("K", "N", "T", "R", "S", 
                    "I", "M", "Q", "H", "P", 
                    "L", "E",  "D", "A", "G", 
                    "V", "Y", "C", "W", "F")
  # expQt$K <- GetExpQt(phy=phy, Q=Q_codon_array[,,"K"], scale.factor=NULL, rates=rates)
  # expQt$N <- GetExpQt(phy=phy, Q=Q_codon_array[,,"N"], scale.factor=NULL, rates=rates)
  # expQt$T <- GetExpQt(phy=phy, Q=Q_codon_array[,,"T"], scale.factor=NULL, rates=rates)
  # expQt$R <- GetExpQt(phy=phy, Q=Q_codon_array[,,"R"], scale.factor=NULL, rates=rates)
  # expQt$S <- GetExpQt(phy=phy, Q=Q_codon_array[,,"S"], scale.factor=NULL, rates=rates)
  # expQt$I <- GetExpQt(phy=phy, Q=Q_codon_array[,,"I"], scale.factor=NULL, rates=rates)
  # expQt$M <- GetExpQt(phy=phy, Q=Q_codon_array[,,"M"], scale.factor=NULL, rates=rates)
  # expQt$Q <- GetExpQt(phy=phy, Q=Q_codon_array[,,"Q"], scale.factor=NULL, rates=rates)
  # expQt$H <- GetExpQt(phy=phy, Q=Q_codon_array[,,"H"], scale.factor=NULL, rates=rates)
  # expQt$P <- GetExpQt(phy=phy, Q=Q_codon_array[,,"P"], scale.factor=NULL, rates=rates)
  # expQt$L <- GetExpQt(phy=phy, Q=Q_codon_array[,,"L"], scale.factor=NULL, rates=rates)
  # expQt$E <- GetExpQt(phy=phy, Q=Q_codon_array[,,"E"], scale.factor=NULL, rates=rates)
  # expQt$D <- GetExpQt(phy=phy, Q=Q_codon_array[,,"D"], scale.factor=NULL, rates=rates)
  # expQt$A <- GetExpQt(phy=phy, Q=Q_codon_array[,,"A"], scale.factor=NULL, rates=rates)
  # expQt$G <- GetExpQt(phy=phy, Q=Q_codon_array[,,"G"], scale.factor=NULL, rates=rates)
  # expQt$V <- GetExpQt(phy=phy, Q=Q_codon_array[,,"V"], scale.factor=NULL, rates=rates)
  # expQt$Y <- GetExpQt(phy=phy, Q=Q_codon_array[,,"Y"], scale.factor=NULL, rates=rates)
  # expQt$C <- GetExpQt(phy=phy, Q=Q_codon_array[,,"C"], scale.factor=NULL, rates=rates)
  # expQt$W <- GetExpQt(phy=phy, Q=Q_codon_array[,,"W"], scale.factor=NULL, rates=rates)
  # expQt$F <- GetExpQt(phy=phy, Q=Q_codon_array[,,"F"], scale.factor=NULL, rates=rates)

  #Generate matrix of root frequencies for each optimal AA:
  root.p_array <- matrix(codon.freq.by.aa, nrow=dim(Q_codon_array)[2], ncol=21)
  root.p_array <- t(root.p_array)
  root.p_array <- root.p_array / rowSums(root.p_array)
  rownames(root.p_array) <- .unique.aa

  phy.sort <- reorder(phy, "pruningwise")
  anc.indices <- unique(phy.sort$edge[,1])
  MultiCoreLikelihoodBySite <- function(i){
    tmp <- GetLikelihoodSAC_CodonForSingleCharGivenOptimum(charnum=i, codon.data=codon.data$unique.site.patterns, phy=phy.sort, Q_codon=expQt[[aa.optim_array[i]]], root.p=root.p_array[aa.optim_array[i],], scale.factor=scale.factor, anc.indices=anc.indices, return.all=FALSE)
    return(tmp)
  }
  final.likelihood.vector.mc <- unlist(mclapply(1:nsites.unique, MultiCoreLikelihoodBySite, mc.cores=n.cores.by.gene.by.site))
  #    final.likelihood.vector <- rep(NA, nsites)
  #for (i in sequence(nsites)) {
  #    final.likelihood.vector[i] <- GetLikelihoodSAC_CodonForSingleCharGivenOptimum(charnum=i, codon.data=codon.data$unique.site.patterns, phy=phy.sort, Q_codon=expQt[[aa.optim_array[i]]], root.p=root.p_array[aa.optim_array[i],], scale.factor=scale.factor, anc.indices=anc.indices, return.all=FALSE)
  #}
  #print(cbind(final.likelihood.vector, final.likelihood.vector.mc))
  return(final.likelihood.vector.mc)
}


GetLikelihoodMutSel_CodonForManyCharVaryingBySite <- function(codon.data, phy, root.p_array=NULL, Q_codon, numcode, n.cores.by.gene.by.site=1) {
  nsites.unique <- dim(codon.data$unique.site.patterns)[2] - 1
  final.likelihood.vector <- rep(NA, nsites.unique)

  diag(Q_codon) = 0
  diag(Q_codon) = -rowSums(Q_codon)
  scale.factor <- -sum(diag(Q_codon) * root.p_array, na.rm=TRUE)
  expQt <- GetExpQt(phy=phy, Q=Q_codon, scale.factor=scale.factor, rates=NULL)

  phy.sort <- reorder(phy, "pruningwise")
  anc.indices <- unique(phy.sort$edge[,1])

  MultiCoreLikelihoodBySite <- function(nsite.index){
    tmp <- GetLikelihoodSAC_CodonForSingleCharGivenOptimum(charnum=nsite.index, codon.data=codon.data$unique.site.patterns, phy=phy.sort, Q_codon=expQt, root.p=root.p_array, scale.factor=scale.factor, anc.indices=anc.indices, return.all=FALSE)
    return(tmp)
  }
  final.likelihood.vector <- unlist(mclapply(1:nsites.unique, MultiCoreLikelihoodBySite, mc.cores=n.cores.by.gene.by.site))
  return(final.likelihood.vector)
}


GetLikelihoodNucleotideForManyCharVaryingBySite <- function(nuc.data, phy, nuc.mutation.rates, include.gamma=FALSE, rates.k=NULL, ncats=NULL, root.p_array=NULL, n.cores.by.gene.by.site=1) {
  nsites.unique <- dim(nuc.data$unique.site.patterns)[2]-1
  final.likelihood.vector <- rep(NA, nsites.unique)
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
  phy.sort <- reorder(phy, "pruningwise")
  anc.indices <- unique(phy.sort$edge[,1])

  MultiCoreLikelihoodBySite <- function(nsite.index){
    tmp <- GetLikelihoodSAC_CodonForSingleCharGivenOptimum(charnum=nsite.index, codon.data=nuc.data$unique.site.patterns, phy=phy, Q_codon=expQt, root.p=root.p_array, scale.factor=scale.factor, anc.indices=anc.indices, return.all=FALSE)
    return(tmp)
  }
  final.likelihood.vector <- unlist(mclapply(1:nsites.unique, MultiCoreLikelihoodBySite, mc.cores=n.cores.by.gene.by.site))
  return(final.likelihood.vector)
}


GetLikelihoodSAC_CodonForManyCharGivenAllParamsEvolvingAA <- function(x, codon.data, phy, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix, include.gamma, gamma.type, ncats, k.levels=0, logspace=FALSE, verbose=TRUE, neglnl=FALSE, n.cores.by.gene.by.site=1, estimate.importance=FALSE) {

  if(logspace) {
    x = exp(x)
  }

  if(estimate.importance == TRUE){
    importance.of.aa.dist.in.selective.environment.change = x[length(x)]
    x = x[-length(x)]
  }else{
    importance.of.aa.dist.in.selective.environment.change = 1
  }

  rate.for.selective.environment.change = x[length(x)]
  x = x[-length(x)]

  if(include.gamma == TRUE){
    shape = x[length(x)]
    x = x[-length(x)]
  }

  C.Phi.q.Ne <- x[1]
  C <- 4
  q <- 4e-7
  Ne <- 5e6
  Phi.q.Ne <- C.Phi.q.Ne / C
  Phi.Ne <- Phi.q.Ne / q
  Phi <- Phi.Ne / Ne
  alpha <- x[2]
  beta <- x[3]
  gamma <- volume.fixed.value

  if(k.levels > 0){
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[9:length(x)], model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[6:length(x)], model=nuc.model, base.freqs=NULL)
      poly.params <- x[4:5]
    }
  }else{
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[7:length(x)], model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[4:length(x)], model=nuc.model, base.freqs=NULL)
    }
  }

  nuc.mutation.rates.vector <- c(nuc.mutation.rates, rate.for.selective.environment.change)
  codon_mutation_matrix <- matrix(nuc.mutation.rates.vector[codon.index.matrix], dim(codon.index.matrix))
  codon_mutation_matrix[is.na(codon_mutation_matrix)]=0
  nsites.unique <- dim(codon.data$unique.site.patterns)[2]-1
  nsites <- sum(codon.data$site.pattern.counts)

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
    #ttmmpp <- c(nuc.mutation.rates.vector, nsites.unique, nsites, C, Phi, rates.k, q, Ne, shape, importance.of.aa.dist.in.selective.environment.change)
    #writeLines(text = paste(ttmmpp), con = "~/Desktop/selac_parameter.txt", sep = "\t")
    final.likelihood.mat = matrix(0, nrow=ncats, ncol=nsites.unique)
    for(k.cat in sequence(ncats)){
      if(k.levels > 0){
        aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
      }else{
        aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
      }
      Q_codon_array <- FastCreateEvolveAACodonFixationProbabilityMatrix(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[k.cat], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999, importance.of.aa.dist.in.selective.environment.change) #Cedric: added importance
      final.likelihood.mat[k.cat,] = GetLikelihoodSAC_CodonForManyCharVaryingBySiteEvolvingAA(codon.data, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site, verbose=verbose)
    }
    likelihood <- sum(log(colSums(exp(final.likelihood.mat)*weights.k)) * codon.data$site.pattern.counts)
  }else{
    if(k.levels > 0){
      aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
    }else{
      aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
    }
    Q_codon_array <- FastCreateEvolveAACodonFixationProbabilityMatrix(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999, importance.of.aa.dist.in.selective.environment.change) #Cedric: added importance
    final.likelihood = GetLikelihoodSAC_CodonForManyCharVaryingBySiteEvolvingAA(codon.data, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site, verbose=verbose)
    likelihood <- sum(final.likelihood * codon.data$site.pattern.counts)
  }

  if(neglnl) {
    likelihood <- -1 * likelihood
  }
  if(verbose > 1) {
    results.vector <- c(likelihood, C*Phi*q, alpha, beta, gamma, Ne, ape::write.tree(phy))
    names(results.vector) <- c("likelihood", "C.Phi.q.Ne", "alpha", "beta", "gamma", "Ne", "phy")
    print(results.vector)
  }else if(verbose){
    results.vector <- c(likelihood, alpha, beta, gamma)
    names(results.vector) <- c("likelihood", "alpha", "beta", "gamma")
    print(results.vector)
    
  }
  if(is.na(likelihood) || is.nan(likelihood)){
    return(1000000)
  }else{
    return(likelihood)
  }
}


GetLikelihoodSAC_CodonForManyCharGivenAllParams <- function(x, codon.data, phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix, include.gamma, gamma.type, ncats, k.levels=0, logspace=FALSE, verbose=TRUE, neglnl=FALSE, n.cores.by.gene.by.site=1) {
  if(logspace) {
    x = exp(x)
  }
  if(include.gamma == TRUE){
    shape = x[length(x)]
    x = x[-length(x)]
  }

  C.Phi.q.Ne <- x[1]
  C <- 4
  q <- 4e-7
  Ne <- 5e6
  Phi.q.Ne <- C.Phi.q.Ne / C
  Phi.Ne <- Phi.q.Ne / q
  Phi <- Phi.Ne / Ne
  alpha <- x[2]
  beta <- x[3]
  gamma <- volume.fixed.value

  if(k.levels > 0){
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[9:length(x)], model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[6:length(x)], model=nuc.model, base.freqs=NULL)
      poly.params <- x[4:5]
    }
  }else{
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[7:length(x)], model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[4:length(x)], model=nuc.model, base.freqs=NULL)
    }
  }

  #codon_mutation_matrix = CreateCodonMutationMatrix(nuc.mutation.rates) #We now make an index matrix first then just place the nucleotide rates into it:
  #codon_mutation_matrix = c(as.vector(nuc.mutation.rates), 0)[codon.index.matrix]
  codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
  codon_mutation_matrix[is.na(codon_mutation_matrix)]=0
  nsites.unique <- dim(codon.data$unique.site.patterns)[2]-1
  nsites <- sum(codon.data$site.pattern.counts)

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
    final.likelihood.mat = matrix(0, nrow=ncats, ncol=nsites.unique)
    for(k.cat in sequence(ncats)){
      if(k.levels > 0){
        aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
      }else{
        aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
      }
      Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[k.cat], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
      final.likelihood.mat[k.cat,] = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
    }
    likelihood <- sum(log(colSums(exp(final.likelihood.mat)*weights.k)) * codon.data$site.pattern.counts)
  }else{
    if(k.levels > 0){
      aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
    }else{
      aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
    }
    Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
    final.likelihood = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
    likelihood <- sum(final.likelihood * codon.data$site.pattern.counts)
  }

  if(neglnl) {
    likelihood <- -1 * likelihood
  }
  if(verbose > 1) {
    results.vector <- c(likelihood, C*Phi*q, alpha, beta, gamma, Ne, ape::write.tree(phy))
    names(results.vector) <- c("likelihood", "C.Phi.q.Ne", "alpha", "beta", "gamma", "Ne", "phy")
    print(results.vector)
  }else if(verbose){
    results.vector <- c(likelihood, alpha, beta, gamma)
    names(results.vector) <- c("likelihood", "alpha", "beta", "gamma")
    print(results.vector)
    
  }
  if(is.na(likelihood) || is.nan(likelihood)){
    return(1000000)
  }else{
    return(likelihood)
  }
}


GetLikelihoodMutSel_CodonForManyCharGivenAllParams <- function(x, codon.data, phy, root.p_array=NULL, numcode, nuc.model, logspace=FALSE, verbose=TRUE, neglnl=FALSE, n.cores.by.gene.by.site=1) {
  if(logspace) {
    x = exp(x)
  }
  if(nuc.model == "JC") {
    base.freqs <- c(x[2:4], 1-sum(x[2:4]))
    nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
    x = x[-c(2:4)]
  }
  if(nuc.model == "GTR") {
    base.freqs <- c(x[2:4], 1-sum(x[2:4]))
    nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[5:9], model=nuc.model, base.freqs=base.freqs)
    x = x[-c(2:9)]
  }
  if(nuc.model == "UNREST") {
    nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[2:12], model=nuc.model, base.freqs=NULL)
    x = x[-c(2:12)]
  }

  #During the early stages of the optimization process it will try weird values for the base frequencies.
  if(any(base.freqs < 0)){
    return(1000000)
  }
  if(!is.null(root.p_array[1])){
    codon.eq.freq <- root.p_array
  }else{
    #.codon.sets <- .codon.sets
    n.codons <- dim(.codon.sets)[1]
    codon.eq.freq <- numeric(n.codons)
    fitness.pars <- c(x[-1],0)
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
  }

  Q_codon = CreateCodonMutationMatrixMutSel(omega.par=x[1], fitness.pars=fitness.pars.ordered, nuc.mutation.rates=nuc.mutation.rates, numcode=numcode)
  final.likelihood <- GetLikelihoodMutSel_CodonForManyCharVaryingBySite(codon.data, phy, root.p_array=codon.eq.freq, Q_codon=Q_codon, numcode=numcode, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
  likelihood <- sum(final.likelihood * codon.data$site.pattern.counts)

  if(neglnl) {
    likelihood <- -1 * likelihood
  }
  if(verbose > 1) {
    results.vector <- c(likelihood, x, ape::write.tree(phy))
    names(results.vector) <- c("likelihood", paste0("param", sequence(length(x))), "phy")
    print(results.vector)
  } else if(verbose){
    results.vector <- c(likelihood, x)
    names(results.vector) <- c("likelihood", paste0("param", sequence(length(x))))
    print(results.vector)
    
  }
  return(likelihood)
}


GetLikelihoodGY94_YN98_CodonForManyCharGivenAllParams <- function(x, codon.data, phy, root.p_array=NULL, model.type="GY94", numcode, logspace=FALSE, verbose=TRUE, neglnl=FALSE, n.cores.by.gene.by.site=1) {
    if(logspace) {
        x = exp(x)
    }
    if(!is.null(root.p_array)){
        codon.freqs <- root.p_array
    }else{
        codon.freqs.tabled <- table(as.matrix(codon.data$unique.site.patterns[,2:dim(codon.data$unique.site.patterns)[2]]))
        codon.freqs <- numeric(64)
        for(codon.index in 1:length(codon.freqs.tabled)){
            codon.freqs[as.numeric(names(codon.freqs.tabled))[codon.index]] <- codon.freqs.tabled[codon.index]
        }
        codon.freqs <- codon.freqs[1:64]/sum(codon.freqs[1:64])
    }
    if(model.type == "GY94"){
        aa.distances <- CreateAADistanceMatrix()
        Q_codon = CreateCodonMutationMatrixGY94(x=x, aa.distances=aa.distances, codon.freqs=codon.freqs, numcode=numcode)
    }else{
        Q_codon = CreateCodonMutationMatrixYN98(x=x, codon.freqs=codon.freqs, numcode=numcode)
    }
    final.likelihood <- GetLikelihoodMutSel_CodonForManyCharVaryingBySite(codon.data, phy, root.p_array=codon.freqs, Q_codon=Q_codon, numcode=numcode, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
    likelihood <- sum(final.likelihood * codon.data$site.pattern.counts)

    if(neglnl) {
        likelihood <- -1 * likelihood
    }
    if(verbose > 1) {
      results.vector <- c(likelihood, x, ape::write.tree(phy))
      names(results.vector) <- c("likelihood", paste0("param", sequence(length(x))), "phy")
      print(results.vector)
    }else if(verbose){
      results.vector <- c(likelihood, x)
      names(results.vector) <- c("likelihood", paste0("param",seq_along(x)))
      print(results.vector)
      
    }
    return(likelihood)
}


GetLikelihoodNucleotideForManyCharGivenAllParams <- function(x, nuc.data, phy, root.p_array=NULL, numcode=1, nuc.model, include.gamma=FALSE, rates.k=NULL, gamma.type="quadrature", ncats=NULL, logspace=FALSE, verbose=TRUE, neglnl=FALSE, n.cores.by.gene.by.site=1) {
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
  nsites.unique <- dim(nuc.data$unique.site.patterns)[2]-1
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
    if(gamma.type == "lognormal"){
      rates.and.weights <- LogNormalQuad(shape=shape, ncats=ncats)
      rates.k <- rates.and.weights[1:ncats]
      weights.k <- rates.and.weights[(ncats+1):(ncats*2)]
    }
    final.likelihood.mat = matrix(0, nrow=ncats, ncol=nsites.unique)
    for(k in sequence(ncats)){
      final.likelihood.mat[k,] = GetLikelihoodNucleotideForManyCharVaryingBySite(nuc.data=nuc.data, phy=phy, nuc.mutation.rates=nuc.mutation.rates, rates.k=rates.k[k], root.p_array=root.p_array, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
    }
    likelihood <- sum(log(colSums(exp(final.likelihood.mat)*weights.k)) * nuc.data$site.pattern.counts)
  }else{
    final.likelihood = GetLikelihoodNucleotideForManyCharVaryingBySite(nuc.data=nuc.data, phy=phy, nuc.mutation.rates=nuc.mutation.rates, rates.k=NULL, root.p_array=root.p_array, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
    likelihood <- sum(final.likelihood * nuc.data$site.pattern.counts)
  }

  if(neglnl) {
    likelihood <- -1 * likelihood
  }
  if(is.na(likelihood)){
    return(1000000)
  }
  if(verbose > 1) {
    results.vector <- c(likelihood, x, ape::write.tree(phy))
    names(results.vector) <- c("likelihood", paste0("param", sequence(length(x))), "phy")
    print(results.vector)
  }else if(verbose){
    results.vector <- c(likelihood, x)
    names(results.vector) <- c("likelihood", paste0("param", seq_along(x)))
    print(results.vector)
    
  }
  return(likelihood)
}


GetOptimalAAPerSite <- function(x, codon.data, phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix, include.gamma=FALSE, gamma.type="quadrature", ncats=4, k.levels=0, logspace=FALSE, verbose=TRUE, neglnl=FALSE, n.cores.by.gene.by.site=1) {
  if(logspace) {
    x = exp(x)
  }
  if(include.gamma == TRUE){
    shape = x[length(x)]
    x = x[-length(x)]
  }

  C.Phi.q.Ne <- x[1]
  C <- 4
  q <- 4e-7
  Ne <- 5e6
  Phi.q.Ne <- C.Phi.q.Ne / C
  Phi.Ne <- Phi.q.Ne / q
  Phi <- Phi.Ne / Ne
  alpha <- x[2]
  beta <- x[3]
  gamma <- volume.fixed.value

  if(k.levels > 0){
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[9:length(x)], model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[6:length(x)], model=nuc.model, base.freqs=NULL)
      poly.params <- x[4:5]
    }
  }else{
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[7:length(x)], model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[4:length(x)], model=nuc.model, base.freqs=NULL)
    }
  }

  if(!is.null(codon.data$unique.site.patterns)){
    codon.data.list <- codon.data
    nsites.unique <- dim(codon.data$unique.site.patterns)[2]-1
  }else{
    nsites.unique <- dim(codon.data)[2]-1
    codon.data.list <- NULL
    codon.data.list$unique.site.patterns <- codon.data
    codon.data.list$site.pattern.counts <- rep(1, nsites.unique)
  }
  nsites <- sum(codon.data.list$site.pattern.counts)

  #codon_mutation_matrix = c(as.vector(nuc.mutation.rates), 0)[codon.index.matrix]
  codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
  codon_mutation_matrix[is.na(codon_mutation_matrix)]=0

  optimal.vector.by.site <- rep(NA, nsites.unique)
  #unique.aa <- GetMatrixAANames(numcode)
  optimal.aa.likelihood.mat <- matrix(0, nrow=length(.unique.aa), ncol=nsites.unique)

  for(i in 1:length(.unique.aa)){
    if(.unique.aa[i]=="*"){
      optimal.aa.likelihood.mat[i,] <- rep(-1000000, nsites.unique)
    }else{
      aa.optim_array = rep(.unique.aa[i], nsites.unique)
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
        final.likelihood.mat = matrix(0, nrow=ncats, ncol=nsites.unique)
        for(k.cat in sequence(ncats)){
          if(k.levels > 0){
            aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
          }else{
            aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
          }
          Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[k.cat], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999)
          tmp = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data.list, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
          tmp[is.na(tmp)] = -1000000
          final.likelihood.mat[k.cat,] = tmp
        }
        optimal.aa.likelihood.mat[i,] <- log(colSums(exp(final.likelihood.mat)*weights.k))
      }else{
        if(k.levels > 0){
          aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
        }else{
          aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
        }
        Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999)
        tmp = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data.list, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
        tmp[is.na(tmp)] = -1000000
        final.likelihood = tmp
        optimal.aa.likelihood.mat[i,] <- final.likelihood
      }
    }
  }
  for(j in 1:nsites.unique){
    optimal.vector.by.site[j] <- .unique.aa[which.is.max(optimal.aa.likelihood.mat[,j])]
  }
  return(optimal.vector.by.site)
}



GetAveAAPerSite <- function(x, codon.data, phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix, include.gamma=FALSE, gamma.type="quadrature", ncats=4, k.levels=0, logspace=FALSE, verbose=TRUE, neglnl=FALSE, n.cores.by.gene.by.site=1) {
  if(logspace) {
    x = exp(x)
  }
  if(include.gamma == TRUE){
    shape = x[length(x)]
    x = x[-length(x)]
  }

  C.Phi.q.Ne <- x[1]
  C <- 4
  q <- 4e-7
  Ne <- 5e6
  Phi.q.Ne <- C.Phi.q.Ne / C
  Phi.Ne <- Phi.q.Ne / q
  Phi <- Phi.Ne / Ne
  alpha <- x[2]
  beta <- x[3]
  gamma <- volume.fixed.value

  if(k.levels > 0){
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[9:length(x)], model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[6:length(x)], model=nuc.model, base.freqs=NULL)
      poly.params <- x[4:5]
    }
  }else{
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[7:length(x)], model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[4:length(x)], model=nuc.model, base.freqs=NULL)
    }
  }

  if(!is.null(codon.data$unique.site.patterns)){
    codon.data.list <- codon.data
    nsites.unique <- dim(codon.data$unique.site.patterns)[2]-1
  }else{
    nsites.unique <- dim(codon.data)[2]-1
    codon.data.list <- NULL
    codon.data.list$unique.site.patterns <- codon.data
    codon.data.list$site.pattern.counts <- rep(1, nsites.unique)
  }
  nsites <- sum(codon.data.list$site.pattern.counts)

  #codon_mutation_matrix = c(as.vector(nuc.mutation.rates), 0)[codon.index.matrix]
  codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
  codon_mutation_matrix[is.na(codon_mutation_matrix)]=0

  likelihood.by.site <- rep(NA, nsites.unique)
  #unique.aa <- GetMatrixAANames(numcode)
  optimal.aa.likelihood.mat <- matrix(0, nrow=length(.unique.aa), ncol=nsites.unique)

  for(i in 1:length(.unique.aa)){
    if(.unique.aa[i]=="*"){
      optimal.aa.likelihood.mat[i,] <- rep(-1000000, nsites.unique)
    }else{
      aa.optim_array = rep(.unique.aa[i], nsites.unique)
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
        final.likelihood.mat = matrix(0, nrow=ncats, ncol=nsites.unique)
        for(k.cat in sequence(ncats)){
          if(k.levels > 0){
            aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
          }else{
            aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
          }
          Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[k.cat], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999)
          tmp = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data.list, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
          tmp[is.na(tmp)] = -1000000
          final.likelihood.mat[k.cat,] = tmp
        }
        optimal.aa.likelihood.mat[i,] <- log(colSums(exp(final.likelihood.mat)*weights.k))
      }else{
        if(k.levels > 0){
          aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
        }else{
          aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
        }
        Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999)
        tmp = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data.list, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
        tmp[is.na(tmp)] = -1000000
        final.likelihood = tmp
        optimal.aa.likelihood.mat[i,] <- final.likelihood
      }
    }
  }
  #Take average
  for(j in 1:nsites.unique){
    #Exclude stop codons, otherwise aritificially inflates likelihood due to our arbitrary setting to -100000. If stop codons are to be included, the following line would have to be modified.
    likelihood.by.site[j] <- log(mean(exp(optimal.aa.likelihood.mat[which(!.unique.aa=="*"),j])))
  }
  final.likelihood <- sum(likelihood.by.site * codon.data$site.pattern.counts)

  if(neglnl) {
    final.likelihood <- -1 * final.likelihood
  }

  return(final.likelihood)
}


#OptimizeModelPars <- function(x, codon.site.data, codon.site.counts, data.type, codon.model, n.partitions, nsites.vector, index.matrix, phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix=NULL, edge.length="optimize", include.gamma=FALSE, gamma.type, ncats, k.levels, logspace=FALSE, verbose=TRUE, parallel.type="by.gene", n.cores=NULL, neglnl=FALSE) {
#    if(logspace) {
#        x <- exp(x)
#    }
#    par.mat <- index.matrix
#    par.mat[] <- c(x, 0)[index.matrix]
#    if(is.null(aa.optim_array)){
#        if(data.type == "nucleotide"){
#            if(nuc.model == "JC"){
#                max.par = 0
#            }
#            if(nuc.model == "GTR"){
#                max.par = 5
#            }
#            if(nuc.model == "UNREST"){
#                max.par = 11
#            }
#            if(include.gamma == TRUE){
#                max.par = max.par + 1
#            }
#            if(is.null(n.cores)){
#                likelihood.vector <- c()
#                for(partition.index in sequence(n.partitions)){
#                    nuc.data = NULL
#                    nuc.data$unique.site.patterns = codon.site.data[[partition.index]]
#                    nuc.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                    likelihood.vector = c(likelihood.vector, GetLikelihoodNucleotideForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, root.p_array=codon.freq.by.gene[[partition.index]], numcode=numcode, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=NULL))
#                }
#                likelihood = sum(likelihood.vector)
#            }else{
#                if(parallel.type=="by.gene"){
#                    MultiCoreLikelihood <- function(partition.index){
#                        nuc.data = NULL
#                        nuc.data$unique.site.patterns = codon.site.data[[partition.index]]
#                        nuc.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                        likelihood.tmp = GetLikelihoodNucleotideForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, root.p_array=codon.freq.by.gene[[partition.index]], numcode=numcode, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=NULL)
#                        return(likelihood.tmp)
#                    }
#                    #This orders the nsites per partition in decreasing order (to increase efficiency):
#                    partition.order <- 1:n.partitions
#                    likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores)))
#                }
#                if(parallel.type == "by.site"){
#                    likelihood.vector <- c()
#                    for(partition.index in sequence(n.partitions)){
#                        nuc.data = NULL
#                        nuc.data$unique.site.patterns = codon.site.data[[partition.index]]
#                        nuc.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                        likelihood.vector = c(likelihood.vector, GetLikelihoodNucleotideForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, root.p_array=codon.freq.by.gene[[partition.index]], numcode=numcode, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=n.cores))
#                    }
#                    likelihood = sum(likelihood.vector)
#                }
#            }
#        }else{
#            if(codon.model == "GY94"){
#                max.par = 2
#                if(is.null(n.cores)){
#                    likelihood.vector <- c()
#                    for(partition.index in sequence(n.partitions)){
#                        codon.data = NULL
#                        codon.data$unique.site.patterns = codon.site.data[[partition.index]]
#                        codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                        likelihood.vector = c(likelihood.vector, GetLikelihoodGY94_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, root.p_array=NULL, numcode=numcode, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=NULL))
#                    }
#                    likelihood = sum(likelihood.vector)
#                }else{
#                    if(parallel.type=="by.gene"){
#                        MultiCoreLikelihood <- function(partition.index){
#                            codon.data = NULL
#                            codon.data$unique.site.patterns = codon.site.data[[partition.index]]
#                            codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                            likelihood.tmp = GetLikelihoodGY94_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, root.p_array=NULL, numcode=numcode, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=NULL)
#                            return(likelihood.tmp)
#                        }
#                        #This orders the nsites per partition in decreasing order (to increase efficiency):
#                        partition.order <- 1:n.partitions
#                        likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores)))
#                    }
#                    if(parallel.type == "by.site"){
#                        likelihood.vector <- c()
#                        for(partition.index in sequence(n.partitions)){
#                            codon.data = NULL
#                            codon.data$unique.site.patterns = codon.site.data[[partition.index]]
#                            codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                            likelihood.vector = c(likelihood.vector, GetLikelihoodGY94_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, root.p_array=NULL, numcode=numcode, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=n.cores))
#                        }
#                        likelihood = sum(likelihood.vector)
#                    }
#                }
#            }else{
#                #To do: figure out way to allow for the crazy 60 fitness par model.
#                if(nuc.model == "JC"){
#                    #base.freq + nuc.rates + omega + fitness.pars
#                    max.par = 3 + 0 + 1 + 19
#                }
#                if(nuc.model == "GTR"){
#                    max.par = 3 + 5 + 1 + 19
#                }
#                if(nuc.model == "UNREST"){
#                    max.par = 0 + 11 + 1 + 19
#                }
#                if(is.null(n.cores)){
#                    likelihood.vector <- c()
#                    for(partition.index in sequence(n.partitions)){
#                        codon.data = NULL
#                        codon.data$unique.site.patterns = codon.site.data[[partition.index]]
#                        codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                        likelihood.vector = c(likelihood.vector, GetLikelihoodMutSel_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, root.p_array=codon.freq.by.gene[[partition.index]], numcode=numcode, nuc.model=nuc.model, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=NULL))
#                    }
#                    likelihood = sum(likelihood.vector)
#                }else{
#                    if(parallel.type=="by.gene"){
#                        MultiCoreLikelihood <- function(partition.index){
#                            codon.data = NULL
#                            codon.data$unique.site.patterns = codon.site.data[[partition.index]]
#                            codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                            likelihood.tmp = GetLikelihoodMutSel_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, root.p_array=codon.freq.by.gene[[partition.index]], numcode=numcode, nuc.model=nuc.model, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=NULL)
#                            return(likelihood.tmp)
#                        }
#                        #This orders the nsites per partition in decreasing order (to increase efficiency):
#                        partition.order <- 1:n.partitions
#                        likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores)))
#                    }
#                    if(parallel.type == "by.site"){
#                        likelihood.vector <- c()
#                        for(partition.index in sequence(n.partitions)){
#                            codon.data = NULL
#                            codon.data$unique.site.patterns = codon.site.data[[partition.index]]
#                            codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                            likelihood.vector = c(likelihood.vector, GetLikelihoodMutSel_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, root.p_array=codon.freq.by.gene[[partition.index]], numcode=numcode, nuc.model=nuc.model, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=n.cores))
#                        }
#                        likelihood = sum(likelihood.vector)
#                    }
#                }
#            }
#        }
#    }else{
#        if(nuc.model == "JC"){
#            max.par = 6
#        }
#        if(nuc.model == "GTR"){
#            max.par = 6 + 5
#        }
#        if(nuc.model == "UNREST"){
#            max.par = 3 + 11
#        }
#        if(include.gamma == TRUE){
#            max.par = max.par + 1
#        }
#        if(k.levels > 0){
#            max.par = max.par + 2
#        }
#        if(is.null(n.cores)){
#            likelihood.vector <- c()
#            for(partition.index in sequence(n.partitions)){
#                codon.data = NULL
#                codon.data$unique.site.patterns = codon.site.data[[partition.index]]
#                codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                likelihood.vector = c(likelihood.vector, GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], codon.freq.by.aa=codon.freq.by.aa[[partition.index]], codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=NULL))
#            }
#            likelihood = sum(likelihood.vector)
#        }else{
#            if(parallel.type == "by.gene"){
#                MultiCoreLikelihood <- function(partition.index){
#                    codon.data = NULL
#                    codon.data$unique.site.patterns = codon.site.data[[partition.index]]
#                    codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                    likelihood.tmp = GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], codon.freq.by.aa=codon.freq.by.aa[[partition.index]], codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=NULL)
#                    return(likelihood.tmp)
#                }
#                #This orders the nsites per partition in decreasing order (to increase efficiency):
#                partition.order <- 1:n.partitions
#                likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores)))
#            }
#            if(parallel.type == "by.site"){
#                likelihood.vector <- c()
#                for(partition.index in sequence(n.partitions)){
#                    codon.data = NULL
#                    codon.data$unique.site.patterns = codon.site.data[[partition.index]]
#                    codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
#                    likelihood.vector = c(likelihood.vector, GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], codon.freq.by.aa=codon.freq.by.aa[[partition.index]], codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=n.cores))
#                }
#                likelihood = sum(likelihood.vector)
#            }
#        }
#    }
#    return(likelihood)
#}


##Redundant to code above. This is a work in progress. Will likely change quite a bit in the future to speed things up.
OptimizeEdgeLengths <- function(x, par.mat, codon.site.data, codon.site.counts, data.type, codon.model, n.partitions, nsites.vector, index.matrix, phy, aa.optim_array=NULL, root.p_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix=NULL, edge.length="optimize", include.gamma=FALSE, gamma.type, ncats, k.levels, logspace=FALSE, verbose=TRUE, n.cores.by.gene, n.cores.by.gene.by.site=1, estimate.importance=FALSE, neglnl=FALSE, HMM=FALSE) {

  if(logspace) {
    x <- exp(x)
  }
  phy$edge.length = x
  if(is.null(aa.optim_array)){
    if(HMM == TRUE) {
      max.par <- dim(par.mat)[2]
      MultiCoreLikelihood <- function(partition.index){
        codon.data = NULL
        codon.data$unique.site.patterns = codon.site.data[[partition.index]]
        codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
        likelihood.tmp = GetLikelihoodSAC_CodonForManyCharGivenAllParamsEvolvingAA(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, codon.freq.by.aa=NULL, codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=estimate.importance)
        return(likelihood.tmp)
      }
      #This orders the nsites per partition in decreasing order (to increase efficiency):
      partition.order <- 1:n.partitions
      likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
    }else{
      if(data.type == "nucleotide"){
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
        MultiCoreLikelihood <- function(partition.index){
          nuc.data = NULL
          nuc.data$unique.site.patterns = codon.site.data[[partition.index]]
          nuc.data$site.pattern.counts = codon.site.counts[[partition.index]]
          likelihood.tmp = GetLikelihoodNucleotideForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, root.p_array=root.p_array[[partition.index]], numcode=numcode, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
          return(likelihood.tmp)
        }
        #This orders the nsites per partition in decreasing order (to increase efficiency):
        partition.order <- 1:n.partitions
        likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
        #There are four options: "none", "GY94", "FMutSel0", "selac".
      }else{
        if(codon.model == "GY94"){
          max.par = 2
          MultiCoreLikelihood <- function(partition.index){
            codon.data = NULL
            codon.data$unique.site.patterns = codon.site.data[[partition.index]]
            codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
            likelihood.tmp = GetLikelihoodGY94_YN98_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, root.p_array=NULL, model.type=codon.model, numcode=numcode, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
            return(likelihood.tmp)
          }
          #This orders the nsites per partition in decreasing order (to increase efficiency):
          partition.order <- 1:n.partitions
          likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
        }
        if(codon.model == "YN98"){
            max.par = 2
            MultiCoreLikelihood <- function(partition.index){
                codon.data = NULL
                codon.data$unique.site.patterns = codon.site.data[[partition.index]]
                codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
                likelihood.tmp = GetLikelihoodGY94_YN98_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, root.p_array=NULL, model.type=codon.model, numcode=numcode, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
                return(likelihood.tmp)
            }
            #This orders the nsites per partition in decreasing order (to increase efficiency):
            partition.order <- 1:n.partitions
            likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
        }
        if(codon.model == "FMutSel0" | codon.model == "FMutSel"){
          if(codon.model == "FMutSel0"){
              if(nuc.model == "JC"){
                  #base.freq + nuc.rates + omega + fitness.pars
                  max.par = 3 + 0 + 1 + 19
              }
              if(nuc.model == "GTR"){
                  max.par = 3 + 5 + 1 + 19
              }
              if(nuc.model == "UNREST"){
                  max.par = 0 + 11 + 1 + 19
              }
          }else{
              if(nuc.model == "JC"){
                  #base.freq + nuc.rates + omega + fitness.pars
                  max.par = 3 + 0 + 1 + 60
              }
              if(nuc.model == "GTR"){
                  max.par = 3 + 5 + 1 + 60
              }
              if(nuc.model == "UNREST"){
                  max.par = 0 + 11 + 1 + 60
              }
          }
          MultiCoreLikelihood <- function(partition.index){
            codon.data = NULL
            codon.data$unique.site.patterns = codon.site.data[[partition.index]]
            codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
            likelihood.tmp = GetLikelihoodMutSel_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, root.p_array=NULL, numcode=numcode, nuc.model=nuc.model, logspace=logspace, verbose=verbose, neglnl=neglnl,  n.cores.by.gene.by.site=n.cores.by.gene.by.site)
            return(likelihood.tmp)
          }
          #This orders the nsites per partition in decreasing order (to increase efficiency):
          partition.order <- 1:n.partitions
          likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
        }
        if(codon.model == "selac"){
          #Triggered the averging of the AA method#
          if(nuc.model == "JC"){
            max.par = 6
          }
          if(nuc.model == "GTR"){
            max.par = 6 + 5
          }
          if(nuc.model == "UNREST"){
            max.par = 3 + 11
          }
          if(include.gamma == TRUE){
            max.par = max.par + 1
          }
          if(k.levels > 0){
            max.par = max.par + 2
          }
          MultiCoreLikelihood <- function(partition.index){
            codon.data = NULL
            codon.data$unique.site.patterns = codon.site.data[[partition.index]]
            codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
            likelihood.tmp = GetAveAAPerSite(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], codon.freq.by.aa=codon.freq.by.aa[[partition.index]], codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
            return(likelihood.tmp)
          }
          #This orders the nsites per partition in decreasing order (to increase efficiency):
          partition.order <- 1:n.partitions
          likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
        }

      }
    }
  }else{
    if(nuc.model == "JC"){
      max.par = 6
    }
    if(nuc.model == "GTR"){
      max.par = 6 + 5
    }
    if(nuc.model == "UNREST"){
      max.par = 3 + 11
    }
    if(include.gamma == TRUE){
      max.par = max.par + 1
    }
    if(k.levels > 0){
      max.par = max.par + 2
    }
    MultiCoreLikelihood <- function(partition.index){
      codon.data = NULL
      codon.data$unique.site.patterns = codon.site.data[[partition.index]]
      codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
      likelihood.tmp = GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], codon.freq.by.aa=codon.freq.by.aa[[partition.index]], codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
      return(likelihood.tmp)
    }
    #This orders the nsites per partition in decreasing order (to increase efficiency):
    partition.order <- 1:n.partitions
    likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
  }
  return(likelihood)
}


#OptimizeModelParsAlphaBetaFixed <- function(x, alpha.beta, codon.site.data, codon.site.counts, data.type, codon.model, n.partitions, nsites.vector, index.matrix, phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix=NULL, edge.length="optimize", include.gamma=FALSE, gamma.type, ncats, k.levels, logspace=FALSE, verbose=TRUE, n.cores=NULL, neglnl=FALSE) {
#    if(logspace) {
#        x <- exp(x)
#    }

#    par.mat.tmp <- index.matrix
#    par.mat.tmp[] <- c(x, 0)[index.matrix]
#if(include.gamma == TRUE){
#par.mat <- matrix(c(par.mat.tmp[1], alpha.beta[1], alpha.beta[2], par.mat.tmp[2:length(par.mat.tmp)], alpha.beta[3]),1, length(par.mat.tmp)+3)
#}else{
#    par.mat <- matrix(c(par.mat.tmp[1], alpha.beta[1], alpha.beta[2], par.mat.tmp[2:length(par.mat.tmp)]),1, length(par.mat.tmp)+2)
#}

#    if(nuc.model == "JC"){
#        max.par = 6
#    }
#    if(nuc.model == "GTR"){
#        max.par = 6 + 5
#    }
#    if(nuc.model == "UNREST"){
#        max.par = 3 + 11
#    }
#    if(include.gamma == TRUE){
#        max.par = max.par + 1
#    }
#    if(k.levels > 0){
#        max.par = max.par + 2
#    }
#    codon.data = NULL
#    codon.data$unique.site.patterns = codon.site.data
#    codon.data$site.pattern.counts = codon.site.counts
#    likelihood.vector = sum(GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=NULL))
#    likelihood = sum(likelihood.vector)

#    return(likelihood)
#}


#OptimizeAlphaBetaOnly <- function(x, par.mat, codon.site.data, codon.site.counts, data.type, codon.model, n.partitions, nsites.vector, index.matrix, phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix=NULL, edge.length="optimize", include.gamma=FALSE, gamma.type, ncats, k.levels, logspace=FALSE, verbose=TRUE, parallel.type="by.gene", n.cores=NULL, neglnl=FALSE) {
#    if(logspace) {
#        x <- exp(x)
#    }
#    par.mat.tmp <- par.mat
#if(include.gamma == TRUE){
#   par.mat <- cbind(par.mat.tmp[,1], x[1], x[2], par.mat.tmp[,2:dim(par.mat.tmp)[2]], x[3])
#}else{
#    par.mat <- cbind(par.mat.tmp[,1], x[1], x[2], par.mat.tmp[,2:dim(par.mat.tmp)[2]])
#}
#    if(nuc.model == "JC"){
#        max.par = 6
#    }
#    if(nuc.model == "GTR"){
#        max.par = 6 + 5
#    }
#    if(nuc.model == "UNREST"){
#        max.par = 3 + 11
#    }
#    if(include.gamma == TRUE){
#        max.par = max.par + 1
#    }
#    if(k.levels > 0){
#        max.par = max.par + 2
#    }
#    if(parallel.type == "by.gene"){
#        MultiCoreLikelihood <- function(partition.index){
#            codon.data = NULL
#            codon.data$unique.site.patterns = codon.site.data[[partition.index]]
#            codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
#           likelihood.tmp = GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], codon.freq.by.aa=codon.freq.by.aa[[partition.index]], codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=NULL)
#            return(likelihood.tmp)
#        }
#This orders the nsites per partition in decreasing order (to increase efficiency):
#        partition.order <- 1:n.partitions
#        likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores)))
#    }
#    if(parallel.type == "by.site"){
#        likelihood.vector <- c()
#        for(partition.index in sequence(n.partitions)){
#            codon.data = NULL
#            codon.data$unique.site.patterns = codon.site.data[[partition.index]]
#            codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
#            likelihood.vector = c(likelihood.vector, GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], codon.freq.by.aa=codon.freq.by.aa[[partition.index]], codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, parallel.type=parallel.type, n.cores=n.cores))
#        }
#        likelihood = sum(likelihood.vector)
#    }
#    return(likelihood)
#}


OptimizeModelParsAlphaBetaGtrFixed <- function(x, alpha.beta.gtr, codon.site.data, codon.site.counts, data.type, codon.model, n.partitions, nsites.vector, index.matrix, phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix=NULL, edge.length="optimize", include.gamma=FALSE, gamma.type, ncats, k.levels, logspace=FALSE, verbose=TRUE, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=1, estimate.importance=FALSE, neglnl=FALSE, HMM=FALSE) {
  if(logspace) {
    x <- exp(x)
  }

  if(HMM == TRUE) {
    if(estimate.importance == TRUE){
      max.par <- length(c(x[1], alpha.beta.gtr, x[2], x[3]))
      par.mat <- matrix(c(x[1], alpha.beta.gtr, x[2], x[3]), 1, max.par)
    }else{
      max.par <- length(c(x[1], alpha.beta.gtr, x[2]))
      par.mat <- matrix(c(x[1], alpha.beta.gtr, x[2]), 1, max.par)
    }
    codon.data = NULL
    codon.data$unique.site.patterns = codon.site.data
    codon.data$site.pattern.counts = codon.site.counts
    likelihood.vector = sum(GetLikelihoodSAC_CodonForManyCharGivenAllParamsEvolvingAA(x=log(par.mat), codon.data=codon.data, phy=phy, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=estimate.importance))
    likelihood = sum(likelihood.vector)
  }else{
    if(is.null(aa.optim_array)){
      if(nuc.model == "JC"){
        max.par = 6
      }
      if(nuc.model == "GTR"){
        max.par = 6 + 5
      }
      if(nuc.model == "UNREST"){
        max.par = 3 + 11
      }
      if(include.gamma == TRUE){
        max.par = max.par + 1
      }
      if(k.levels > 0){
        max.par = max.par + 2
      }

      #THIS ASSUMES A SEPARATE GAMMA PER GENE
      #if(include.gamma == TRUE){
      #    par.mat <- matrix(c(x[1], alpha.beta.gtr, x[2]), 1, max.par)
      #}else{
      #    par.mat <- matrix(c(x[1], alpha.beta.gtr), 1, max.par)
      #}

      par.mat <- matrix(c(x[1], alpha.beta.gtr), 1, max.par)

      codon.data = NULL
      codon.data$unique.site.patterns = codon.site.data
      codon.data$site.pattern.counts = codon.site.counts
      likelihood.vector = sum(GetAveAAPerSite(x=log(par.mat), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site))
      likelihood = sum(likelihood.vector)
    }else{
      if(nuc.model == "JC"){
        max.par = 6
      }
      if(nuc.model == "GTR"){
        max.par = 6 + 5
      }
      if(nuc.model == "UNREST"){
        max.par = 3 + 11
      }
      if(include.gamma == TRUE){
        max.par = max.par + 1
      }
      if(k.levels > 0){
        max.par = max.par + 2
      }

      #THIS ASSUMES A SEPARATE GAMMA PER GENE
      #if(include.gamma == TRUE){
      #    par.mat <- matrix(c(x[1], alpha.beta.gtr, x[2]), 1, max.par)
      #}else{
      #    par.mat <- matrix(c(x[1], alpha.beta.gtr), 1, max.par)
      #}

      par.mat <- matrix(c(x[1], alpha.beta.gtr), 1, max.par)

      codon.data = NULL
      codon.data$unique.site.patterns = codon.site.data
      codon.data$site.pattern.counts = codon.site.counts
      likelihood.vector = sum(GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site))
      likelihood = sum(likelihood.vector)
    }
  } ## end else for HMM == TRUE
  return(likelihood)
}


OptimizeAlphaBetaGtrOnly <- function(x, fixed.pars, codon.site.data, codon.site.counts, data.type, codon.model, n.partitions, nsites.vector, index.matrix, phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix=NULL, edge.length="optimize", include.gamma=FALSE, gamma.type=gamma.type, ncats, k.levels, logspace=FALSE, verbose=TRUE, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=1, neglnl=FALSE, HMM=FALSE, estimate.importance=FALSE) {
  if(logspace) {
    x <- exp(x)
  }
  if(HMM == TRUE) {
    if(estimate.importance == TRUE){
      par.mat <- c()
      for(row.index in 1:dim(fixed.pars)[1]){
        par.mat <- rbind(par.mat, c(fixed.pars[row.index,1], x, fixed.pars[row.index,2], fixed.pars[row.index,3]))
      }
    }else{
      par.mat <- c()
      for(row.index in 1:dim(fixed.pars)[1]){
        par.mat <- rbind(par.mat, c(fixed.pars[row.index,1], x, fixed.pars[row.index,2]))
      }
    }
    max.par <- dim(par.mat)[2]
    MultiCoreLikelihood <- function(partition.index){
      codon.data = NULL
      codon.data$unique.site.patterns = codon.site.data[[partition.index]]
      codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
      likelihood.tmp = GetLikelihoodSAC_CodonForManyCharGivenAllParamsEvolvingAA(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, codon.freq.by.aa=NULL, codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=estimate.importance)
      return(likelihood.tmp)
    }

    #This orders the nsites per partition in decreasing order (to increase efficiency):
    partition.order <- 1:n.partitions
    likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))

  }else{
    if(is.null(aa.optim_array)){
      if(nuc.model == "JC"){
        max.par = 6
      }
      if(nuc.model == "GTR"){
        max.par = 6 + 5
      }
      if(nuc.model == "UNREST"){
        max.par = 3 + 11
      }
      if(include.gamma == TRUE){
        max.par = max.par + 1
      }
      if(k.levels > 0){
        max.par = max.par + 2
      }

      #THIS ASSUMES SEPARATE GAMMA PER GENE:
      #    if(include.gamma == TRUE){
      #    par.mat <- c()
      #    for(row.index in 1:dim(fixed.pars)[1]){
      #        par.mat <- rbind(par.mat, c(fixed.pars[row.index,1], x, fixed.pars[row.index,2]))
      #    }
      #}else{
      #    par.mat <- c()
      #    for(row.index in 1:dim(fixed.pars)[1]){
      #        par.mat <- rbind(par.mat, c(fixed.pars[row.index,1], x))
      #    }
      #}

      par.mat <- c()
      for(row.index in 1:dim(fixed.pars)[1]){
        par.mat <- rbind(par.mat, c(fixed.pars[row.index,1], x))
      }

      MultiCoreLikelihood <- function(partition.index){
        codon.data = NULL
        codon.data$unique.site.patterns = codon.site.data[[partition.index]]
        codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
        likelihood.tmp = GetAveAAPerSite(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], codon.freq.by.aa=codon.freq.by.aa[[partition.index]], codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
        return(likelihood.tmp)
      }
      #This orders the nsites per partition in decreasing order (to increase efficiency):
      partition.order <- 1:n.partitions
      likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
    }else{
      if(nuc.model == "JC"){
        max.par = 6
      }
      if(nuc.model == "GTR"){
        max.par = 6 + 5
      }
      if(nuc.model == "UNREST"){
        max.par = 3 + 11
      }
      if(include.gamma == TRUE){
        max.par = max.par + 1
      }
      if(k.levels > 0){
        max.par = max.par + 2
      }

      #THIS ASSUMES SEPARATE GAMMA PER GENE:
      #    if(include.gamma == TRUE){
      #    par.mat <- c()
      #    for(row.index in 1:dim(fixed.pars)[1]){
      #        par.mat <- rbind(par.mat, c(fixed.pars[row.index,1], x, fixed.pars[row.index,2]))
      #    }
      #}else{
      #    par.mat <- c()
      #    for(row.index in 1:dim(fixed.pars)[1]){
      #        par.mat <- rbind(par.mat, c(fixed.pars[row.index,1], x))
      #    }
      #}

      par.mat <- c()
      for(row.index in 1:dim(fixed.pars)[1]){
        par.mat <- rbind(par.mat, c(fixed.pars[row.index,1], x))
      }

      MultiCoreLikelihood <- function(partition.index){
        codon.data = NULL
        codon.data$unique.site.patterns = codon.site.data[[partition.index]]
        codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
        likelihood.tmp = GetLikelihoodSAC_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array[[partition.index]], codon.freq.by.aa=codon.freq.by.aa[[partition.index]], codon.freq.by.gene=codon.freq.by.gene[[partition.index]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
        return(likelihood.tmp)
      }
      #This orders the nsites per partition in decreasing order (to increase efficiency):
      partition.order <- 1:n.partitions
      likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
    }
  }
  return(likelihood)
}


OptimizeModelParsLarge <- function(x, codon.site.data, codon.site.counts, data.type, codon.model, n.partitions, nsites.vector, index.matrix, phy, aa.optim_array=NULL, root.p_array=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix=NULL, edge.length="optimize", include.gamma=FALSE, gamma.type, ncats, k.levels, logspace=FALSE, verbose=TRUE, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=1, neglnl=FALSE) {

    poo <- x
    if(logspace) {
        x <- exp(x)
    }

    if(class(index.matrix)=="numeric"){
        index.matrix <- matrix(index.matrix, 1, length(index.matrix))
    }
    par.mat <- index.matrix
    par.mat[] <- c(x, 0)[index.matrix]

    if(data.type == "nucleotide"){
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
        likelihood.vector <- c()
        for(partition.index in sequence(n.partitions)){
            nuc.data = NULL
            nuc.data$unique.site.patterns = codon.site.data
            nuc.data$site.pattern.counts = codon.site.counts
            likelihood.vector = c(likelihood.vector, GetLikelihoodNucleotideForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), nuc.data=nuc.data, phy=phy, root.p_array=root.p_array, numcode=numcode, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site))
        }
        likelihood = sum(likelihood.vector)
    }else{
        if(codon.model == "GY94" | codon.model == "YN98"){
            max.par = 2
            MultiCoreLikelihood <- function(partition.index){
                codon.data = NULL
                codon.data$unique.site.patterns = codon.site.data[[partition.index]]
                codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
                likelihood.tmp <- c()
                try(likelihood.tmp <- GetLikelihoodGY94_YN98_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, root.p_array=NULL, model.type=codon.model, numcode=numcode, logspace=logspace, verbose=verbose, neglnl=neglnl, n.cores.by.gene.by.site=n.cores.by.gene.by.site))
                if(length(likelihood.tmp)==0){
                    return(10000000)
                }else{
                    return(likelihood.tmp)
                }
            }
            #This orders the nsites per partition in decreasing order (to increase efficiency):
            partition.order <- 1:n.partitions
            likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
        }
        if(codon.model == "FMutSel0" | codon.model == "FMutSel") {
            if(codon.model == "FMutSel0"){
                #To do: figure out way to allow for the crazy 60 fitness par model.
                if(nuc.model == "JC"){
                    #base.freq + nuc.rates + omega + fitness.pars
                    max.par = 3 + 0 + 1 + 19
                }
                if(nuc.model == "GTR"){
                    max.par = 3 + 5 + 1 + 19
                }
                if(nuc.model == "UNREST"){
                    max.par = 0 + 11 + 1 + 19
                }
            }else{
                #To do: figure out way to allow for the crazy 60 fitness par model.
                if(nuc.model == "JC"){
                    #base.freq + nuc.rates + omega + fitness.pars
                    max.par = 3 + 0 + 1 + 60
                }
                if(nuc.model == "GTR"){
                    max.par = 3 + 5 + 1 + 60
                }
                if(nuc.model == "UNREST"){
                    max.par = 0 + 11 + 1 + 60
                }
            }
            save(poo, phy, index.matrix, codon.site.data, codon.site.counts, file="checkpoint.fmutsel.Rsave")
            MultiCoreLikelihood <- function(partition.index){
                codon.data = NULL
                codon.data$unique.site.patterns = codon.site.data[[partition.index]]
                codon.data$site.pattern.counts = codon.site.counts[[partition.index]]
                likelihood.tmp <- c()
                try(likelihood.tmp <- GetLikelihoodMutSel_CodonForManyCharGivenAllParams(x=log(par.mat[partition.index,1:max.par]), codon.data=codon.data, phy=phy, root.p_array=NULL, numcode=numcode, nuc.model=nuc.model, logspace=logspace, verbose=verbose, neglnl=neglnl,  n.cores.by.gene.by.site=n.cores.by.gene.by.site))
                if(length(likelihood.tmp)==0){
                    return(10000000)
                }else{
                    return(likelihood.tmp)
                }
            }
            #This orders the nsites per partition in decreasing order (to increase efficiency):
            partition.order <- 1:n.partitions
            likelihood <- sum(unlist(mclapply(partition.order[order(nsites.vector, decreasing=TRUE)], MultiCoreLikelihood, mc.cores=n.cores.by.gene)))
        }
    }
    return(likelihood)
}


ComputeStartingBranchLengths <- function(phy, data, data.type="codon", recalculate.starting.brlen){
  if(recalculate.starting.brlen || is.null(phy$edge.length)) {
    if(is.null(phy$edge.length)){
      phy$edge.length = rep(1, length(phy$edge[,1]))
    }
    data.mat <- DNAbinToNucleotideCharacter(data)
    new.tip<-list(edge=matrix(c(2L,1L),1,2),tip.label="FAKEY_MCFAKERSON", edge.length=1, Nnode=1L)
    class(new.tip) <- "phylo"
    phy.with.outgroup <- bind.tree(phy, new.tip,where="root")
    new.tip.data <- matrix(c("FAKEY_MCFAKERSON", rep("-", dim(data.mat)[2]-1)), dim(data.mat)[2], 1)
    new.tip.df <- as.data.frame(t(new.tip.data))
    rownames(new.tip.df) <- "FAKEY_MCFAKERSON"
    colnames(new.tip.df) <- colnames(data.mat)
    data.with.outgroup <- rbind(data.mat, new.tip.df)
    dats.mat <- as.matrix(data.with.outgroup[,-1])
    if(data.type=="codon"){
      third.position <- seq(3,dim(dats.mat)[2], by=3)
      dat <- phyDat(dats.mat[,third.position], type="DNA")
      mpr.tre <- acctran(phy.with.outgroup, dat)
      mpr.tre$edge.length <- mpr.tre$edge.length/dim(dats.mat[,third.position])[2]
    }else{
      dat <- phyDat(dats.mat, type="DNA")
      mpr.tre <- acctran(phy.with.outgroup, dat)
      mpr.tre$edge.length <- mpr.tre$edge.length/dim(dats.mat)[2]
    }
    mpr.tre.pruned <- drop.tip(mpr.tre, "FAKEY_MCFAKERSON")
    mpr.tre.pruned$edge.length[mpr.tre.pruned$edge.length == 0] <- 1e-7
  } else {
    mpr.tre.pruned <- phy
  }
  return(mpr.tre.pruned)
}


GetFitnessStartingValues <- function(codon.freqs, n.pars = 21){
  initial.vals <- c()
  for(i in 1:n.pars){
    #This is taken directly from PAML code:
    initial.vals <- c(initial.vals, (codon.freqs[i]+.001)/(codon.freqs[n.pars]+.002*runif(1)))
  }
  if(n.pars == 64){
    initial.vals = initial.vals[-c(49,51,57,64)]
  }
  initial.vals[initial.vals < 0.0001] <- 0.001
  return(initial.vals)
}


GetCAI <- function(codon.data, aa.optim, numcode=1, w){
  ref.codon.freqs <- GetCodonFreqsByAA(codon.data[1,-1], aa.optim, numcode=1)
  w.array <- matrix(ref.codon.freqs, nrow=64, ncol=21)
  w.array <- t(w.array)
  w.array <- w.array / apply(w.array, 1, max)
  w.array[17,] <- 0
  #unique.aa <- GetMatrixAANames(numcode=numcode)
  rownames(w.array) <- .unique.aa
  wi <- apply(w.array, 2, max)
  wi[wi < 1e-4] <- 0.01
  wi = w
  cai <- exp((1/(dim(codon.data)[2]-1)) * sum(log(wi[as.numeric(codon.data[1,])][-1])))
  return(cai)
}


DiscreteGamma <- function (shape, ncats){
  quantiles <- qgamma((1:(ncats - 1))/ncats, shape = shape, rate = shape)
  return(diff(c(0, pgamma(quantiles * shape, shape + 1, rate=shape), 1)) * ncats)
}


LogNormalQuad <- function(shape, ncats){
  s = shape
  m = -(s^2)/2
  pp <- gauss.quad.prob(ncats, dist="normal", mu=m, sigma=s)
  return(c(exp(pp$nodes/m), pp$weights))
}


LaguerreQuad <- function(shape, ncats) {
  # Determine rates based on alpha and the number of bins
  # bins roots normalized to 1 of the General Laguerre Quadrature
  # first ncats elements are rates with mean 1
  # second ncats elements are probabilities with sum 1
  roots <- findRoots(shape - 1, ncats)
  weights <- numeric(ncats)
  f <- prod(1 + (shape - 1)/(1:ncats))

  for (i in 1:ncats) {
    weights[i] <- f*roots[i]/((ncats + 1)^2*Laguerre(roots[i], shape - 1, ncats + 1)^2)
  }
  roots <- roots/shape
  return(c(roots, weights))
}


findRoots <- function(shape, ncats) {
  # Determine rates based on Gamma's alpha and the number of bins
  # bins roots normalized to 1 of the General Laguerre Polynomial (GLP)
  coeff  <- integer(ncats + 1)
  for (i in 0:ncats) {
    coeff[i + 1] <- (-1)^i*nChooseK(ncats + shape, ncats - i)/factorial(i)
  }
  return(sort(Re(polyroot(coeff))))
}


Laguerre <- function(x, shape, degree) {
  y <- 0
  for (i in 0:degree) {
    y <- y + (-1)^i*choose(degree + shape, degree - i)*x^i/factorial(i)
  }
  return(y)
}


#Took this from R.basic -- the C version did not work when LaguerreQuad was called internally. Adding this function fixed this issue (JMB 9-29-2016).
nChooseK <- function(n, k, log=FALSE) {
  nChooseK0 <- function(n, k) {
    if((n == k) || (k==0))
      return(1);
    m <- min(k, n-k);
    prod(seq(from=n, to=(n-m+1), by=-1)/(seq(from=m, to=1, by=-1)));
  }
  # Process the arguments
  if (is.logical(log)) {
    if (log == TRUE)
      log <- exp(1)
    else
      log <- NULL;
  }
  # Repeat n or k to make the of equal length.
  nn <- length(n);
  nk <- length(k);
  if (nn > nk) {
    k <- rep(k, length.out=nn);
    nk <- nn;
  } else if (nn < nk) {
    n <- rep(n, length.out=nk);
    nn <- nk;
  }
  if (is.null(log)) {
    gamma(n+1) / (gamma(n-k+1) * gamma(k+1));
  } else {
    (lgamma(n+1) - (lgamma(n-k+1) + lgamma(k+1))) / log(log);
  }
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
  diag(x) <- 0
  x<-x/max(x)
  g <- igraph::graph.adjacency(x, weighted=TRUE, mode="directed")
  g.layout <- igraph::layout.fruchterman.reingold(g)
  plot(g, layout=g.layout, edge.width=10*igraph::get.edge.attribute(g, "weight"), edge.curved=TRUE)
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


DNAbinToCodonCharacter <- function(x, frame=0, corHMM.format=TRUE) {
  bound.characters <- sapply(as.character(x), paste, collapse="")
  #following fn is derived from code for uco in seqinr
  SplitToCodons <- function(seq.string, frame) {
    seq.string<-strsplit(seq.string, split="")[[1]]
    if (any(seq.string %in% LETTERS)) {
      seq.string <- tolower(seq.string)
    }
    return(sapply(splitseq(seq = seq.string, frame = frame, word = 3), CodonStringToCharacter))
  }
  split.characters <- t(sapply(bound.characters, SplitToCodons, frame=frame))
  colnames(split.characters) <- sequence(dim(split.characters)[2])
  if(corHMM.format) {
    split.characters <- cbind(data.frame(Taxa=rownames(split.characters)), data.frame(split.characters))
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


DNAbinToNucleotideCharacter <- function(x, frame=0, corHMM.format=TRUE) {
  bound.characters <- sapply(as.character(x), paste, collapse="")
  #following fn is derived from code for uco in seqinr
  SplitToCodons <- function(seq.string, frame) {
    seq.string<-strsplit(seq.string, split="")[[1]]
    if (any(seq.string %in% LETTERS)) {
      seq.string <- tolower(seq.string)
    }
    return(sapply(splitseq(seq = seq.string, frame = frame, word = 1), NucleotideStringToCharacter))
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
  #codon.sets <- CreateCodonSets()
  #codon.set.translate <- apply(.codon.sets, 2, n2s)
  #codon.name <- apply(.codon.set.translate, 1, paste, collapse="")
  codon.aa <- sapply(.codon.name,TranslateCodon, numcode=numcode)
  names(codon.aa ) = NULL
  #unique.aa <- unique(codon.aa)
  return(.unique.aa)
}


GetCodonFreqsByAA <- function(codon.data, aa.opt.vector, numcode){
  #codon.sets <- CreateCodonSets()
  #codon.set.translate <- apply(.codon.sets, 2, n2s)
  #codon.name <- apply(.codon.set.translate, 1, paste, collapse="")
  aa.translations <- .aa.translation[[numcode=numcode]][codon.data=codon.data[,1]]
  names(aa.translations) = NULL
  #unique.aa <- unique(aa.translation)
  codon.freqs <- c()
  for(aa.id.index in sequence(21)) {
    cols <- which(aa.opt.vector == .unique.aa[aa.id.index])
    codon.freqs.tmp <- rep(0, 64)
    for(col.index in sequence(length(cols))) {
      for(row.index in sequence(dim(codon.data)[1])) {
        if(codon.data[row.index, cols[col.index]]<65){
          codon.freqs.tmp[codon.data[row.index, cols[col.index]]] <- codon.freqs.tmp[codon.data[row.index, cols[col.index]]] + 1
        }
      }
    }
    codon.freqs <- c(codon.freqs, codon.freqs.tmp)
  }
  return(codon.freqs)
}


GetCodonFreqsByGene <- function(codon.data){
  codon.freqs.tabled <- table(as.matrix(codon.data[,2:dim(codon.data)[2]]))
  codon.freqs <- numeric(64)
  for(codon.index in 1:length(codon.freqs)){
    codon.freqs[as.numeric(names(codon.freqs.tabled))[codon.index]] <- codon.freqs.tabled[codon.index]
  }
  codon.freqs <- codon.freqs[1:64]/sum(codon.freqs[1:64])
  return(codon.freqs)
}


GetCodonFreqsByGeneHMM <- function(codon.data){
  codon.freqs.tabled <- table(as.matrix(codon.data[,2:dim(codon.data)[2]]))
  codon.freqs <- numeric(64)
  for(codon.index in 1:length(codon.freqs)){
    codon.freqs[as.numeric(names(codon.freqs.tabled))[codon.index]] <- codon.freqs.tabled[codon.index]
  }
  codon.freqs <- codon.freqs[1:64]/sum(codon.freqs[1:64])
  codon.freqs.full <- c()
  for(amino.acid.index in 1:21) {
    codon.freqs.full <- c(codon.freqs.full, codon.freqs)
  }
  codon.freqs.full <- codon.freqs.full/sum(codon.freqs.full)
  return(codon.freqs.full)
}


GetAAFreqsByGene <- function(codon.data, aa.opt.vector, numcode){
  #codon.sets <- CreateCodonSets()
  #codon.set.translate <- apply(.codon.sets, 2, n2s)
  #codon.name <- apply(.codon.set.translate, 1, paste, collapse="")
  aa.translations <- .aa.translation[[numcode]][.codon.name]
  names(aa.translations) = NULL
  #unique.aa <- unique(aa.translation)
  eq.freqs <- c()
  for(aa.id.index in sequence(21)) {
    cols <- which(aa.opt.vector == .unique.aa[aa.id.index])
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
  root.p_array <- matrix(eq.freqs, nrow=64, ncol=21)
  root.p_array <- t(root.p_array)
  eq.freqs <- rowSums(root.p_array)
  eq.freqs <- as.vector(eq.freqs / sum(eq.freqs))

  return(eq.freqs)
}


GetMaxName <- function(x) {
  x.tmp <- x
  x.tmp <- x[-which(x == "NA")]
  if(length(x.tmp)==0){
    return(names(table(x))[(which.is.max(table(x)))])
  }else{
    return(names(table(x.tmp))[(which.is.max(table(x.tmp)))])
  }
}



######################################################################################################################################
######################################################################################################################################
### Likelihood calculator -- Two step process
######################################################################################################################################
######################################################################################################################################

#Step 1: We perform exponentiation as few times as possible; in the case of selac, for example, we do it for each of the possible optimal amino acids and store the matrices.
#code <- "
#Eigen::MatrixXd m = Rcpp::as<Eigen::MatrixXd>(a);
#return(Rcpp::wrap(m.exp()));"
#eigenExpM <- cxxfunction(signature(a="numeric"), code, plugin="RcppEigen")

#Use specialized expm, copied from package expm
internal_expm <- function (x) {
  stopifnot(is.numeric(x) || (isM <- inherits(x, "dMatrix")) || 
              inherits(x, "mpfrMatrix"))
  if (length(d <- dim(x)) != 2) 
    stop("argument is not a matrix")
  if (d[1] != d[2]) 
    stop("matrix not square")
  method <- "Higham08"
  preconditioning = "2bal"
  A<-x
  n <- d[1]
  if (n <= 1) 
    return(exp(A))
  # if (balancing) {
  baP <- expm::balance(A, "P")
  baS <- expm::balance(baP$z, "S")
  A <- baS$z
  # }
  nA <- Matrix::norm(A, "1")
  I <- if (is(A, "Matrix")) 
    Matrix::Diagonal(n)
  else diag(n)
  if (nA <= 2.1) {
    t <- c(0.015, 0.25, 0.95, 2.1)
    l <- which.max(nA <= t)
    C <- rbind(c(120, 60, 12, 1, 0, 0, 0, 0, 0, 0), 
               c(30240, 15120, 3360, 420, 30, 1, 0, 0, 0, 0), 
               c(17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1, 0,  0), 
               c(17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1))
    A2 <- A %*% A
    P <- I
    U <- C[l, 2] * I
    V <- C[l, 1] * I
    for (k in 1:l) {
      P <- P %*% A2
      U <- U + C[l, (2 * k) + 2] * P
      V <- V + C[l, (2 * k) + 1] * P
    }
    U <- A %*% U
    X <- solve(V - U, V + U)
  }
  else {
    s <- log2(nA/5.4)
    B <- A
    if (s > 0) {
      s <- ceiling(s)
      B <- B/(2^s)
    }
    c. <- c(64764752532480000, 32382376266240000, 7771770303897600, 
            1187353796428800, 129060195264000, 10559470521600, 
            670442572800, 33522128640, 1323241920, 40840800, 
            960960, 16380, 182, 1)
    B2 <- B %*% B
    B4 <- B2 %*% B2
    B6 <- B2 %*% B4
    U <- B %*% (B6 %*% (c.[14] * B6 + c.[12] * B4 + c.[10] *  B2) + 
                  c.[8] * B6 + c.[6] * B4 + c.[4] * B2 + c.[2] * I)
    V <- B6 %*% (c.[13] * B6 + c.[11] * B4 + c.[9] * B2) + 
      c.[7] * B6 + c.[5] * B4 + c.[3] * B2 + c.[1] * I
    X <- solve(V - U, V + U)
    if (s > 0) 
      for (t in 1:s) X <- X %*% X
  }
  # if (balancing) {
  d <- baS$scale
  X <- X * (d * rep(1/d, each = n))
  pp <- as.integer(baP$scale)
  if (baP$i1 > 1) {
    for (i in (baP$i1 - 1):1) {
      tt <- X[, i]
      X[, i] <- X[, pp[i]]
      X[, pp[i]] <- tt
      tt <- X[i, ]
      X[i, ] <- X[pp[i], ]
      X[pp[i], ] <- tt
    }
  }
  if (baP$i2 < n) {
    for (i in (baP$i2 + 1):n) {
      tt <- X[, i]
      X[, i] <- X[, pp[i]]
      X[, pp[i]] <- tt
      tt <- X[i, ]
      X[i, ] <- X[pp[i], ]
      X[pp[i], ] <- tt
    }
  }
  # }
  return(X)
}
internal_expmt <- function (A, t_vec) {
  stopifnot(is.numeric(A) || (isM <- inherits(A, "dMatrix")) || 
              inherits(A, "mpfrMatrix"))
  if (length(d <- dim(A)) != 2) 
    stop("argument is not a matrix")
  if (d[1] != d[2]) 
    stop("matrix not square")
  method <- "Higham08"
  preconditioning = "2bal"
  n <- d[1]
  if (n <= 1) 
    return(as.list(exp(A*t_vec)))  # force a list format return
  # if (balancing) {
  baP <- expm::balance(A, "P")
  baS <- expm::balance(baP$z, "S")
  A <- baS$z
  # }
  nA <- Matrix::norm(A, "1")
  I <- if (is(A, "Matrix")) 
    Matrix::Diagonal(n)
  else diag(n)
  res <- as.list(numeric(length(t_vec)))
  
  C <- rbind(c(120, 60, 12, 1, 0, 0, 0, 0, 0, 0), 
             c(30240, 15120, 3360, 420, 30, 1, 0, 0, 0, 0), 
             c(17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1, 0,  0), 
             c(17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1))
  
  c. <- c(64764752532480000, 32382376266240000, 7771770303897600, 
          1187353796428800, 129060195264000, 10559470521600, 
          670442572800, 33522128640, 1323241920, 40840800, 
          960960, 16380, 182, 1)
  t <- c(0.015, 0.25, 0.95, 2.1)
  sA <- log2(nA/5.4)
  A2_base <- A %*% A
  if(any(nA * abs(as.numeric(t_vec)) > 2.1)){
    A4_base <- A2_base %*% A2_base
    A6_base <- A2_base %*% A4_base
  }
  for(res_i in seq_len(length(t_vec))){
    t_i=t_vec[res_i]
    if (nA * abs(t_i) <= 2.1) {
      l <- which.max(nA * abs(t_i)  <= t)
      A2 <- A2_base * t_i * t_i
      P <- I
      U <- C[l, 2] * I
      V <- C[l, 1] * I
      for (k in 1:l) {
        P <- P %*% A2
        U <- U + C[l, (2 * k) + 2] * P
        V <- V + C[l, (2 * k) + 1] * P
      }
      U <- A %*% U * t_i
      X <- solve(V - U, V + U)
    }
    else {
      s <- sA + log2(abs(t_i))
      B <- A * t_i
      if (s > 0) {
        s <- ceiling(s)
        B <- B/(2^s)
        B2 <- A2_base * t_i * t_i / (4^s)
        B4 <- A4_base * t_i ^ 4 / (16^s)
        B6 <- A6_base * t_i ^ 6 / (64^s)
      } else {
        B2 <- A2_base * t_i * t_i
        B4 <- A4_base * t_i ^ 4 
        B6 <- A6_base * t_i ^ 6 
      }
      U <- B %*% (B6 %*% (c.[14] * B6 + c.[12] * B4 + c.[10] *  B2) + 
                    c.[8] * B6 + c.[6] * B4 + c.[4] * B2 + c.[2] * I)
      V <- B6 %*% (c.[13] * B6 + c.[11] * B4 + c.[9] * B2) + 
        c.[7] * B6 + c.[5] * B4 + c.[3] * B2 + c.[1] * I
      X <- solve(V - U, V + U)
      if (s > 0) 
        for (t in 1:s) X <- X %*% X
    }
    # if (balancing) {
    d <- baS$scale
    X <- X * (d * rep(1/d, each = n))
    pp <- as.integer(baP$scale)
    if (baP$i1 > 1) {
      for (i in (baP$i1 - 1):1) {
        tt <- X[, i]
        X[, i] <- X[, pp[i]]
        X[, pp[i]] <- tt
        tt <- X[i, ]
        X[i, ] <- X[pp[i], ]
        X[pp[i], ] <- tt
      }
    }
    if (baP$i2 < n) {
      for (i in (baP$i2 + 1):n) {
        tt <- X[, i]
        X[, i] <- X[, pp[i]]
        X[, pp[i]] <- tt
        tt <- X[i, ]
        X[i, ] <- X[pp[i], ]
        X[pp[i], ] <- tt
      }
    }
    # }
    res[[res_i]]<-X
  }
  return(res)
}

## This is a editted copy from the expm package; will be replaced with HMM variants for fixed values
internal_expAtv <- function(A, v, t=1)
{
  #Hardcoded arguments
  tol=1e-7; btol = 1e-7; m.max = 30; mxrej = 10 #constant
  ## R translation:  Ravi Varadhan, Johns Hopkins University
  ##		   "cosmetic", apply to sparse A: Martin Maechler, ETH Zurich
  d <- dim(A)
  # HMM constant: m <- c(1344,1344)
  n <- d[1]
  # HMM constant: n <- 1344
  if(n <= 1) {
    if(n == 1) return(exp(A*t)*v)
    stop("nrow(A) must be >= 1")
  }
  m <- min(n, m.max)  
  # HMM constant: m <- 30
  gamma <- 0.9        # constant
  delta <- 1.2        # constant
  nA <- Matrix::norm(A, "I")  # varies with Q
  # Next line varies with Q and phy
  if(nA <  1e-6) { ## rescaling, by MMaechler, needed for small norms
    A <- A/nA
    t <- t*nA
    nA <- 1
  }
  rndoff <- nA * .Machine$double.eps # Varies with Q
  
  t_1 <- abs(t) # Varies with phy
  sgn <- sign(t) # Varies with phy
  t_now <- 0
  s_error <- 0
  k1 <- 2
  mb <- m    # HMM constant: mb <- 30
  xm <- 1/m  # HMM constant: xm <- 1/30
  # Next line is constant for all tips in HMM: beta <- 2*sqrt(5)
  beta <- sqrt(sum(v*v))# = norm(v) = |\ v ||
  if(beta == 0) ## border case: v is all 0, and the result is too
    return(v)
  # Next line is constant for all tips in HMM: fact <- 6.74967950018045e+33
  fact <- (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1))
  myRound <- function(tt) {
    s <- 10^(floor(log10(tt)) - 1)
    ceiling(tt/s)*s
  }
  t_new <- myRound( (1/nA)*(fact*tol/(4*beta*nA))^xm )
  # alt for HMM, varies with Q, phy and site
  # t_new <- myRound( (1/nA)*(nA*beta)^(-1/30)*7.48584399202831 )
  # alt constant for HMM tips, varies with Q
  # t_new <- myRound( (nA)^(-31/30)*7.12126158103164 )
  
  V <- matrix(0, n, m+1)    #HMM init: V <- matrix(0,1344,31)
  H <- matrix(0, m+2, m+2)  #HMM initt: H <- matrix(0,32,32) 
  # use  Matrix(V[,j],nrow =n, ncol=1 ) later on?
  # nstep <- n.rej <- 0L      #irrelevant
  w <- v
  # updated in loop:
  # t_now, t_new, V, H, w, beta
  # updated on loop break: 
  # mb
  while (t_now < t_1) {
    # nstep <- nstep + 1L
    t_step <- min(t_1 - t_now, t_new)
    # if(verbose) cat(sprintf("while(t_now = %g < ..): nstep=%d, t_step=%g\n",
    #                         t_now, nstep, t_step))
    V[,1] <- (1/beta)*w
    for (j in 1:m) {
      p <- as.vector(A %*% V[,j])  ## as of commit ab3e84e, this %*% is ~82% of all work!
      for (i in 1:j) {
        H[i,j] <- s <- sum(V[,i] *  p)
        p <- p - s * V[,i]
      }
      s <- sqrt(sum(p*p))
      if (s < btol) {
        k1 <- 0
        mb <- j
        t_step <- t_1 - t_now
        break
      }
      H[j+1, j] <- s
      V[, j+1] <- p / s
    } ## j-loop complete
    if (k1 != 0) {
      H[m+2, m+1] <- 1
      av <- A %*% V[, m+1]  ## as of commit ab3e84e, this %*% is just ~2.7% of all work
      avnorm <- sqrt(sum(av * av))
    }
    i.rej <- 0L
    while (i.rej <= mxrej) {
      mx <- mb + k1; imx <- seq_len(mx) # = 1:mx
      # if(verbose) cat(sprintf("	inner while: k1=%d -> mx=%d\n",
      #                         k1, mx))
      F <- internal_expm(sgn * t_step * H[imx,imx, drop=FALSE])
      if (k1 == 0) {
        err_loc <- btol
        break
      } else {
        phi1 <- abs(beta * F[m+1,1])
        phi2 <- abs(beta * F[m+2,1] * avnorm)
        if(is.nan(phi1) || is.nan(phi2))
          stop("NaN phi values; probably overflow in expm()")
        if (phi1 > 10*phi2) {
          err_loc <- phi2
          xm <- 1/m
        } else if (phi1 > phi2) {
          err_loc <- (phi1 * phi2)/(phi1 - phi2)
          xm <- 1/m
        } else {
          err_loc <- phi1
          xm <- 1/(m-1)
        }
      }
      if (err_loc <= delta * t_step*tol) break
      else {
        if (i.rej == mxrej)
          stop(gettextf('The requested tolerance (tol=%g) is too small for mxrej=%d.',
                        tol, mxrej))
        t_step <- gamma * t_step * (t_step * tol / err_loc)^xm
        s <- 10^(floor(log10(t_step))-1)
        t_step <- s * ceiling(t_step / s)
        i.rej <- i.rej + 1L
      }
    }## end{ while (i.rej < mx..) }
    # n.rej <- n.rej + i.rej
    mx <- mb + max(0, k1-1); imx <- seq_len(mx) # = 1:mx
    w <- as.vector(V[, imx] %*% (beta*F[imx,1, drop=FALSE]))
    beta <- sqrt(sum(w*w))
    t_now <- t_now + t_step
    t_new <- myRound(gamma * t_step * (t_step*tol/err_loc)^xm)
    # err_loc <- max(err_loc, rndoff)
    # s_error <- s_error + err_loc
  }# end{ while }
  return(w)
}

## HMM variant of expAtv for evaluating multiple t for fixed v 
exp_A_tvec_codon <- function(A, codon, tvec=1, v=NULL, subset=NULL )
{
  #Hardcoded arguments
  tol=1e-7; btol = 1e-7; m.max = 30; mxrej = 10 #constant
  
  ## R translation:  Ravi Varadhan, Johns Hopkins University
  ##		   "cosmetic", apply to sparse A: Martin Maechler, ETH Zurich
  d <- dim(A)
  # HMM constant: m <- c(1344,1344)
  n <- d[1]
  # HMM constant: n <- 1344
  if(!missing(codon)) {
    stopifnot(n == 1344)#, "If using codon notation, A must be 1344x1344.")
    stopifnot(length(codon)==1)#, "Only one codon may be processed at a time.")
    stopifnot(codon <65 || codon >0)#, "If using codon notation, codon must be in 1:64.")
    v=Matrix::sparseVector(x=rep(1,20),i=(c(0:15,17:20)*64+codon),length = 1344)
    # Next line is constant for all tips in HMM: beta <- 2*sqrt(5)
    beta <- 2*sqrt(5)
  } else {
    stopifnot(!is.null(v))#, "If not using codon notation, v must be provided.")
    stopifnot(length(v)==n)#, "v must be the same size as a side of A.")
    beta <- sqrt(sum(v*v))# = norm(v) = |\ v ||
    
  }
  stopifnot(is.null(subset))#, "Subset notation not yet implemented.")
  
  if(n <= 1) {
    if(n == 1) return(as.list(exp(A*tvec)*v))
    stop("nrow(A) must be >= 1")
  }
  
  stopifnot(all(tvec>=0) || all(tvec<=0))# , "Mixed sign notation not yet implemented.")
  if( length(tvec)==0) return(list())
  if( length(tvec)==1) return(list(internal_expAtv(A=A,v=as.numeric(v),t=tvec)))
  res = rep(list(v),length(tvec))
  
  m <- min(n, m.max)  
  # HMM constant: m <- 30
  gamma <- 0.9        # constant
  delta <- 1.2        # constant
  nA <- Matrix::norm(A, "I")  # varies with Q
  # Next line varies with Q and phy
  if(nA <  1e-6) { ## rescaling, by MMaechler, needed for small norms
    A <- A/nA
    tvec <- tvec*nA
    nA <- 1
  }
  rndoff <- nA * .Machine$double.eps # Varies with Q
  if(all(tvec>=0)){
    t_1 = max(tvec)
    sgn=1
    t_1_vec=tvec
  } else if(all(tvec<=0)){
    t_1 = min(tvec)
    sgn=-1
    t_1_vec=-tvec
    
  } else stop("This line should be impossible to reach.")
    
  t_now <- 0
  s_error <- 0
  k1 <- 2
  mb <- m    # HMM constant: mb <- 30
  xm <- 1/m  # HMM constant: xm <- 1/30
  mx1=m+2
  mx2=m+1
  imx1=seq_len(mx1)
  imx2=seq_len(mx2)
  if(beta == 0) ## border case: v is all 0, and the result is too
    return(res)
  # Next line is constant for all tips in HMM: fact <- 6.74967950018045e+33
  fact <- (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1))
  
  myRound <- function(tt) {
    s <- 10^(floor(log10(tt)) - 1)
    ceiling(tt/s)*s
  }
  t_new <- myRound( (1/nA)*(fact*tol/(4*beta*nA))^xm )
  # alt for HMM, varies with Q, phy and site
  # t_new <- myRound( (1/nA)*(nA*beta)^(-1/30)*7.48584399202831 )
  # alt constant for HMM tips, varies with Q
  # t_new <- myRound( (nA)^(-31/30)*7.12126158103164 )
  
  V <- matrix(0, n, m+1)    #HMM init: V <- matrix(0,1344,31)
  H <- matrix(0, m+2, m+2)  #HMM initt: H <- matrix(0,32,32) 
  w <- as.numeric(v)
  # updated in loop:
  # t_now, t_new, V, H, w, beta
  # updated on loop break: 
  # mb
  while (t_now < t_1) {
    # nstep <- nstep + 1L
    subset = t_1_vec > t_now
    t_step <- min(t_1_vec[subset] - t_now, t_new)  # interval length to evaluate over
    V[,1] <- (1/beta)*w
    for (j in 1:m) {
      p <- as.vector(A %*% V[,j]) 
      for (i in 1:j) {
        H[i,j] <- s <- sum(V[,i] *  p)
        p <- p - s * V[,i]
      }
      s <- sqrt(sum(p*p))
      if (s < btol) {
        k1 <- 0
        mb <- j
        t_step <- t_1 - t_now
        break
      }
      H[j+1, j] <- s
      V[, j+1] <- p / s
    } ## j-loop complete
    if (k1 != 0) {
      H[m+2, m+1] <- 1
      av <- A %*% V[, m+1]  ## as of commit ab3e84e, this %*% is just ~2.7% of all work
      avnorm <- sqrt(sum(av * av))
    } else {  # cash out, all remaining evaluations can be performed at once!
      mx <- mb; imx = seq_len(mx)
      ivec = which(t_1_vec > t_now)
      F_list <- internal_expmt(A = H[imx,imx,drop=F],t_vec = tvec[ivec])
      res[ivec] <- lapply(F_list, function(F_mat) (beta * V[,imx] %*% F_mat[imx,1]) )
      return(res)
    }
    i.rej <- 0L
    while (i.rej <= mxrej) {
      #mx <- mb + k1; imx <- seq_len(mx) # = 1:mx
      F_edge <- internal_expm(sgn * t_step * H[imx1,imx1, drop=FALSE])
      
      # Check for adjustment due to errors
      phi1 <- abs(beta * F_edge[m+1,1])
      phi2 <- abs(beta * F_edge[m+2,1] * avnorm)
      if(is.nan(phi1) || is.nan(phi2))
        stop("NaN phi values; probably overflow in expm()")
      if (phi1 > 10*phi2) {
        err_loc <- phi2
        xm <- 1/m
      } else if (phi1 > phi2) {
        err_loc <- (phi1 * phi2)/(phi1 - phi2)
        xm <- 1/m
      } else {
        err_loc <- phi1
        xm <- 1/(m-1)
      }
      # Check if error is within tolerance
      if (err_loc <= delta * t_step*tol) break
      else {
        if (i.rej == mxrej)
          stop(gettextf('The requested tolerance (tol=%g) is too small for mxrej=%d.',
                        tol, mxrej))
        # reduce evaluation interval and re-check error estimates
        t_step <- gamma * t_step * (t_step * tol / err_loc)^xm
        s <- 10^(floor(log10(t_step))-1)
        t_step <- s * ceiling(t_step / s)
        i.rej <- i.rej + 1L
      }
    }## end{ while (i.rej < mx..) }
    # n.rej <- n.rej + i.rej
    #mx <- mb + max(0, k1-1); imx <- seq_len(mx) # = 1:mx
    # Check to see if tvec falls in this interval
    if(any(t_1_vec < t_now+t_step & t_1_vec > t_now)){
      # Due to errors in estimates, code rewritten.
      # It should never enter this block any more
      stop("Interval should never surround any tvec point, tvec should only include endpoints.")
      ivec = which(t_1_vec < t_now+t_step & t_1_vec > t_now)
      F_list <- internal_expmt(A = H[imx1,imx1,drop=F],t_vec = tvec[ivec])
      res[ivec] <- lapply(F_list, function(F_mat) as.vector(V[,imx2] %*% (beta * F_mat[imx2,1,drop=F])) )
    }
    
    # Move on to next interval 
    w <- as.vector(V[, imx2] %*% (beta*F_edge[imx2,1, drop=FALSE]))
    beta <- sqrt(sum(w*w))
    t_now <- t_now + t_step
    t_new <- myRound(gamma * t_step * (t_step*tol/err_loc)^xm)
    # Check to see if tvec is on this edge
    if(any(key <- t_1_vec == t_now))
      res[key] <- rep(list(w),sum(key))
  }# end{ while }
  return(res)
}


GetExpQt <- function(phy, Q, scale.factor, rates=NULL){

  if(!is.null(scale.factor)){
    Q.scaled = Q * (1/scale.factor)
  }else{
    Q.scaled = Q
  }
  if(!is.null(rates)){
    Q.scaled = Q.scaled * rates
  }
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  expQt <- as.list(numeric(nb.tip + nb.node))
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  #phy <- reorder(phy, "pruningwise")
  #Obtain an object of all the unique ancestors
  anc <- unique(phy$edge[,1])
  desRows <- do.call(c,lapply(anc,
                              function(focal){
                                which(phy$edge[,1]==focal)
                              }))
  
  desNodes <- phy$edge[desRows,2]
  expQt[desNodes] <-  internal_expmt(Q.scaled,phy$edge.length[desRows])
  
  return(expQt)
}


#GetExpQtParallel <- function(phy, Q, scale.factor, rates=NULL, ncores = 1){
#
#    if(!is.null(scale.factor)){
#        Q.scaled = Q * (1/scale.factor)
#    }else{
#        Q.scaled = Q
#    }

#    if(!is.null(rates)){
#        Q.scaled = Q.scaled * rates
#    }
#    nb.tip <- length(phy$tip.label)
#   nb.node <- phy$Nnode
#    expQt <- as.list(numeric(nb.tip + nb.node))
#    TIPS <- 1:nb.tip
#    comp <- numeric(nb.tip + nb.node)
#phy <- reorder(phy, "pruningwise")
#Obtain an object of all the unique ancestors
#    anc <- unique(phy$edge[,1])
#    edge.idx <- rep(NA, length(phy$edge[,1]))

#    k <- 1
#    for (i  in seq(from = 1, length.out = nb.node)) {
#the ancestral node at row i is called focal
#        focal <- anc[i]
#Get descendant information of focal
#        desRows <- which(phy$edge[,1]==focal)
#        desNodes <- phy$edge[desRows,2]
#        for (desIndex in sequence(length(desRows))){
#            edge.idx[k] <- desNodes[desIndex]
#            k <- k + 1
#        }
#    }
#    expQt <- mclapply(X = edge.idx, FUN = function(X){return(eigenExpM(Q.scaled * phy$edge.length[desRows[desIndex]]))}, mc.cores = ncores)
#    return(expQt)
#}


#Step 2: Finish likelihood by taking our already exponentiated Q down the tree and simply re-traverse the tree and multiply by the observed likelihood.
FinishLikelihoodCalculation <- function(phy, liks, Q, root.p, anc){

  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)
  if(any(root.p < 0) | any(is.na(root.p))){
    return(1000000)
  }
  #Obtain an object of all the unique ancestors
  for (i  in seq(from = 1, length.out = nb.node)) {
    #the ancestral node at row i is called focal
    focal <- anc[i]
    #Get descendant information of focal
    desRows <- which(phy$edge[,1]==focal)
    desNodes <- phy$edge[desRows,2]
    v <- 1
    for (desIndex in desNodes){
      if(desIndex <= nb.tip){
        if(sum(liks[desIndex,]) < 2){
          v <- v * (Q[[desIndex]] %*% liks[desIndex,])
        }
      }else{
        v <- v * (Q[[desIndex]] %*% liks[desIndex,])
      }
    }
    comp[focal] <- sum(v)
    liks[focal,] <- v/comp[focal]
  }
  #Specifies the root:
  root <- nb.tip + 1L
  #If any of the logs have NAs restart search:
  if(is.nan(sum(log(comp[-TIPS]))) || is.na(sum(log(comp[-TIPS])))){
    return(1000000)
  }
  else{
    loglik <- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
    if(is.infinite(loglik)){return(1000000)}
  }
  loglik
}


FinishLikelihoodCalculationHMM <- function(phy, liks, Q, root.p, anc){

  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  TIPS <- 1:nb.tip
  comp <- numeric(nb.tip + nb.node)

  if(any(root.p < 0) | any(is.na(root.p))){
    return(1000000)
  }
  #Obtain an object of all the unique ancestors
  for (i  in seq(from = 1, length.out = nb.node)) {
    #the ancestral node at row i is called focal
    focal <- anc[i]
    #Get descendant information of focal
    desRows<-which(phy$edge[,1]==focal)
    desNodes<-phy$edge[desRows,2]
    v <- 1
    for (desIndex in desNodes){
      if(desIndex <= nb.tip){
        if(sum(liks[desIndex,]) < 65){
          v <- v * (Q[[desIndex]] %*% liks[desIndex,])
        }
      }else{
        v <- v * (Q[[desIndex]] %*% liks[desIndex,])
      }
    }
    comp[focal] <- sum(v)
    liks[focal,] <- v/comp[focal]
  }
  #Specifies the root:
  root <- nb.tip + 1L
  #If any of the logs have NAs restart search:
  if(is.nan(sum(log(comp[-TIPS]))) || is.na(sum(log(comp[-TIPS])))){
    return(1000000)
  }
  else{
    loglik<- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,])))
    if(is.infinite(loglik)){return(1000000)}
  }
  loglik
}


#Step 2: Finish likelihood by taking our already exponentiated Q down the tree and simply re-traverse the tree and multiply by the observed likelihood.
#This version was rewritten by Cedric Landerer but is not any faster than the version I wrote.
#code <- '
#arma::mat m = Rcpp::as<arma::mat>(a);
#arma::vec v = Rcpp::as<arma::vec>(e);
#arma::vec out = Rcpp::as<arma::vec>(q);
#return Rcpp::wrap( out % (m * v) );
#'
#rcppArma <- cxxfunction(signature(a="numeric",e="numeric", q="numeric"), code, plugin="RcppArmadillo")

#FinishLikelihoodCalculationOLD <- function(phy, liks, Q, root.p, anc){
#   root <- length(phy$tip.label) + 1
#   nb.node <- phy$Nnode
#   nb.tip <- Ntip(phy)
#   comp <- numeric(nb.node)
#
#   if(any(root.p < 0) | any(is.na(root.p))){
#       return(10000000000)
#   }
#
#   #Obtain an object of all the unique ancestors
#   for (i in 1:nb.node) {
#       #the ancestral node at row i is called focal
#       focal <- anc[i]
#       #Get descendant information of focal
#       desRows <- which(phy$edge[,1]==focal)
#       desNodes <- phy$edge[desRows,2]
#       #desNodes <- desNodes[desNodes <= nb.tip]
#       v <- rep(1, 64)
#       for (desIndex in desNodes){
#           v <- rcppArma(Q[[desIndex]], liks[desIndex,], v)
#       }
#       comp[i] <- sum(v)
#       liks[focal,] <- v/comp[i]
#   }
#   #If any of the logs have NAs restart search:
#   log.comp <- log(comp)
#   if(any(!is.finite(log.comp))){
#       loglik <- 10000000000
#   }else{
#       loglik <- -(sum(log.comp) + log(sum(root.p * liks[root,])))
#       if(!is.finite(loglik)){
#           loglik <- 10000000000
#       }
#   }
#   return(loglik)
#}


######################################################################################################################################
######################################################################################################################################
### Likelihood calculator -- ODE solver
######################################################################################################################################

TreeTraversalODE <- function(phy, Q_codon_array_vectored, liks.HMM, bad.likelihood=-100000, root.p) {

  ##start with first method and move to next if problems encountered
  ## when solving ode, such as negative pr values < neg.pr.threshold
  ode.method.vec <- c("ode45", "lsoda")
  num.ode.method <- length(ode.method.vec)

  rtol = 1e-7 #default 1e-6 returns a negative value under long branch testing conditions
  atol = 1e-6 #default 1e-6

  neg.pr.threshold <- -10*atol


  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode

  anc <- unique(phy$edge[,1])
  TIPS <- 1:nb.tip

  comp <- numeric(nb.tip + nb.node)

  for (i in seq(from = 1, length.out = nb.node)) {
    focal <- anc[i]
    desRows <- which(phy$edge[,1]==focal) ##des = descendant
    desNodes <- phy$edge[desRows,2]
    state.pr.vector = rep(1, dim(liks.HMM)[2]) ##

    for (desIndex in sequence(length(desRows))){
      yini <- liks.HMM[desNodes[desIndex],]
      times=c(0, phy$edge.length[desRows[desIndex]])

      ode.not.solved <- TRUE
      ode.solver.attempt <- 0

      while(ode.not.solved && ode.solver.attempt < num.ode.method){
        ode.solver.attempt <- ode.solver.attempt+1
        ode.method <-  ode.method.vec[ode.solver.attempt]

        subtree.pr.ode.obj <- ode(
          y=yini, times=times, func = "selacHMM",
          parms=Q_codon_array_vectored, initfunc="initmod_selacHMM",
          dllname = "selac",
          method=ode.method, rtol=rtol, atol=atol
        )

        ## CHECK TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
        ## $istate should be = 0 [documentation in doc/deSolve.Rnw indicates
        ## it should be 2]
        ## Values < 0 indicate problems
        ## TODO: take advantage of while() around ode solving created
        ## for when we hit negative values
        istate <- attributes(subtree.pr.ode.obj)$istate[1]

        if(istate < 0){
          ## For \code{lsoda, lsodar, lsode, lsodes, vode, rk, rk4, euler} these are
          error.text <- switch(as.character(istate),
                               "-1"="excess work done",
                               "-2"="excess accuracy requested",
                               "-3"="illegal input detected",
                               "-4"="repeated error test failures",
                               "-5"="repeated convergence failures",
                               "-6"="error weight became zero",
                               paste("unknown error. ode() istate value: ", as.character(istate))
          )

          warning(print(paste("selac.R: Integration of descendent index", desIndex, ": ode solver returned state = ",  istate, " : ", error.text)))

          if(ode.solver.attempt < num.ode.method){
            warning.message <- paste("\tTrying ode method ", ode.method.vec[ode.solver.attempt+1])
            warning(warning.message)
          }else{
            warning.message <- paste(warning.message, "No additional ode methods available. Returning bad.likelihood: ", bad.likelihood)
            warning(warning.message)
            return(bad.likelihood)
          }
          #return(bad.likelihood)
        }else{
          ##no integration issues,
          ## object consists of pr values at start and end time
          ## extract final state variable, dropping time entry
          subtree.pr.vector <- subtree.pr.ode.obj[dim(subtree.pr.ode.obj)[[1]],-1]


          ## test for negative entries
          ## if encountered and less than neg.pr.threshold
          ## replace the negative values to 0
          ## if there are values less than neg.pr.threshold, then
          ## resolve equations using more robust method on the list
          ## Alternative: use 'event' option in deSolve as described at
          ## http://stackoverflow.com/questions/34424716/using-events-in-desolve-to-prevent-negative-state-variables-r
          neg.vector.pos <- which(subtree.pr.vector < 0, arr.ind=TRUE)
          num.neg.vector.pos <- length(neg.vector.pos)

          if(num.neg.vector.pos > 0){
            min.vector.val <- min(subtree.pr.vector[neg.vector.pos])
            neg.vector.pos.as.string <- toString(neg.vector.pos)

            warning.message <- paste("WARNING: subtree.pr.vector solved with ode method ", ode.method, " contains ", num.neg.vector.pos, " negative values at positions ", neg.vector.pos.as.string ,  "of a ", length(subtree.pr.vector), " vector." )


            if(min.vector.val > neg.pr.threshold){
              warning.message <- paste(warning.message, "\nMinimum value ", min.vector.val, " >  ", neg.pr.threshold, " the neg.pr.threshold.\nSetting all negative values to 0.")
              warning(warning.message)
              subtree.pr.vector[neg.vector.pos] <- 0

            }else{
              warning.message <- paste(warning.message, "selac.R: minimum value ", min.vector.val, " <  ", neg.pr.threshold, " the neg.pr.threshold.")

              if(ode.solver.attempt < num.ode.method){
                warning.message <- paste(warning.message, " Trying ode method ", ode.method.vec[ode.solver.attempt+1])
                warning(warning.message)

              }else{
                warning.message <- paste(warning.message, "No additional ode methods available. Returning bad.likelihood: ", bad.likelihood)
                warning(warning.message)
                return(bad.likelihood)
              }
            }
          }else{
            ## no negative values in pr.vs.time.matrix
            ode.not.solved <- FALSE
          }
        } ## end else to istate < 0
      } ##end while() for ode solver


      state.pr.vector <- state.pr.vector * subtree.pr.vector
    }
    comp[focal] <- sum(state.pr.vector)
    liks.HMM[focal,] <- state.pr.vector/comp[focal]
  }
  root.node <- nb.tip + 1L

  ##Check for negative transition rates
  ##mikeg:  For now, just issue warning


  neg.nodes <- which(liks.HMM[root.node,] <0)
  if(length(neg.nodes)>0){
    warning(paste("selac.R: encountered " , length(neg.nodes), " negatives values in liks.HMM[", root.node, ", ", neg.nodes, " ] =  ",  liks.HMM[root.node, neg.nodes], " at position ", i, " , desIndex ", desIndex))
  }



  loglik <- -(sum(log(comp[-TIPS])) + log(sum(root.p * liks.HMM[root.node,])))

  ##return bad.likelihood if loglik is bad
  if(!is.finite(loglik)) return(bad.likelihood)

  return(loglik)
}



#' @title Efficient optimization of the SELAC model
#' @aliases selac
#'
#' @description
#' Efficient optimization of model parameters under the SELAC model
#'
#' @param codon.data.path Provides the path to the directory containing the gene specific fasta files of coding data. Must have a ".fasta" line ending.
#' @param n.partitions The number of partitions to analyze. The order is based on the Unix order of the fasta files in the directory.
#' @param phy The phylogenetic tree to optimize the model parameters.
#' @param data.type The data type being tested. Options are "codon" or "nucleotide".
#' @param codon.model The type of codon model to use. There are four options: "none", "GY94", "YN98", "FMutSel0", "FMutSel", "selac".
#' @param edge.length Indicates whether or not edge lengths should be optimized. By default it is set to "optimize", other option is "fixed", which is the user-supplied branch lengths.
#' @param edge.linked A logical indicating whether or not edge lengths should be optimized separately for each gene. By default, a single set of each lengths is optimized for all genes.
#' @param optimal.aa Indicates what type of optimal.aa should be used. There are five options: "none", "majrule", "averaged, "optimize", or "user".
#' @param nuc.model Indicates what type nucleotide model to use. There are three options: "JC", "GTR", or "UNREST".
#' @param include.gamma A logical indicating whether or not to include a discrete gamma model.
#' @param gamma.type Indicates what type of gamma distribution to use. Options are "quadrature" after the Laguerre quadrature approach of Felsenstein 2001 or median approach of Yang 1994 or "lognormal" after a lognormal quadrature approach.
#' @param ncats The number of discrete categories.
#' @param numcode The ncbi genetic code number for translation. By default the standard (numcode=1) genetic code is used.
#' @param diploid A logical indicating whether or not the organism is diploid or not.
#' @param k.levels Provides how many levels in the polynomial. By default we assume a single level (i.e., linear).
#' @param aa.properties User-supplied amino acid distance properties. By default we assume Grantham (1974) properties.
#' @param verbose Logical indicating whether each iteration be printed to the screen.
#' @param n.cores.by.gene The number of cores to dedicate to parallelize analyses across gene.
#' @param n.cores.by.gene.by.site The number of cores to decidate to parallelize analyses by site WITHIN a gene. Note n.cores.by.gene*n.cores.by.gene.by.site is the total number of cores dedicated to the analysis.
#' @param max.tol Supplies the relative optimization tolerance.
#' @param max.tol.edges Supplies the relative optimization tolerance for branch lengths only. Default is that is the same as the max.tol.
#' @param max.evals Supplies the max number of iterations tried during optimization.
#' @param max.restarts Supplies the number of random restarts.
#' @param user.optimal.aa If optimal.aa is set to "user", this option allows for the user-input optimal amino acids. Must be a list. To get the proper order of the partitions see "GetPartitionOrder" documentation.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.
#' @param recalculate.starting.brlen Whether to use given branch lengths in the starting tree or recalculate them.
#' @param output.by.restart Logical indicating whether or not each restart is saved to a file. Default is TRUE.
#' @param output.restart.filename Designates the file name for each random restart.
#' @param user.supplied.starting.param.vals Designates user-supplied starting values for C.q.phi.Ne, Grantham alpha, and Grantham beta. Default is NULL.
#' @param tol.step If > 1, makes for coarser tolerance at earlier iterations of the optimizer
#' @param optimizer.algorithm The optimizer used by nloptr.
#' @param start.from.mle If TRUE, will start optimization from the MLE. Default is FALSE.
#' @param mle.matrix The user-supplied matrix of parameter values for when start.from.mle is set to TRUE.
#' @param partition.order Allows for a specialized order of the partitions to be gathered from the working directory.
#'
#' @details
#' Here we optimize parameters across each gene separately while keeping the shared parameters, alpha, beta, edge lengths, and nucleotide substitution parameters constant across genes. We then optimize alpha, beta, gtr, and the edge lengths while keeping the rest of the parameters for each gene fixed. This approach is potentially more efficient than simply optimizing all parameters simultaneously, especially if fitting models across 100's of genes.
#'
#' @examples
#' \dontrun{
#' phy <- ape::read.tree(file=system.file("extdata", "rokasYeast.tre", package="selac"))
#' result <- SelacOptimize(codon.data.path = paste0(find.package("selac"), '/extdata/'), n.partitions=1, phy=phy, max.evals=10)
#' print(result)
#' }
#' @export
SelacOptimize <- function(codon.data.path, n.partitions=NULL, phy, data.type="codon", codon.model="selac", edge.length="optimize", edge.linked=TRUE, optimal.aa="optimize", nuc.model="GTR", include.gamma=FALSE, gamma.type="quadrature", ncats=4, numcode=1, diploid=TRUE, k.levels=0, aa.properties=NULL, verbose=FALSE, n.cores.by.gene=1, n.cores.by.gene.by.site=1, max.tol=1e-3, max.tol.edges=1e-3, max.evals=1000000, max.restarts=3, user.optimal.aa=NULL, fasta.rows.to.keep=NULL, recalculate.starting.brlen=TRUE, output.by.restart=TRUE, output.restart.filename="restartResult", user.supplied.starting.param.vals=NULL, tol.step=1, optimizer.algorithm="NLOPT_LN_SBPLX", start.from.mle=FALSE, mle.matrix=NULL, partition.order=NULL, max.iterations=6) {

  if(!data.type == "codon" & !data.type == "nucleotide"){
    stop("Check that your data type input is correct. Options are codon or nucleotide", call.=FALSE)
  }
  if(!codon.model == "none" & !codon.model == "GY94" & !codon.model == "YN98" & !codon.model == "FMutSel0" & !codon.model == "FMutSel" & !codon.model == "selac"){
    stop("Check that your codon model is correct. Options are GY94, FMutSel0, or selac", call.=FALSE)
  }
  if(!edge.length == "optimize" & !edge.length == "fixed"){
    stop("Check that you have a supported edge length option. Options are optimize or fixed.", call.=FALSE)
  }
  if(!optimal.aa == "optimize" & !optimal.aa == "majrule" & !optimal.aa == "averaged" & !optimal.aa == "none" & !optimal.aa == "user"){
    stop("Check that you have a supported optimal amino acid option. Options are optimize, majrule, none, or user", call.=FALSE)
  }
  if(!nuc.model == "JC" & !nuc.model == "GTR" & !nuc.model == "UNREST"){
    stop("Check that you have a supported nucleotide substitution model. Options are JC, GTR, or UNREST.", call.=FALSE)
  }
  if(!gamma.type == "quadrature" & !gamma.type == "median" & !gamma.type == "lognormal"){
    stop("Check that you have a supported gamma type. Options are quadrature after Felsenstein 2001 or median after Yang 1994 or lognormal.", call.=FALSE)
  }

  if(!is.null(user.optimal.aa)){
    if(is.list(user.optimal.aa) == FALSE){
      stop("User-supplied optimal amino acids must be input as a list.", call.=FALSE)
    }
  }

  if(start.from.mle == TRUE){
    partitions <- partition.order
  }else{
    partitions <- system(paste("ls -1 ", codon.data.path, "*.fasta", sep=""), intern=TRUE)
  }

  if(is.null(n.partitions)){
    n.partitions <- length(partitions)
  }else{
    n.partitions = n.partitions
  }

  if(n.partitions<n.cores.by.gene) {
    warning(paste0("You have ", n.partitions, " partition (set with the n.partitions argument) but are asking to run across ", n.cores.by.gene, " cores, so ", n.cores.by.gene - n.partitions, " cores will not be used"))
  }

  
  #fix user-supplied branch lengths that are shorter than the prespecified minimum.
  if(recalculate.starting.brlen==FALSE){ 
  phy$edge.length[phy$edge.length <= 1e-8]<-(1e-8)*1.01
    }
  
  
#check that the taxon names are correctly formatted in comparison to the tree (DE)
{fastas=list.files(codon.data.path,pattern="*.fasta")
  if (identical(sort(unlist(names(read.FASTA(fastas[1])))) , sort(unlist(phy[[4]])))
  ) {
    print("DATA CHECK: Taxa in first alignment identical to taxa in tree. Good.")
    errorStatus<-"Safe"
  } else {
    print("Error: Taxa in alignment are not identical to taxa in tree. Check your input files. Exiting...")
    errorStatus<-"exit"
  }
  stopifnot(errorStatus!="exit")
  
  #check that the taxon names are correctly formated within all alignment files (DE)
  if(
    all(unlist(foreach(j = 1:length(fastas))%do%{                               ##makes pairwise checks of all the alignments to see if the taxon names are identical. If they are ALL identical then no error is thrown. 
      foreach(k = length(fastas):1)%do%{                                         
        identical(                                                               
          sort(names(read.FASTA(fastas[j]))), sort(names(read.FASTA(fastas[k])))
        )
      }
    }))
  ){
    print("DATA CHECK: Taxa in alignments identical to each other. Good.")
    errorStatus<-"Safe"
  } else {
    print("Error: Taxa in alignments are not identical each other. Check your input files. Exiting...")
    errorStatus<-"exit"
  }
  stopifnot(errorStatus!="exit")
} 
 
  
  #checks that the options being used make sense (DE)
  if(codon.model=="none"&optimal.aa=="none"&data.type=="codon"){
  print("You have turned off amino-acid optimization and you're not using a codon model. Please set your data-type to 'nucleotide'. Exiting...")
  errorStatus <- "exit"
  }
  stopifnot(errorStatus!="exit")
  
 ##########
  
  
  cat(paste("Using", n.cores.by.gene * n.cores.by.gene.by.site, "total processors", sep=" "), "\n")

  cat("Initializing data and model parameters...", "\n")

  site.pattern.data.list <- as.list(numeric(n.partitions))
  site.pattern.count.list <- as.list(numeric(n.partitions))
  nsites.vector <- c()
  if(optimal.aa == "none"){
    if(data.type == "nucleotide"){
      empirical.base.freq.list <- as.list(numeric(n.partitions))
      starting.branch.lengths <- matrix(0, n.partitions, length(phy$edge[,1]))
      for (partition.index in sequence(n.partitions)) {
        gene.tmp <- read.dna(partitions[partition.index], format='fasta')
        if(!is.null(fasta.rows.to.keep)){
          gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
        }else{
          gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
        }
        starting.branch.lengths[partition.index,] <- ComputeStartingBranchLengths(phy, gene.tmp, data.type=data.type, recalculate.starting.brlen=recalculate.starting.brlen)$edge.length
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
      codon.freq.by.gene.list <- as.list(numeric(n.partitions))
      empirical.aa.freq.list <- as.list(numeric(n.partitions))
      starting.branch.lengths <- matrix(0, n.partitions, length(phy$edge[,1]))
      for (partition.index in sequence(n.partitions)) {
        gene.tmp <- read.dna(partitions[partition.index], format='fasta')
        if(!is.null(fasta.rows.to.keep)){
          gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
        }else{
          gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
        }
        starting.branch.lengths[partition.index,] <- ComputeStartingBranchLengths(phy, gene.tmp, data.type=data.type, recalculate.starting.brlen=recalculate.starting.brlen)$edge.length
        codon.data <- DNAbinToCodonNumeric(gene.tmp)
        codon.data <- codon.data[phy$tip.label,]
        nsites.vector = c(nsites.vector, dim(codon.data)[2] - 1)
        aa.data <- ConvertCodonNumericDataToAAData(codon.data, numcode=numcode)
        aa.optim <- apply(aa.data[, -1], 2, GetMaxName) #starting values for all, final values for majrule
        empirical.aa.freq.list[[partition.index]] <- GetAAFreqsByGene(codon.data[,-1], aa.optim, numcode=numcode)
        codon.freq.by.gene.list[[partition.index]] <- GetCodonFreqsByGene(codon.data[,-1])
        codon.data <- SitePattern(codon.data, includes.optimal.aa=FALSE)
        site.pattern.data.list[[partition.index]] = codon.data$unique.site.patterns
        site.pattern.count.list[[partition.index]] = codon.data$site.pattern.counts
      }
    }
  }else{
    codon.freq.by.aa.list <- as.list(numeric(n.partitions))
    codon.freq.by.gene.list <- as.list(numeric(n.partitions))
    starting.branch.lengths <- matrix(0, n.partitions, length(phy$edge[,1]))
    aa.optim.list <- as.list(numeric(n.partitions))
    aa.optim.full.list <- as.list(numeric(n.partitions))
    for (partition.index in sequence(n.partitions)) {
      gene.tmp <- read.dna(partitions[partition.index], format='fasta')
      if(!is.null(fasta.rows.to.keep)){
        gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
      }else{
        gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
      }
      starting.branch.lengths[partition.index,] <- ComputeStartingBranchLengths(phy, gene.tmp, data.type=data.type, recalculate.starting.brlen=recalculate.starting.brlen)$edge.length
      codon.data <- DNAbinToCodonNumeric(gene.tmp)
      codon.data <- codon.data[phy$tip.label,]
      nsites.vector = c(nsites.vector, dim(codon.data)[2] - 1)
      aa.data <- ConvertCodonNumericDataToAAData(codon.data, numcode=numcode)
      if(optimal.aa == "user"){
        aa.optim <- user.optimal.aa[[partition.index]]
        aa.optim.full.list[[partition.index]] <- aa.optim
      }else{
        aa.optim <- apply(aa.data[, -1], 2, GetMaxName) #starting values for all, final values for majrule
        aa.optim.full.list[[partition.index]] <- aa.optim
      }
      codon.freq.by.aa.list[[partition.index]] <- GetCodonFreqsByAA(codon.data[,-1], aa.optim, numcode=numcode)
      codon.freq.by.gene.list[[partition.index]] <- GetCodonFreqsByGene(codon.data[,-1])
      aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
      colnames(aa.optim.frame.to.add) <- colnames(codon.data)
      codon.data <- rbind(codon.data, aa.optim.frame.to.add)
      codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
      site.pattern.data.list[[partition.index]] = codon.data$unique.site.patterns
      site.pattern.count.list[[partition.index]] = codon.data$site.pattern.counts
      aa.optim.list[[partition.index]] = codon.data$optimal.aa
    }
  }

  opts <- list("algorithm" = optimizer.algorithm, "maxeval" = max.evals, "ftol_rel" = max.tol)
  opts.edge <- list("algorithm" = optimizer.algorithm, "maxeval" = max.evals, "ftol_rel" = max.tol.edges)


  results.final <- c()
  if(nuc.model == "JC"){
    nuc.ip = NULL
    max.par.model.count = 0
    parameter.column.names <- c()
  }
  if(nuc.model == "GTR"){
    nuc.ip = rep(1, 5)
    max.par.model.count = 5
    parameter.column.names <- c("C_A", "G_A", "T_A", "G_C", "T_C")
  }
  if(nuc.model == "UNREST"){
    nuc.ip = rep(1, 11)
    max.par.model.count = 11
    parameter.column.names <- c("C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T")
  }

  if(optimal.aa=="none") {
    if(data.type == "nucleotide"){
      codon.index.matrix = NA
      if(include.gamma == TRUE){
        ip = c(1,nuc.ip)
        upper = c(5, rep(21, length(ip)-1))
        lower = rep(-21, length(ip))
        max.par.model.count = max.par.model.count + 1
        parameter.column.names <- c("shape.gamma", parameter.column.names)
      }else{
        ip = nuc.ip
        upper = rep(21, length(ip))
        lower = rep(-21, length(ip))
      }
      index.matrix = matrix(0, n.partitions, length(ip))
      index.matrix[1,] = 1:ncol(index.matrix)
      ip.vector = ip
      upper.vector = upper
      lower.vector = lower
      if(n.partitions > 1){
        for(partition.index in 2:n.partitions){
          ip.vector = c(ip.vector, ip)
          upper.vector = c(upper.vector, upper)
          lower.vector = c(lower.vector, lower)
          index.matrix.tmp = numeric(max.par.model.count)
          index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
          index.matrix[partition.index,] <- index.matrix.tmp
        }
      }
      number.of.current.restarts <- 1
      best.lik <- 1000000
      while(number.of.current.restarts < (max.restarts+1)){
        cat(paste("Finished. Performing analysis...", sep=""), "\n")
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
        if(edge.length == "optimize"){
          cat("       Optimizing edge lengths", "\n")
          phy$edge.length <- colMeans(starting.branch.lengths)
          #opts.edge <- opts
          upper.edge <- rep(log(10), length(phy$edge.length))
          lower.edge <- rep(log(1e-8), length(phy$edge.length))
          results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=empirical.base.freq.list, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          print(results.edge.final$objective)
          print(exp(results.edge.final$solution))
          phy$edge.length <- exp(results.edge.final$solution)
        }
        cat("       Optimizing model parameters", "\n")
        ParallelizedOptimizedByGene <- function(n.partition){
          optim.by.gene <- nloptr(x0=log(mle.pars.mat[n.partition,]), eval_f = OptimizeModelParsLarge, ub=upper.vector[1:dim(mle.pars.mat)[2]], lb=lower.vector[1:dim(mle.pars.mat)[2]], opts=opts, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=empirical.base.freq.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, neglnl=TRUE)
          tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
          return(tmp.pars)
        }
        results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
        parallelized.parameters <- t(matrix(unlist(results.set),dim(index.matrix)[2]+1,n.partitions))
        results.final <- NULL
        results.final$objective <- sum(parallelized.parameters[,1])
        results.final$solution <- c(t(parallelized.parameters[,-1]))
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
        print(results.final$objective)
        print(mle.pars.mat)

        current.likelihood <- results.final$objective
        cat(paste("Current likelihood", current.likelihood, sep=" "), "\n")
        lik.diff <- 10
        iteration.number <- 1
        while(lik.diff != 0 & iteration.number<=max.iterations){
          cat(paste("Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
          if(edge.length == "optimize"){
            cat("       Optimizing edge lengths", "\n")
            #opts.edge <- opts
            opts.edge$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))

            results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=empirical.base.freq.list, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            print(results.edge.final$objective)
            print(exp(results.edge.final$solution))
            phy$edge.length <- exp(results.edge.final$solution)
          }
          cat("       Optimizing model parameters", "\n")
          opts.params <- opts
          opts.params$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))

          ParallelizedOptimizedByGene <- function(n.partition){
            optim.by.gene <- nloptr(x0=log(mle.pars.mat[n.partition,]), eval_f = OptimizeModelParsLarge, ub=upper.vector[1:dim(mle.pars.mat)[2]], lb=lower.vector[1:dim(mle.pars.mat)[2]], opts=opts.params, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=empirical.base.freq.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, neglnl=TRUE)
            tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
            return(tmp.pars)
          }
          results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
          parallelized.parameters <- t(matrix(unlist(results.set),dim(index.matrix)[2]+1,n.partitions))
          results.final <- NULL
          results.final$objective <- sum(parallelized.parameters[,1])
          results.final$solution <- c(t(parallelized.parameters[,-1]))
          mle.pars.mat <- index.matrix
          mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
          print(results.final$objective)
          print(mle.pars.mat)

          lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
          current.likelihood <- results.final$objective
          cat(paste("Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
          iteration.number <- iteration.number + 1
        }
        #Output for use in sims#
        if(output.by.restart == TRUE){
          obj.tmp = list(np=max(index.matrix) + length(phy$edge.length) + sum(nsites.vector), loglik = -results.final$objective, AIC = -2*(-results.final$objective)+2*(max(index.matrix) + length(phy$edge.length) + sum(nsites.vector)), mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=NULL, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, empirical.base.freqs=empirical.base.freq.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=ip.vector)
          class(obj.tmp) = "selac"
          save(obj.tmp,file=paste(paste(codon.data.path, output.restart.filename, sep=""), number.of.current.restarts, "Rsave", sep="."))
        }
        ########################
        if(results.final$objective < best.lik){
          best.ip <- ip.vector
          best.lik <- results.final$objective
          best.solution <- mle.pars.mat
          best.edge.lengths <- phy$edge.length
        }
        number.of.current.restarts <- number.of.current.restarts + 1
      }

      loglik <- -(best.lik) #to go from neglnl to lnl
      mle.pars.mat <- best.solution
      if(edge.length == "optimize"){
        phy$edge.length <- best.edge.lengths
      }
      cat("Finished. Summarizing results...", "\n")
      colnames(mle.pars.mat) <- parameter.column.names

      if(edge.length == "optimize"){
        np <- max(index.matrix) + length(phy$edge.length)
      }else{
        np <- max(index.matrix)
      }
      obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, AICc = NULL, mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, data.type=data.type, codon.model=codon.model, nsites=nsites.vector, aa.optim=NULL, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, empirical.base.freqs=empirical.base.freq.list, max.tol=max.tol, max.evals=max.evals)
      class(obj) = "selac"
    }else{
      if(codon.model == "GY94"){
        max.par.model.count <- 2
        ip = c(1,1)
        parameter.column.names <- c("V", "kappa")
        upper = rep(log(99), length(ip))
        lower = rep(-21, length(ip))

        codon.index.matrix = NA

        index.matrix = matrix(0, n.partitions, length(ip))
        index.matrix[1,] = 1:ncol(index.matrix)
        ip.vector = ip
        upper.vector = upper
        lower.vector = lower
        if(n.partitions > 1){
            for(partition.index in 2:n.partitions){
                #ip.vector = c(ip.vector, 1)
                #upper.vector = c(upper.vector, log(99))
                #lower.vector = c(lower.vector, -10)
                #index.matrix.tmp = numeric(max.par.model.count)
                #index.matrix.tmp[2] = 2
                #index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
                #index.matrix[partition.index,] <- index.matrix.tmp
                index.matrix[partition.index,] <- 1:ncol(index.matrix)
            }
        }
      }
      if(codon.model == "YN98"){
          max.par.model.count <- 2
          ip = c(1,1)
          parameter.column.names <- c("omega", "kappa")
          upper = rep(log(99), length(ip))
          lower = rep(-21, length(ip))

          codon.index.matrix = NA

          index.matrix = matrix(0, n.partitions, length(ip))
          index.matrix[1,] = 1:ncol(index.matrix)
          ip.vector = ip
          upper.vector = upper
          lower.vector = lower
          if(n.partitions > 1){
              for(partition.index in 2:n.partitions){
                  #ip.vector = c(ip.vector, 1)
                  #upper.vector = c(upper.vector, log(99))
                  #lower.vector = c(lower.vector, -10)
                  #index.matrix.tmp = numeric(max.par.model.count)
                  #index.matrix.tmp[2] = 2
                  #index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
                  #index.matrix[partition.index,] <- index.matrix.tmp
                  index.matrix[partition.index,] <- 1:ncol(index.matrix)
              }
          }
      }
      if(codon.model == "FMutSel0"){
          empirical.aa.freq.unlist <- matrix(unlist(empirical.aa.freq.list), ncol = 21, byrow = TRUE)
          empirical.aa.freq <- colSums(empirical.aa.freq.unlist)/ sum(colSums(empirical.aa.freq.unlist))
          fitness.pars <- GetFitnessStartingValues(codon.freqs=empirical.aa.freq)[-c(17,21)]
          aa.ordered <- .unique.aa
          aa.ordered <- aa.ordered[-c(17,21)]
          if(nuc.model == "UNREST"){
              max.par.model.count <- max.par.model.count + 1 + 19
              ip = c(0.4, nuc.ip, fitness.pars)
              parameter.column.names <- c("omega", parameter.column.names, paste("fitness", aa.ordered, sep="_"))
              upper = rep(log(99), length(ip))
              lower = rep(-10, length(ip))
          }else{
              max.par.model.count <- max.par.model.count + 3 + 1 + 19
              ip = c(0.4, .25, .25, .25, nuc.ip, fitness.pars)
              parameter.column.names <- c("omega", "freqA", "freqC", "freqG", parameter.column.names, paste("fitness", aa.ordered, sep="_"))
              upper = c(log(99), 0, 0, 0, rep(log(99), length(ip)-4))
              lower = rep(-10, length(ip))
          }

          codon.index.matrix = NA

          index.matrix = matrix(0, n.partitions, length(ip))
          index.matrix[1,] = 1:ncol(index.matrix)
          ip.vector = ip
          upper.vector = upper
          lower.vector = lower
          if(n.partitions > 1){
              for(partition.index in 2:n.partitions){
                  #ip.vector = c(ip.vector, 0.4)
                  #upper.vector = c(upper.vector, log(99))
                  #lower.vector = c(lower.vector, -10)
                  #index.matrix.tmp = numeric(max.par.model.count)
                  #index.matrix.tmp[2:max.par.model.count] = 2:max.par.model.count
                  #index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
                  #index.matrix[partition.index,] <- index.matrix.tmp
                  index.matrix[partition.index,] <- 1:ncol(index.matrix)
              }
          }
      }
      if(codon.model == "FMutSel"){
          empirical.codon.freq.unlist <- matrix(unlist(codon.freq.by.gene.list), ncol = 64, byrow = TRUE)
          empirical.codon.freq <- colSums(empirical.codon.freq.unlist)/ sum(colSums(empirical.codon.freq.unlist))
          fitness.pars <- GetFitnessStartingValues(codon.freqs=empirical.codon.freq, n.pars=64)
          codon.ordered <- .codon.name
          codon.ordered <- codon.ordered[-c(49,51,57,64)]
          if(nuc.model == "UNREST"){
              max.par.model.count <- max.par.model.count + 1 + 60
              ip = c(0.4, nuc.ip, fitness.pars)
              parameter.column.names <- c("omega", parameter.column.names, paste("fitness", codon.ordered, sep="_"))
              upper = c(rep(log(99), length(ip)-3))
              lower = rep(-10, length(ip))
          }else{
              max.par.model.count <- max.par.model.count + 3 + 1 + 60
              ip = c(0.4, .25, .25, .25, nuc.ip, fitness.pars)
              parameter.column.names <- c("omega", "freqA", "freqC", "freqG", parameter.column.names, paste("fitness", codon.ordered, sep="_"))
              upper = c(log(99), 0, 0, 0, rep(log(99), length(ip)-4))
              lower = rep(-10, length(ip))
          }

          codon.index.matrix = NA

          index.matrix = matrix(0, n.partitions, length(ip))
          index.matrix[1,] = 1:ncol(index.matrix)
          ip.vector = ip
          upper.vector = upper
          lower.vector = lower
          if(n.partitions > 1){
              for(partition.index in 2:n.partitions){
                  #ip.vector = c(ip.vector, 0.4)
                  #upper.vector = c(upper.vector, log(99))
                  #lower.vector = c(lower.vector, -10)
                  #index.matrix.tmp = numeric(max.par.model.count)
                  #index.matrix.tmp[2:max.par.model.count] = 2:max.par.model.count
                  #index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
                  #index.matrix[partition.index,] <- index.matrix.tmp
                  index.matrix[partition.index,] <- 1:ncol(index.matrix)
              }
          }
      }

      number.of.current.restarts <- 1
      best.lik <- 10000000
      while(number.of.current.restarts < (max.restarts+1)){
        cat(paste("Finished. Performing analysis...", sep=""), "\n")
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
        print(mle.pars.mat)
        if(edge.length == "optimize"){
          cat("       Optimizing edge lengths", "\n")
          phy$edge.length <- colMeans(starting.branch.lengths)
          #opts.edge <- opts
          upper.edge <- rep(log(10), length(phy$edge.length))
          lower.edge <- rep(log(1e-8), length(phy$edge.length))
          results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          print(results.edge.final$objective)
          print(exp(results.edge.final$solution))
          phy$edge.length <- exp(results.edge.final$solution)
        }
        cat("       Optimizing model parameters", "\n")
        #ParallelizedOptimizedByGene <- function(n.partition){
        optim.by.gene <- nloptr(x0=log(ip.vector), eval_f = OptimizeModelParsLarge, ub=upper.vector, lb=lower.vector, opts=opts, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, neglnl=TRUE)
          #tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
          #return(tmp.pars)
          #}
        #results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
        #results.set <- lapply(1:n.partitions, ParallelizedOptimizedByGene)
        #parallelized.parameters <- t(matrix(unlist(results.set),dim(index.matrix)[2]+1,n.partitions))
        #results.final <- NULL
        #results.final$objective <- sum(parallelized.parameters[,1])
        #results.final$solution <- c(t(parallelized.parameters[,-1]))
        results.final$objective <- optim.by.gene$objective
        results.final$solution <- optim.by.gene$solution
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
        print(results.final$objective)
        print(mle.pars.mat)

        current.likelihood <- results.final$objective
        cat(paste("Current likelihood", current.likelihood, sep=" "), "\n")
        lik.diff <- 10
        iteration.number <- 1
        while(lik.diff != 0 & iteration.number<=max.iterations){
          cat(paste("Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
          if(edge.length == "optimize"){
            cat("       Optimizing edge lengths", "\n")
            #opts.edge <- opts
            opts.edge$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))

            results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            print(results.edge.final$objective)
            print(exp(results.edge.final$solution))
            phy$edge.length <- exp(results.edge.final$solution)
          }
          cat("       Optimizing model parameters", "\n")
          opts.params <- opts
          opts.params$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
          print(length(results.final$solution))
          print(results.final$solution)
          #ParallelizedOptimizedByGene <- function(n.partition){
          optim.by.gene <- nloptr(x0=results.final$solution, eval_f = OptimizeModelParsLarge, ub=upper.vector, lb=lower.vector, opts=opts, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=NULL, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, neglnl=TRUE)
          #tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
          # return(tmp.pars)
          #}
          #results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
          #results.set <- lapply(1:n.partitions, ParallelizedOptimizedByGene)
          #parallelized.parameters <- t(matrix(unlist(results.set),dim(index.matrix)[2]+1,n.partitions))
          #results.final <- NULL
          #results.final$objective <- sum(parallelized.parameters[,1])
          #results.final$solution <- c(t(parallelized.parameters[,-1]))
          results.final$objective <- optim.by.gene$objective
          results.final$solution <- optim.by.gene$solution
          mle.pars.mat <- index.matrix
          mle.pars.mat[] <- c(exp(results.final$solution), 0)[index.matrix]
          print(results.final$objective)
          print(mle.pars.mat)

          lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
          current.likelihood <- results.final$objective
          cat(paste("Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
          iteration.number <- iteration.number + 1
        }
        #Output for use in sims#
        if(output.by.restart == TRUE){
          obj.tmp = list(np=max(index.matrix) + length(phy$edge.length) + sum(nsites.vector), loglik = -results.final$objective, AIC = -2*(-results.final$objective)+2*(max(index.matrix) + length(phy$edge.length) + sum(nsites.vector)), mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=NULL, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=NULL, empirical.aa.freqs=empirical.aa.freq.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=ip.vector)
          class(obj.tmp) = "selac"
          save(obj.tmp,file=paste(paste(codon.data.path, output.restart.filename, sep=""), number.of.current.restarts, "Rsave", sep="."))
        }
        ########################
        if(results.final$objective < best.lik){
          best.ip <- ip.vector
          best.lik <- results.final$objective
          best.solution <- mle.pars.mat
          best.edge.lengths <- phy$edge.length
        }
        number.of.current.restarts <- number.of.current.restarts + 1
      }
      loglik <- -(best.lik) #to go from neglnl to lnl
      mle.pars.mat <- best.solution
      if(edge.length == "optimize"){
        phy$edge.length <- best.edge.lengths
      }
      cat("Finished. Summarizing results...", "\n")
      colnames(mle.pars.mat) <- parameter.column.names

      if(edge.length == "optimize"){
        np <- max(index.matrix) + length(phy$edge.length)
      }else{
        np <- max(index.matrix)
      }
      obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, AICc = NULL, mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, data.type=data.type, codon.model=codon.model, nsites=nsites.vector, aa.optim=NULL, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, empirical.aa.freqs=empirical.aa.freq.list, max.tol=max.tol, max.evals=max.evals)
      class(obj) = "selac"
    }
  }
  if(optimal.aa=="majrule" | optimal.aa=="optimize" | optimal.aa=="averaged" | optimal.aa=="user") {
    codon.index.matrix = CreateCodonMutationMatrixIndex()
    cpv.starting.parameters <- GetAADistanceStartingParameters(aa.properties=aa.properties)
    if(max.restarts > 1){
      selac.starting.vals <- matrix(0, max.restarts+1, 3)
      selac.starting.vals[,1] <- runif(n = max.restarts+1, min = (10^-10)*5e6, max = (10^-5)*5e6)
      selac.starting.vals[,2] <- runif(n = max.restarts+1, min = 0.01, max = 3)
      selac.starting.vals[,3] <- runif(n = max.restarts+1, min = 0.01, max = 1)
    }else{
      if(is.null(user.supplied.starting.param.vals)){
        selac.starting.vals <- matrix(c(2, 1.8292716544, 0.1017990371), 1, 3)
        selac.starting.vals <- rbind(selac.starting.vals, c(2, 1.8292716544, 0.1017990371))
      }else{
        selac.starting.vals <- matrix(c(user.supplied.starting.param.vals[1], user.supplied.starting.param.vals[2], user.supplied.starting.param.vals[3]), 1, 3)
        selac.starting.vals <- rbind(selac.starting.vals, c(user.supplied.starting.param.vals[1], user.supplied.starting.param.vals[2], user.supplied.starting.param.vals[3]))
      }
    }
    if(include.gamma == TRUE){
      #Gamma variation is turned ON:
      if(nuc.model == "JC"){
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "shape.gamma")
          upper = c(log(50),  21, 21, 0, 0, 0, 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 1
        }else{
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "shape.gamma")
          upper = c(log(50), 21, 21, 0, 0, 0, 10, 10, 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 2 + 1
        }
      }
      if(nuc.model == "GTR") {
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "C_A", "G_A", "T_A", "G_C", "T_C", "shape.gamma")
          upper = c(log(50), 21, 21, 0, 0, 0, rep(21, length(nuc.ip)), 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 1
        }else{
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, nuc.ip, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "C_A", "G_A", "T_A", "G_C", "T_C", "shape.gamma")
          upper = c(log(50), 21, 21, 0, 0, 0, 10, 10, rep(21, length(nuc.ip)), 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 2	+ 1
        }
      }
      if(nuc.model == "UNREST") {
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], nuc.ip, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "shape.gamma")
          upper = c(log(50), 21, 21, rep(21, length(nuc.ip)), 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 1
        }else{
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 1, 1, nuc.ip, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "a0", "a1", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "shape.gamma")
          upper = c(log(50), 21, 21, 10, 10, rep(21, length(nuc.ip)), 5)
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 2 + 1
        }
      }
      index.matrix = matrix(0, n.partitions, max.par.model.count)
      index.matrix[1,] = 1:ncol(index.matrix)
      ip.vector = ip
      if(n.partitions > 1){
        # Gamma variation is turned ON:
        for(partition.index in 2:n.partitions){
          if(nuc.model == "JC"){
            #ip.vector = c(ip.vector, ip[1], ip[8])
            if(start.from.mle == TRUE){
              ip.vector = c(ip.vector, mle.matrix[partition.index,1])
            }else{
              ip.vector = c(ip.vector, ip[1])
            }
          }else{
            if(nuc.model == "GTR"){
              index.matrix.tmp = numeric(max.par.model.count)
              if(k.levels == 0){
                #index.matrix.tmp[c(2:11)] = c(2:11)
                #ip.vector = c(ip.vector, ip[1], ip[12])
                index.matrix.tmp[c(2:12)] = c(2:12)
                if(start.from.mle == TRUE){
                  ip.vector = c(ip.vector, mle.matrix[partition.index,1])
                }else{
                  ip.vector = c(ip.vector, ip[1])
                }
              }else{
                #index.matrix.tmp[c(2:13)] = c(2:13)
                #ip.vector = c(ip.vector, ip[1], ip[14])
                index.matrix.tmp[c(2:14)] = c(2:14)
                if(start.from.mle == TRUE){
                  ip.vector = c(ip.vector, mle.matrix[partition.index,1])
                }else{
                  ip.vector = c(ip.vector, ip[1])
                }
              }
            }else{
              index.matrix.tmp = numeric(max.par.model.count)
              if(k.levels == 0){
                #index.matrix.tmp[c(2:14)] = c(2:14)
                #ip.vector = c(ip.vector, ip[1], ip[15])
                index.matrix.tmp[c(2:15)] = c(2:15)
                if(start.from.mle == TRUE){
                  ip.vector = c(ip.vector, mle.matrix[partition.index,1])
                }else{
                  ip.vector = c(ip.vector, ip[1])
                }
              }else{
                #index.matrix.tmp[c(2:16)] = c(2:16)
                #ip.vector = c(ip.vector, ip[1], ip[17])
                index.matrix.tmp[c(2:17)] = c(2:17)
                if(start.from.mle == TRUE){
                  ip.vector = c(ip.vector, ip[partition.index,1])
                }else{
                  ip.vector = c(ip.vector, ip[1])
                }
              }
            }
          }
          index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
          index.matrix[partition.index,] <- index.matrix.tmp
        }
      }
    }else{
      # Gamma variation is turned OFF:
      if(nuc.model == "JC"){
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25)
          }
          parameter.column.names <- c("C.q.phi.Ne",  "alpha", "beta", "freqA", "freqC", "freqG")
          upper = c(log(50), 21, 21, 0, 0, 0)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 0
        }else{
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1")
          upper = c(log(50), 21, 21, 0, 0, 0, 10, 10)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 0 + 2
        }
      }
      if(nuc.model == "GTR") {
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "C_A", "G_A", "T_A", "G_C", "T_C")
          upper = c(log(50), 21, 21, 0, 0, 0, rep(21, length(nuc.ip)))
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 0
        }else{
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, nuc.ip)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "C_A", "G_A", "T_A", "G_C", "T_C")
          upper = c(log(50), 21, 21, 0, 0, 0, 10, 10, rep(21, length(nuc.ip)))
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 0 + 2
        }
      }
      if(nuc.model == "UNREST") {
        if(k.levels == 0){
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], nuc.ip)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T")
          upper = c(log(50), 21, 21, rep(21, length(nuc.ip)))
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 0
        }else{
          if(start.from.mle == TRUE){
            ip = mle.matrix[1,]
          }else{
            ip = c(selac.starting.vals[1,1:3], 1, 1, nuc.ip)
          }
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "a0", "a1", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T")
          upper = c(log(50), 21, 21, 10, 10, rep(21, length(nuc.ip)))
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 0 + 2
        }
      }
      index.matrix = matrix(0, n.partitions, max.par.model.count)
      index.matrix[1,] = 1:ncol(index.matrix)
      ip.vector = ip
      if(n.partitions > 1){
        for(partition.index in 2:n.partitions){
          if(nuc.model == "JC"){
            if(start.from.mle == TRUE){
              ip.vector = c(ip.vector, mle.matrix[partition.index,1])
            }else{
              ip.vector = c(ip.vector, ip[1])
            }
          }else{
            if(nuc.model == "GTR"){
              if(start.from.mle == TRUE){
                ip.vector = c(ip.vector, mle.matrix[partition.index,1])
              }else{
                ip.vector = c(ip.vector, ip[1])
              }
              index.matrix.tmp = numeric(max.par.model.count)
              if(k.levels == 0){
                index.matrix.tmp[c(2:11)] = c(2:11)
              }else{
                index.matrix.tmp[c(2:13)] = c(2:13)
              }
            }else{
              if(start.from.mle == TRUE){
                ip.vector = c(ip.vector, mle.matrix[partition.index,1])
              }else{
                ip.vector = c(ip.vector, ip[1])
              }
              index.matrix.tmp = numeric(max.par.model.count)
              if(k.levels == 0){
                index.matrix.tmp[c(2:14)] = c(2:14)
              }else{
                index.matrix.tmp[c(2:16)] = c(2:16)
              }
            }
          }
          index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
          index.matrix[partition.index,] <- index.matrix.tmp
        }
      }
    }

    #THIS IS FOR THERE IS A SEPARATE GAMMA PER GENE:
    #if(include.gamma == TRUE){
    #    index.matrix.red <- t(matrix(1:(n.partitions*2), 2, n.partitions))
    #}else{
    #    index.matrix.red <- t(matrix(1:n.partitions, 1, n.partitions))
    #}

    #This is so we can break out alpha, beta, GTR, and gamma which are shared among ALL genes:
    index.matrix.red <- t(matrix(1:n.partitions, 1, n.partitions))

    if(optimal.aa == "optimize"){
      number.of.current.restarts <- 1
      aa.optim.original <- aa.optim.list
      best.lik <- 1000000
      while(number.of.current.restarts < (max.restarts+1)){
        cat(paste("Finished. Performing random restart ", number.of.current.restarts,"...", sep=""), "\n")
        aa.optim.list <- aa.optim.original
        cat("       Doing first pass using majority-rule optimal amino acid...", "\n")
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
        print(mle.pars.mat)
        if(edge.length == "optimize"){
          cat("              Optimizing edge lengths", "\n")
          phy$edge.length <- colMeans(starting.branch.lengths)
          #phy$edge.length <- colMeans(starting.branch.lengths) / (1/selac.starting.vals[number.of.current.restarts, 2])
          #opts.edge <- opts
          upper.edge <- rep(log(50), length(phy$edge.length))
          lower.edge <- rep(log(1e-8), length(phy$edge.length))
          results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=NULL, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          print(results.edge.final$objective)
          print(exp(results.edge.final$solution))
          phy$edge.length <- exp(results.edge.final$solution)
        }
        cat("              Optimizing model parameters", "\n")

        #if(include.gamma == TRUE){
        #    alpha.beta.gtr <- mle.pars.mat[1,c(2:(max.par.model.count-1))]
        #    upper.bounds.shared <- upper[c(2:(max.par.model.count-1))]
        #    lower.bounds.shared <- lower[c(2:(max.par.model.count-1))]
        #}else{
        alpha.beta.gtr <- mle.pars.mat[1,c(2:max.par.model.count)]
        upper.bounds.shared <- upper[c(2:max.par.model.count)]
        lower.bounds.shared <- lower[c(2:max.par.model.count)]
        #}

        ParallelizedOptimizedByGene <- function(n.partition){
          #if(include.gamma == TRUE){
          #    tmp.par.mat <- mle.pars.mat[,c(1, max.par.model.count)]
          #    upper.bounds.gene <- upper[c(1,max.par.model.count)]
          #    lower.bounds.gene <- lower[c(1,max.par.model.count)]
          #}else{
          tmp.par.mat <- as.matrix(mle.pars.mat[,1])
          upper.bounds.gene <- upper[1]
          lower.bounds.gene <- lower[1]
          #}
          optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f = OptimizeModelParsAlphaBetaGtrFixed, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, alpha.beta.gtr=alpha.beta.gtr, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, aa.optim_array=aa.optim.list[[n.partition]], codon.freq.by.aa=codon.freq.by.aa.list[[n.partition]], codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
          return(tmp.pars)
        }
        results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
        #if(include.gamma == TRUE){
        #The number of columns is 3: [1] log-likelihood, [2] C.q.phi, [3] phi gamma:
        #parallelized.parameters <- t(matrix(unlist(results.set), 3, n.partitions))
        #}else{
        #The number of columns is 2: [1] log-likelihood, [2] C.q.phi:
        parallelized.parameters <- t(matrix(unlist(results.set), 2, n.partitions))
        #}
        results.final <- NULL
        results.final$objective <- sum(parallelized.parameters[,1])
        results.final$solution <- c(t(parallelized.parameters[,-1]))
        mle.pars.mat.red <- index.matrix.red
        mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]
        print(mle.pars.mat.red)
        optim.alpha.beta.gtr.all.genes <- nloptr(x0=log(alpha.beta.gtr), eval_f = OptimizeAlphaBetaGtrOnly, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix.red, phy=phy, aa.optim_array=aa.optim.list, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
        results.final$objective <- optim.alpha.beta.gtr.all.genes$objective
        alpha.beta.gtr <- exp(optim.alpha.beta.gtr.all.genes$solution)
        #if(include.gamma == TRUE){
        #    mle.pars.mat <- c()
        #    for(row.index in 1:dim(mle.pars.mat.red)[1]){
        #        mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2]))
        #    }
        #}else{
        mle.pars.mat <- c()
        for(row.index in 1:dim(mle.pars.mat.red)[1]){
          mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr))
        }
        #}
        print(results.final$objective)
        print(mle.pars.mat)

        current.likelihood <- results.final$objective
        cat(paste("       Current likelihood", current.likelihood, sep=" "), "\n")
        lik.diff <- 10
        iteration.number <- 1
        while(lik.diff != 0 & iteration.number <= max.iterations){
          cat(paste("       Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
          cat("              Optimizing amino acids", "\n")
          aa.optim.list <- as.list(numeric(n.partitions))
          ParallelizedOptimizeAAByGene <- function(n.partition){
            gene.tmp <- read.dna(partitions[n.partition], format='fasta')
            if(!is.null(fasta.rows.to.keep)){
              gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
            }else{
              gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
            }
            codon.data <- DNAbinToCodonNumeric(gene.tmp)
            codon.data <- codon.data[phy$tip.label,]
            tmp.aa.optim.full <- GetOptimalAAPerSite(x=log(mle.pars.mat[n.partition,]), codon.data=codon.data, phy=phy, aa.optim_array=aa.optim.list[[n.partition]], codon.freq.by.aa=codon.freq.by.aa.list[[n.partition]], codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, neglnl=TRUE, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
            return(tmp.aa.optim.full)
          }
          aa.optim.full.list <- mclapply(1:n.partitions, ParallelizedOptimizeAAByGene, mc.cores=n.cores.by.gene)
          for(partition.index in sequence(n.partitions)) {
            gene.tmp <- read.dna(partitions[partition.index], format='fasta')
            if(!is.null(fasta.rows.to.keep)){
              gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
            }else{
              gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
            }
            codon.data <- DNAbinToCodonNumeric(gene.tmp)
            codon.data <- codon.data[phy$tip.label,]
            codon.freq.by.aa.list[[partition.index]] <- GetCodonFreqsByAA(codon.data[,-1], aa.optim.full.list[[partition.index]], numcode=numcode)
            aa.optim.frame.to.add <- matrix(c("optimal", aa.optim.full.list[[partition.index]]), 1, dim(codon.data)[2])
            colnames(aa.optim.frame.to.add) <- colnames(codon.data)
            codon.data <- rbind(codon.data, aa.optim.frame.to.add)
            codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
            site.pattern.data.list[[partition.index]] = codon.data$unique.site.patterns
            site.pattern.count.list[[partition.index]] = codon.data$site.pattern.counts
            aa.optim.list[[partition.index]] = codon.data$optimal.aa
          }
          if(edge.length == "optimize"){
            cat("              Optimizing edge lengths", "\n")
            #opts.edge <- opts
            opts.edge$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))

            results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=NULL, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            print(results.edge.final$objective)
            print(exp(results.edge.final$solution))
            phy$edge.length <- exp(results.edge.final$solution)
          }
          cat("              Optimizing model parameters", "\n")

          ParallelizedOptimizedByGene <- function(n.partition){
            #if(include.gamma == TRUE){
            #    tmp.par.mat <- mle.pars.mat[,c(1, max.par.model.count)]
            #    upper.bounds.gene <- upper[c(1,max.par.model.count)]
            #    lower.bounds.gene <- lower[c(1,max.par.model.count)]
            #}else{
            tmp.par.mat <- as.matrix(mle.pars.mat[,1])
            upper.bounds.gene <- upper[1]
            lower.bounds.gene <- lower[1]
            #}
            opts.params <- opts
            opts.params$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
            optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f=OptimizeModelParsAlphaBetaGtrFixed, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts.params, alpha.beta.gtr=alpha.beta.gtr, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, aa.optim_array=aa.optim.list[[n.partition]], codon.freq.by.aa=codon.freq.by.aa.list[[n.partition]], codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
            return(tmp.pars)
          }
          results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
          #if(include.gamma == TRUE){
          #The number of columns is 3: [1] log-likelihood, [2] C.q.phi.Ne, [3] phi gamma:
          #parallelized.parameters <- t(matrix(unlist(results.set), 3, n.partitions))
          #}else{
          #The number of columns is 2: [1] log-likelihood, [2] C.q.phi.Ne:
          parallelized.parameters <- t(matrix(unlist(results.set), 2, n.partitions))
          #}
          results.final <- NULL
          results.final$objective <- sum(parallelized.parameters[,1])
          results.final$solution <- c(t(parallelized.parameters[,-1]))
          mle.pars.mat.red <- index.matrix.red
          mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]

          optim.alpha.beta.gtr.all.genes <- nloptr(x0=log(alpha.beta.gtr), eval_f = OptimizeAlphaBetaGtrOnly, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix.red, phy=phy, aa.optim_array=aa.optim.list, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          results.final$objective <- optim.alpha.beta.gtr.all.genes$objective
          alpha.beta.gtr <- exp(optim.alpha.beta.gtr.all.genes$solution)
          #if(include.gamma == TRUE){
          #    mle.pars.mat <- c()
          #    for(row.index in 1:dim(mle.pars.mat.red)[1]){
          #        mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2]))
          #    }
          #}else{
          mle.pars.mat <- c()
          for(row.index in 1:dim(mle.pars.mat.red)[1]){
            mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr))
          }
          #}
          print(results.final$objective)
          print(mle.pars.mat)
          lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
          current.likelihood <- results.final$objective
          cat(paste("       Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
          iteration.number <- iteration.number + 1
        }
        #Output for use in sims#
        if(output.by.restart == TRUE){
          obj.tmp = list(np=max(index.matrix) + length(phy$edge.length) + sum(nsites.vector), loglik = -results.final$objective, AIC = -2*(-results.final$objective)+2*(max(index.matrix) + length(phy$edge.length) + sum(nsites.vector)), mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=aa.optim.full.list, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=ip.vector)
          class(obj.tmp) = "selac"
          save(obj.tmp, file=paste(paste(codon.data.path, output.restart.filename, sep=""), number.of.current.restarts, "Rsave", sep="."))
        }
        ########################
        if(results.final$objective < best.lik){
          best.ip <- ip.vector
          best.lik <- results.final$objective
          best.solution <- mle.pars.mat
          best.edge.lengths <- phy$edge.length
          best.aa.optim.list <- aa.optim.full.list
          best.codon.freq.by.aa <- codon.freq.by.aa.list
          best.codon.freq.by.gene <- codon.freq.by.gene.list
        }
        number.of.current.restarts <- number.of.current.restarts + 1
        ip.vector[c(index.matrix[,1])] <- selac.starting.vals[number.of.current.restarts, 1]
        ip.vector[2:3] <- selac.starting.vals[number.of.current.restarts, 2:3]
        aa.optim.list <- aa.optim.original
      }
      selac.starting.vals <- best.ip
      loglik <- -(best.lik) #to go from neglnl to lnl
      mle.pars.mat <- best.solution
      aa.optim.full.list <- best.aa.optim.list
      codon.freq.by.aa.list <- best.codon.freq.by.aa
      codon.freq.by.gene.list <- best.codon.freq.by.gene

      if(edge.length == "optimize"){
        phy$edge.length <- best.edge.lengths
      }
    }else{
      if(optimal.aa == "averaged"){
        aa.optim.list = NULL
      }
      number.of.current.restarts <- 1
      best.lik <- 1000000
      while(number.of.current.restarts < (max.restarts+1)){
        if(optimal.aa == "user"){
          cat(paste("Finished. Performing random restart ", number.of.current.restarts," using user-supplied optimal amino acids...", sep=""), "\n")
        }else{
          cat(paste("Finished. Performing random restart ", number.of.current.restarts," using majority-rule optimal amino acids...", sep=""), "\n")
        }
        mle.pars.mat <- index.matrix
        mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
        cat("       Doing first pass...", "\n")
        print(mle.pars.mat)
        if(edge.length == "optimize"){
          cat("              Optimizing edge lengths", "\n")
          phy$edge.length <- colMeans(starting.branch.lengths)
          #phy$edge.length <- colMeans(starting.branch.lengths) / (1/selac.starting.vals[number.of.current.restarts, 2])
          #opts.edge <- opts
          upper.edge <- rep(log(50), length(phy$edge.length))
          lower.edge <- rep(log(1e-8), length(phy$edge.length))
          results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=NULL, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          print(results.edge.final$objective)
          print(exp(results.edge.final$solution))
          phy$edge.length <- exp(results.edge.final$solution)
        }
        cat("              Optimizing model parameters", "\n")

        #if(include.gamma == TRUE){
        #    alpha.beta.gtr <- mle.pars.mat[1,c(2:(max.par.model.count-1))]
        #    upper.bounds.shared <- upper[c(2:(max.par.model.count-1))]
        #    lower.bounds.shared <- lower[c(2:(max.par.model.count-1))]
        #}else{
        alpha.beta.gtr <- mle.pars.mat[1,c(2:max.par.model.count)]
        upper.bounds.shared <- upper[c(2:max.par.model.count)]
        lower.bounds.shared <- lower[c(2:max.par.model.count)]
        #}

        ParallelizedOptimizedByGene <- function(n.partition){
          #if(include.gamma == TRUE){
          #    tmp.par.mat <- mle.pars.mat[,c(1, max.par.model.count)]
          #    upper.bounds.gene <- upper[c(1,max.par.model.count)]
          #    lower.bounds.gene <- lower[c(1,max.par.model.count)]
          #}else{
          tmp.par.mat <- as.matrix(mle.pars.mat[,1])
          upper.bounds.gene <- upper[1]
          lower.bounds.gene <- lower[1]
          #}
          optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f = OptimizeModelParsAlphaBetaGtrFixed, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, alpha.beta.gtr=alpha.beta.gtr, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, aa.optim_array=aa.optim.list[[n.partition]], codon.freq.by.aa=codon.freq.by.aa.list[[n.partition]], codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
          return(tmp.pars)
        }
        results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
        #if(include.gamma == TRUE){
        #The number of columns is 3: [1] log-likelihood, [2] C.q.phi.Ne, [3] phi gamma:
        #    parallelized.parameters <- t(matrix(unlist(results.set), 3, n.partitions))
        #}else{
        #The number of columns is 2: [1] log-likelihood, [2] C.q.phi.Ne:
        parallelized.parameters <- t(matrix(unlist(results.set), 2, n.partitions))
        #}
        results.final <- NULL
        results.final$objective <- sum(parallelized.parameters[,1])
        results.final$solution <- c(t(parallelized.parameters[,-1]))
        mle.pars.mat.red <- index.matrix.red
        mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]
        print(mle.pars.mat.red)
        optim.alpha.beta.gtr.all.genes <- nloptr(x0=log(alpha.beta.gtr), eval_f = OptimizeAlphaBetaGtrOnly, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix.red, phy=phy, aa.optim_array=aa.optim.list, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
        results.final$objective <- optim.alpha.beta.gtr.all.genes$objective
        alpha.beta.gtr <- exp(optim.alpha.beta.gtr.all.genes$solution)
        #if(include.gamma == TRUE){
        #    mle.pars.mat <- c()
        #    for(row.index in 1:dim(mle.pars.mat.red)[1]){
        #        mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2]))
        #    }
        #}else{
        mle.pars.mat <- c()
        for(row.index in 1:dim(mle.pars.mat.red)[1]){
          mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr))
        }
        #}
        print(results.final$objective)
        print(mle.pars.mat)

        current.likelihood <- results.final$objective
        cat(paste("       Current likelihood", current.likelihood, sep=" "), "\n")
        lik.diff <- 10
        iteration.number <- 1
        while(lik.diff != 0 & iteration.number <= max.iterations){
          cat(paste("       Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
          if(edge.length == "optimize"){
            cat("              Optimizing edge lengths", "\n")
            #opts.edge <- opts
            opts.edge$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))

            results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=aa.optim.list, root.p_array=NULL, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            print(results.edge.final$objective)
            print(exp(results.edge.final$solution))
            phy$edge.length <- exp(results.edge.final$solution)
          }
          cat("              Optimizing model parameters", "\n")

          ParallelizedOptimizedByGene <- function(n.partition){
            #if(include.gamma == TRUE){
            #    tmp.par.mat <- mle.pars.mat[,c(1, max.par.model.count)]
            #    upper.bounds.gene <- upper[c(1,max.par.model.count)]
            #    lower.bounds.gene <- lower[c(1,max.par.model.count)]
            #}else{
            tmp.par.mat <- as.matrix(mle.pars.mat[,1])
            upper.bounds.gene <- upper[1]
            lower.bounds.gene <- lower[1]
            #}
            opts.params <- opts
            opts.params$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
            optim.by.gene <- nloptr(x0=log(tmp.par.mat[n.partition,]), eval_f=OptimizeModelParsAlphaBetaGtrFixed, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts.params, alpha.beta.gtr=alpha.beta.gtr, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, aa.optim_array=aa.optim.list[[n.partition]], codon.freq.by.aa=codon.freq.by.aa.list[[n.partition]], codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
            tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
            return(tmp.pars)
          }
          results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
          #if(include.gamma == TRUE){
          #The number of columns is 3: [1] log-likelihood, [2] C.q.phi, [3] phi gamma:
          #    parallelized.parameters <- t(matrix(unlist(results.set), 3, n.partitions))
          #}else{
          #The number of columns is 2: [1] log-likelihood, [2] C.q.phi:
          parallelized.parameters <- t(matrix(unlist(results.set), 2, n.partitions))
          #}
          results.final <- NULL
          results.final$objective <- sum(parallelized.parameters[,1])
          results.final$solution <- c(t(parallelized.parameters[,-1]))
          mle.pars.mat.red <- index.matrix.red
          mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]

          optim.alpha.beta.gtr.all.genes <- nloptr(x0=log(alpha.beta.gtr), eval_f = OptimizeAlphaBetaGtrOnly, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix.red, phy=phy, aa.optim_array=aa.optim.list, codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=FALSE, neglnl=TRUE, HMM=FALSE)
          results.final$objective <- optim.alpha.beta.gtr.all.genes$objective
          alpha.beta.gtr <- exp(optim.alpha.beta.gtr.all.genes$solution)
          #if(include.gamma == TRUE){
          #   mle.pars.mat <- c()
          #   for(row.index in 1:dim(mle.pars.mat.red)[1]){
          #        mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2]))
          #    }
          #}else{
          mle.pars.mat <- c()
          for(row.index in 1:dim(mle.pars.mat.red)[1]){
            mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr))
          }
          #}
          print(results.final$objective)
          print(mle.pars.mat)
          lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
          current.likelihood <- results.final$objective
          cat(paste("       Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
          iteration.number <- iteration.number + 1
        }
        #Output for use in sims#
        if(output.by.restart == TRUE){
          obj.tmp = list(np=max(index.matrix) + length(phy$edge.length) + sum(nsites.vector), loglik = -results.final$objective, AIC = -2*(-results.final$objective)+2*(max(index.matrix) + length(phy$edge.length) + sum(nsites.vector)), mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=aa.optim.full.list, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=ip.vector)
          class(obj.tmp) = "selac"
          save(obj.tmp, file=paste(paste(codon.data.path, output.restart.filename, sep=""), number.of.current.restarts, "Rsave", sep="."))
        }
        ########################
        if(results.final$objective < best.lik){
          best.ip <- ip.vector
          best.lik <- results.final$objective
          best.solution <- mle.pars.mat
          best.edge.lengths <- phy$edge.length
          best.codon.freq.by.aa <- codon.freq.by.aa.list
          best.codon.freq.by.gene <- codon.freq.by.gene.list
        }
        number.of.current.restarts <- number.of.current.restarts + 1
        print(ip.vector)
        ip.vector[c(index.matrix[,1])] <- selac.starting.vals[number.of.current.restarts, 1]
        ip.vector[2:3] <- selac.starting.vals[number.of.current.restarts, 2:3]
        print(ip.vector)
      }
      selac.starting.vals <- best.ip
      loglik <- -(best.lik) #to go from neglnl to lnl
      mle.pars.mat <- best.solution
      codon.freq.by.aa.list <- best.codon.freq.by.aa
      codon.freq.by.gene.list <- best.codon.freq.by.gene

      if(edge.length == "optimize"){
        phy$edge.length <- best.edge.lengths
      }
    }
    cat("Finished. Summarizing results...", "\n")
    colnames(mle.pars.mat) <- parameter.column.names

    if(edge.length == "optimize"){
      if(optimal.aa == "user" | optimal.aa == "majrule" | optimal.aa == "averaged"){
        np <- max(index.matrix) + length(phy$edge.length)
      }else{
        np <- max(index.matrix) + length(phy$edge.length) + sum(nsites.vector)
      }
    }else{
      if(optimal.aa == "user" | optimal.aa == "majrule" | optimal.aa == "averaged"){
        np <- max(index.matrix)
      }else{
        np <- max(index.matrix) + sum(nsites.vector)
      }
    }

    #Counting parameters: Do we count the nsites too? Yup.
    obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=aa.optim.full.list, aa.optim.type=optimal.aa, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], codon.freq.by.aa=codon.freq.by.aa.list, codon.freq.by.gene=codon.freq.by.gene.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=selac.starting.vals)
    class(obj) = "selac"
  }
  return(obj)
}


#' @title Efficient optimization of a Hidden Markov SELAC model
#'
#' @description
#' Efficient optimization of model parameters under a HMM SELAC model
#'
#' @param codon.data.path Provides the path to the directory containing the gene specific fasta files of coding data. Must have a ".fasta" line ending.
#' @param n.partitions The number of partitions to analyze. The order is based on the Unix order of the fasta files in the directory.
#' @param phy The phylogenetic tree to optimize the model parameters.
#' @param data.type The data type being tested. Options are "codon" or "nucleotide".
#' @param codon.model The type of codon model to use. There are four options: "none", "GY94", "FMutSel0", "selac".
#' @param edge.length Indicates whether or not edge lengths should be optimized. By default it is set to "optimize", other option is "fixed", which user-supplied branch lengths.
#' @param edge.linked A logical indicating whether or not edge lengths should be optimized separately for each gene. By default, a single set of each lengths is optimized for all genes.
#' @param optimal.aa Indicates what type of optimal.aa should be used. There are four options: "none", "majrule", "optimize", or "user".
#' @param nuc.model Indicates what type nucleotide model to use. There are three options: "JC", "GTR", or "UNREST".
#' @param estimate.aa.importance Indicates whether gene specific importance of distance parameter is to be estimate.
#' @param include.gamma A logical indicating whether or not to include a discrete gamma model.
#' @param gamma.type Indicates what type of gamma distribution to use. Options are "quadrature" after the Laguerre quadrature approach of Felsenstein 2001 or median approach of Yang 1994.
#' @param ncats The number of discrete categories.
#' @param numcode The ncbi genetic code number for translation. By default the standard (numcode=1) genetic code is used.
#' @param diploid A logical indicating whether or not the organism is diploid or not.
#' @param k.levels Provides how many levels in the polynomial. By default we assume a single level (i.e., linear).
#' @param aa.properties User-supplied amino acid distance properties. By default we assume Grantham (1974) properties.
#' @param verbose Logical indicating whether each iteration be printed to the screen.
#' @param n.cores.by.gene The number of cores to dedicate to parallelize analyses across gene.
#' @param n.cores.by.gene.by.site The number of cores to decidate to parallelize analyses by site WITHIN a gene. Note n.cores.by.gene*n.cores.by.gene.by.site is the total number of cores dedicated to the analysis.
#' @param max.tol Supplies the relative optimization tolerance.
#' @param max.tol.edges Supplies the relative optimization tolerance for branch lengths only. Default is that is the same as the max.tol.
#' @param max.evals Supplies the max number of iterations tried during optimization.
#' @param max.restarts Supplies the number of random restarts.
#' @param user.optimal.aa If optimal.aa is set to "user", this option allows for the user-input optimal amino acids. Must be a list. To get the proper order of the partitions see "GetPartitionOrder" documentation.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.
#' @param recalculate.starting.brlen Whether to use given branch lengths in the starting tree or recalculate them.
#' @param output.by.restart Logical indicating whether or not each restart is saved to a file. Default is TRUE.
#' @param output.restart.filename Designates the file name for each random restart.
#' @param user.supplied.starting.param.vals Designates user-supplied starting values for C.q.phi.Ne, Grantham alpha, and Grantham beta. Default is NULL.
#' @param tol.step If > 1, makes for coarser tolerance at earlier iterations of the optimizer
#' @param optimizer.algorithm The optimizer used by nloptr.
#'
#' @details
#' A hidden Markov model which no longers optimizes the optimal amino acids, but instead allows for the optimal sequence to vary along branches, clades, taxa, etc. Like the original function, we optimize parameters across each gene separately while keeping the shared parameters, alpha, beta, edge lengths, and nucleotide substitution parameters constant across genes. We then optimize alpha, beta, gtr, and the edge lengths while keeping the rest of the parameters for each gene fixed. This approach is potentially more efficient than simply optimizing all parameters simultaneously, especially if fitting models across 100's of genes.
SelacHMMOptimize <- function(codon.data.path, n.partitions=NULL, phy, data.type="codon", codon.model="selac", edge.length="optimize", edge.linked=TRUE, nuc.model="GTR", estimate.aa.importance=FALSE, include.gamma=FALSE, gamma.type="quadrature", ncats=4, numcode=1, diploid=TRUE, k.levels=0, aa.properties=NULL, verbose=FALSE, n.cores.by.gene=1, n.cores.by.gene.by.site=1, max.tol=1e-3, max.tol.edges=1e-3, max.evals=1000000, max.restarts=3, user.optimal.aa=NULL, fasta.rows.to.keep=NULL, recalculate.starting.brlen=TRUE, output.by.restart=TRUE, output.restart.filename="restartResult", user.supplied.starting.param.vals=NULL, tol.step=1, optimizer.algorithm="NLOPT_LN_SBPLX", max.iterations=6) {

  if(!data.type == "codon"){
    stop("Check that your data type input is correct. Options currently are codon only", call.=FALSE)
  }
  if(!codon.model == "selac"){
    stop("Check that your codon model is correct. Options are currently selac only", call.=FALSE)
  }
  if(!edge.length == "optimize" & !edge.length == "fixed"){
    stop("Check that you have a supported edge length option. Options are optimize or fixed.", call.=FALSE)
  }
  if(!nuc.model == "JC" & !nuc.model == "GTR" & !nuc.model == "UNREST"){
    stop("Check that you have a supported nucleotide substitution model. Options are JC, GTR, or UNREST.", call.=FALSE)
  }
  if(!gamma.type == "quadrature" & !gamma.type == "median"){
    stop("Check that you have a supported gamma type. Options are quadrature after Felsenstein 2001 or median after Yang 1994.", call.=FALSE)
  }

  cat(paste("Using", n.cores.by.gene * n.cores.by.gene.by.site, "total processors", sep=" "), "\n")

  cat("Initializing data and model parameters...", "\n")

  partitions <- system(paste("ls -1 ", codon.data.path, "*.fasta", sep=""), intern=TRUE)

  estimate.importance <- estimate.aa.importance

  if(is.null(n.partitions)){
    n.partitions <- length(partitions)
  }else{
    n.partitions = n.partitions
  }

  site.pattern.data.list <- as.list(numeric(n.partitions))
  site.pattern.count.list <- as.list(numeric(n.partitions))
  nsites.vector <- c()
  codon.freq.by.aa.list <- as.list(numeric(n.partitions))
  codon.freq.by.gene.list <- as.list(numeric(n.partitions))
  starting.branch.lengths <- matrix(0, n.partitions, length(phy$edge[,1]))
  for (partition.index in sequence(n.partitions)) {
    gene.tmp <- read.dna(partitions[partition.index], format='fasta')
    if(!is.null(fasta.rows.to.keep)){
      gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
    }else{
      gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
    }
    starting.branch.lengths[partition.index,] <- ComputeStartingBranchLengths(phy, gene.tmp, data.type=data.type, recalculate.starting.brlen=recalculate.starting.brlen)$edge.length
    codon.data <- DNAbinToCodonNumeric(gene.tmp)
    codon.data <- codon.data[phy$tip.label,]
    nsites.vector = c(nsites.vector, dim(codon.data)[2] - 1)
    codon.freq.by.aa.list[[partition.index]] <- NULL
    codon.freq.by.gene.list[[partition.index]] <- GetCodonFreqsByGeneHMM(codon.data[,-1])
    codon.data <- SitePattern(codon.data, includes.optimal.aa=FALSE)
    site.pattern.data.list[[partition.index]] = codon.data$unique.site.patterns
    site.pattern.count.list[[partition.index]] = codon.data$site.pattern.counts
  }

  opts <- list("algorithm" = optimizer.algorithm, "maxeval" = max.evals, "ftol_rel" = max.tol)
  opts.edge <- list("algorithm" = optimizer.algorithm, "maxeval" = max.evals, "ftol_rel" = max.tol.edges)

  results.final <- c()
  if(nuc.model == "JC"){
    nuc.ip = NULL
    max.par.model.count = 0
    parameter.column.names <- c("sel.reg")
  }
  if(nuc.model == "GTR"){
    nuc.ip = rep(1, 5)
    max.par.model.count = 5
    parameter.column.names <- c("C_A", "G_A", "T_A", "G_C", "T_C", "sel.reg")
  }
  if(nuc.model == "UNREST"){
    nuc.ip = rep(1, 11)
    max.par.model.count = 11
    parameter.column.names <- c("C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "sel.reg")
  }

  codon.index.matrix <- CreateCodonMutationMatrixIndexEvolveAA()
  cpv.starting.parameters <- GetAADistanceStartingParameters(aa.properties=aa.properties)
  if(max.restarts > 1){
    selac.starting.vals <- matrix(0, max.restarts+1, 3)
    selac.starting.vals[,1] <- runif(n = max.restarts+1, min = (10^-10)*5e6, max = (10^-5)*5e6)
    selac.starting.vals[,2] <- runif(n = max.restarts+1, min = 0.01, max = 3)
    selac.starting.vals[,3] <- runif(n = max.restarts+1, min = 0.01, max = 1)
  }else{
    if(is.null(user.supplied.starting.param.vals)){
      selac.starting.vals <- matrix(c(2, 1.8292716544, 0.1017990371), 1, 3)
      selac.starting.vals <- rbind(selac.starting.vals, c(2, 1.8292716544, 0.1017990371))
    }else{
      selac.starting.vals <- matrix(c(user.supplied.starting.param.vals[1], user.supplied.starting.param.vals[2], user.supplied.starting.param.vals[3]), 1, 3)
      selac.starting.vals <- rbind(selac.starting.vals, c(user.supplied.starting.param.vals[1], user.supplied.starting.param.vals[2], user.supplied.starting.param.vals[3]))
    }
  }
  if(include.gamma == TRUE){
    #Gamma variation is turned ON:
    if(nuc.model == "JC"){
      if(k.levels == 0){
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "shape.gamma", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 5, 21, 4)
          lower = c(rep(-21, length(ip)-1), -4)
          max.par.model.count = 6 + 0 + 1 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "shape.gamma", "sel.reg")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 5, 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 1 + 1
        }
      }else{
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, 1, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "shape.gamma", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 21, 21, 5, 21, 4)
          lower = c(rep(-21, length(ip)-1), -4)
          max.par.model.count = 6 + 0 + 2 + 1 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, 1, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "shape.gamma", "sel.reg")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 21, 21, 5, 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 2 + 1 + 1
        }
      }
    }
    if(nuc.model == "GTR") {
      if(k.levels == 0){
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip, 1, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "C_A", "G_A", "T_A", "G_C", "T_C", "shape.gamma", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), 0, 0, 0, rep(21, length(nuc.ip)), 5, 21, 4)
          lower = c(rep(-21, length(ip)-1), -4)
          max.par.model.count = 6 + 5 + 1 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip, 1, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "C_A", "G_A", "T_A", "G_C", "T_C", "shape.gamma", "sel.reg")
          upper = c(log(50), log(20), log(2), 0, 0, 0, rep(21, length(nuc.ip)), 5, 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 1 + 1
        }
      }else{
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, nuc.ip, 1, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "C_A", "G_A", "T_A", "G_C", "T_C", "shape.gamma", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 21, 21, rep(21, length(nuc.ip)), 5, 21, 4)
          lower = c(rep(-21, length(ip)-1), -4)
          max.par.model.count = 6 + 5 + 2	+ 1 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, nuc.ip, 1, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "C_A", "G_A", "T_A", "G_C", "T_C", "shape.gamma", "sel.reg")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 21, 21, rep(21, length(nuc.ip)), 5, 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 2	+ 1 + 1
        }
      }
    }
    if(nuc.model == "UNREST") {
      if(k.levels == 0){
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], nuc.ip, 1, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "shape.gamma", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), rep(21, length(nuc.ip)), 5, 21, 4)
          lower = c(rep(-21, length(ip)-1), -4)
          max.par.model.count = 3 + 11 + 1 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], nuc.ip, 1, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "shape.gamma", "sel.reg")
          upper = c(log(50), log(20), log(20), rep(21, length(nuc.ip)), 5, 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 1 + 1
        }
      }else{
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], 1, 1, nuc.ip, 1, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "a0", "a1", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "shape.gamma", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), 21, 21, rep(21, length(nuc.ip)), 5, 21, 4)
          lower = c(rep(-21, length(ip)-1), -4)
          max.par.model.count = 3 + 11 + 2 + 1 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], 1, 1, nuc.ip, 1, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "a0", "a1", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "shape.gamma", "sel.reg")
          upper = c(log(50), log(20), log(20), 21, 21, rep(21, length(nuc.ip)), 5, 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 2 + 1 + 1
        }
      }
    }
    index.matrix = matrix(0, n.partitions, max.par.model.count)
    index.matrix[1,] = 1:ncol(index.matrix)
    ip.vector = ip
    if(n.partitions > 1){
      # Gamma variation is turned ON:
      for(partition.index in 2:n.partitions){
        if(nuc.model == "JC"){
          if(estimate.importance == TRUE){
            ip.vector = c(ip.vector, ip[1], ip[8], ip[9])
          }else{
            ip.vector = c(ip.vector, ip[1], ip[8])
          }
        }else{
          if(nuc.model == "GTR"){
            index.matrix.tmp = numeric(max.par.model.count)
            if(k.levels == 0){
              if(estimate.importance == TRUE){
                index.matrix.tmp[c(2:12)] = c(2:12)
                ip.vector = c(ip.vector, ip[1], ip[13], ip[14])
              }else{
                index.matrix.tmp[c(2:12)] = c(2:12)
                ip.vector = c(ip.vector, ip[1], ip[13])
              }
            }else{
              if(estimate.importance == TRUE){
                index.matrix.tmp[c(2:14)] = c(2:14)
                ip.vector = c(ip.vector, ip[1], ip[15], ip[16])
              }else{
                index.matrix.tmp[c(2:14)] = c(2:14)
                ip.vector = c(ip.vector, ip[1], ip[15])
              }
            }
          }else{
            index.matrix.tmp = numeric(max.par.model.count)
            if(k.levels == 0){
              if(estimate.importance == TRUE){
                index.matrix.tmp[c(2:15)] = c(2:15)
                ip.vector = c(ip.vector, ip[1], ip[16], ip[17])
              }else{
                index.matrix.tmp[c(2:15)] = c(2:15)
                ip.vector = c(ip.vector, ip[1], ip[16])
              }
            }else{
              if(estimate.importance == TRUE){
                index.matrix.tmp[c(2:17)] = c(2:17)
                ip.vector = c(ip.vector, ip[1], ip[18], ip[19])
              }else{
                index.matrix.tmp[c(2:17)] = c(2:17)
                ip.vector = c(ip.vector, ip[1], ip[18])
              }
            }
          }
        }
        index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
        index.matrix[partition.index,] <- index.matrix.tmp
      }
    }
  }else{
    # Gamma variation is turned OFF:
    if(nuc.model == "JC"){
      if(k.levels == 0){
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 21, 4)
          lower = c(rep(-21, length(ip)-1), -4)
          max.par.model.count = 6 + 0 + 0 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "sel.reg")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 0 + 1
        }
      }else{
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 21, 21, 21, 4)
          lower = c(rep(-21, length(ip)-1), -4)
          max.par.model.count = 6 + 0 + 0 + 2 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "sel.reg")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 21, 21, 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 0 + 0 + 2 + 1
        }
      }
    }
    if(nuc.model == "GTR") {
      if(k.levels == 0){
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "C_A", "G_A", "T_A", "G_C", "T_C", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), 0, 0, 0, rep(21, length(nuc.ip)), 21, 4)
          lower = c(rep(-21, length(ip)-1), -4)
          max.par.model.count = 6 + 5 + 0 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, nuc.ip, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "C_A", "G_A", "T_A", "G_C", "T_C", "sel.reg")
          upper = c(log(50), log(20), log(20), 0, 0, 0, rep(21, length(nuc.ip)), 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 0 + 1
        }
      }else{
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, nuc.ip, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "C_A", "G_A", "T_A", "G_C", "T_C", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 21, 21, rep(21, length(nuc.ip)), 21, 4)
          lower = c(rep(-21, length(ip)-1), -4)
          max.par.model.count = 6 + 5 + 0 + 2 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], 0.25, 0.25, 0.25, 1, 1, nuc.ip, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "freqA", "freqC", "freqG", "a0", "a1", "C_A", "G_A", "T_A", "G_C", "T_C", "sel.reg")
          upper = c(log(50), log(20), log(20), 0, 0, 0, 21, 21, rep(21, length(nuc.ip)), 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 6 + 5 + 0 + 2 + 1
        }
      }
    }
    if(nuc.model == "UNREST") {
      if(k.levels == 0){
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], nuc.ip, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), rep(21, length(nuc.ip)), 21, 4)
          lower = c(rep(-21, length(ip)-1), -4)
          max.par.model.count = 3 + 11 + 0 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], nuc.ip, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "sel.reg")
          upper = c(log(50), log(20), log(20), rep(21, length(nuc.ip)), 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 0 + 1
        }
      }else{
        if(estimate.importance == TRUE){
          ip = c(selac.starting.vals[1,1:3], 1, 1, nuc.ip, 0.01, 1)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "a0", "a1", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "sel.reg", "importance")
          upper = c(log(50), log(20), log(20), 21, 21, rep(21, length(nuc.ip)), 21, 4)
          lower = rep(-21, length(ip)-1, -4)
          max.par.model.count = 3 + 11 + 0 + 2 + 1 + 1
        }else{
          ip = c(selac.starting.vals[1,1:3], 1, 1, nuc.ip, 0.01)
          parameter.column.names <- c("C.q.phi.Ne", "alpha", "beta", "a0", "a1", "C_A", "G_A", "T_A", "A_C", "G_C", "T_C", "A_G", "C_G", "A_T", "C_T", "G_T", "sel.reg")
          upper = c(log(50), log(20), log(20), 21, 21, rep(21, length(nuc.ip)), 21)
          lower = rep(-21, length(ip))
          max.par.model.count = 3 + 11 + 0 + 2 + 1
        }
      }
    }
    index.matrix = matrix(0, n.partitions, max.par.model.count)
    index.matrix[1,] = 1:ncol(index.matrix)
    ip.vector = ip
    if(n.partitions > 1){
      for(partition.index in 2:n.partitions){
        if(nuc.model == "JC"){
          if(estimate.importance == TRUE){
            ip.vector = c(ip.vector, ip[1], ip[7], ip[8])
          }else{
            ip.vector = c(ip.vector, ip[1], ip[7])
          }
        }else{
          if(nuc.model == "GTR"){
            index.matrix.tmp = numeric(max.par.model.count)
            if(k.levels == 0){
              if(estimate.importance == TRUE){
                index.matrix.tmp[c(2:11)] = c(2:11)
                ip.vector = c(ip.vector, ip[1], ip[12], ip[13])
              }else{
                index.matrix.tmp[c(2:11)] = c(2:11)
                ip.vector = c(ip.vector, ip[1], ip[12])
              }
            }else{
              if(estimate.importance == TRUE){
                index.matrix.tmp[c(2:13)] = c(2:13)
                ip.vector = c(ip.vector, ip[1], ip[14], ip[15])
              }else{
                index.matrix.tmp[c(2:13)] = c(2:13)
                ip.vector = c(ip.vector, ip[1], ip[14])
              }
            }
          }else{
            index.matrix.tmp = numeric(max.par.model.count)
            if(k.levels == 0){
              if(estimate.importance == TRUE){
                index.matrix.tmp[c(2:14)] = c(2:14)
                ip.vector = c(ip.vector, ip[1], ip[15], ip[16])
              }else{
                index.matrix.tmp[c(2:14)] = c(2:14)
                ip.vector = c(ip.vector, ip[1], ip[15])
              }
            }else{
              if(estimate.importance == TRUE){
                index.matrix.tmp[c(2:16)] = c(2:16)
                ip.vector = c(ip.vector, ip[1], ip[17], ip[18])
              }else{
                index.matrix.tmp[c(2:16)] = c(2:16)
                ip.vector = c(ip.vector, ip[1], ip[17])
              }
            }
          }
        }
        index.matrix.tmp[index.matrix.tmp==0] = seq(max(index.matrix)+1, length.out=length(index.matrix.tmp[index.matrix.tmp==0]))
        index.matrix[partition.index,] <- index.matrix.tmp
      }
    }
  }

  #THIS IS FOR THERE IS A SEPARATE GAMMA PER GENE:
  if(estimate.importance == TRUE){
    index.matrix.red <- t(matrix(1:(n.partitions*3), 3, n.partitions))
  }else{
    index.matrix.red <- t(matrix(1:(n.partitions*2), 2, n.partitions))
  }
  #This is so we can break out alpha, beta, GTR, and gamma which are shared among ALL genes:
  #index.matrix.red <- t(matrix(1:n.partitions, 1, n.partitions))

  number.of.current.restarts <- 1
  best.lik <- 1000000
  while(number.of.current.restarts < (max.restarts+1)){
    cat(paste("Finished. Performing random restart ", number.of.current.restarts," ...", sep=""), "\n")
    mle.pars.mat <- index.matrix
    mle.pars.mat[] <- c(ip.vector, 0)[index.matrix]
    cat("       Doing first pass...", "\n")
    print(mle.pars.mat)
    if(edge.length == "optimize"){
      cat("              Optimizing edge lengths", "\n")
      phy$edge.length <- colMeans(starting.branch.lengths)
      #phy$edge.length <- colMeans(starting.branch.lengths) / (1/selac.starting.vals[number.of.current.restarts, 2])
      #opts.edge <- opts
      upper.edge <- rep(log(10), length(phy$edge.length))
      lower.edge <- rep(log(1e-8), length(phy$edge.length))
      results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=estimate.importance, neglnl=TRUE, HMM=TRUE)
      print("here!!")
      print(results.edge.final$objective)
      print(exp(results.edge.final$solution))
      phy$edge.length <- exp(results.edge.final$solution)
    }
    cat("              Optimizing model parameters", "\n")

    if(estimate.importance == TRUE){
      alpha.beta.gtr <- mle.pars.mat[1,c(2:(max.par.model.count-2))]
      upper.bounds.shared <- upper[c(2:(max.par.model.count-2))]
      lower.bounds.shared <- lower[c(2:(max.par.model.count-2))]
    }else{
      alpha.beta.gtr <- mle.pars.mat[1,c(2:(max.par.model.count-1))]
      upper.bounds.shared <- upper[c(2:(max.par.model.count-1))]
      lower.bounds.shared <- lower[c(2:(max.par.model.count-1))]
    }

    ParallelizedOptimizedByGene <- function(n.partition){
      if(estimate.importance == TRUE){
        # make sure this is always a vector no matter the number of partitions
        tmp.par.mat <- mle.pars.mat[n.partition,c(1, max.par.model.count-1, max.par.model.count)]
        upper.bounds.gene <- upper[c(1, max.par.model.count-1, max.par.model.count)]
        lower.bounds.gene <- lower[c(1, max.par.model.count-1, max.par.model.count)]
      }else{
        tmp.par.mat <- mle.pars.mat[,c(1, max.par.model.count)]
        upper.bounds.gene <- upper[c(1, max.par.model.count)]
        lower.bounds.gene <- lower[c(1, max.par.model.count)]
      }
      optim.by.gene <- nloptr(x0=log(tmp.par.mat), eval_f = OptimizeModelParsAlphaBetaGtrFixed, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, alpha.beta.gtr=alpha.beta.gtr, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=estimate.importance, neglnl=TRUE, HMM=TRUE)
      tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
      return(tmp.pars)
    }
    results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
    if(estimate.importance == TRUE){
      #The number of columns is 3: [1] log-likelihood, [2] C.q.phi.Ne, [3] aa transition [4] importance:
      parallelized.parameters <- t(matrix(unlist(results.set), 4, n.partitions))
    }else{
      #The number of columns is 2: [1] log-likelihood, [2] C.q.phi.Ne [3] aa transition:
      parallelized.parameters <- t(matrix(unlist(results.set), 3, n.partitions))
    }

    results.final <- NULL
    results.final$objective <- sum(parallelized.parameters[,1])
    results.final$solution <- c(t(parallelized.parameters[,-1]))
    mle.pars.mat.red <- index.matrix.red
    mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]
    print(mle.pars.mat.red)
    optim.alpha.beta.gtr.all.genes <- nloptr(x0=log(alpha.beta.gtr), eval_f = OptimizeAlphaBetaGtrOnly, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix.red, phy=phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=estimate.importance, neglnl=TRUE, HMM=TRUE)
    results.final$objective <- optim.alpha.beta.gtr.all.genes$objective
    alpha.beta.gtr <- exp(optim.alpha.beta.gtr.all.genes$solution)
    if(estimate.importance == TRUE){
      mle.pars.mat <- c()
      for(row.index in 1:dim(mle.pars.mat.red)[1]){
        mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2], mle.pars.mat.red[row.index,3]))
      }
    }else{
      mle.pars.mat <- c()
      for(row.index in 1:dim(mle.pars.mat.red)[1]){
        mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2]))
      }
    }
    print(results.final$objective)
    print(mle.pars.mat)

    current.likelihood <- results.final$objective
    cat(paste("       Current likelihood", current.likelihood, sep=" "), "\n")
    lik.diff <- 10
    iteration.number <- 1
    while(lik.diff != 0 & iteration.number <= max.iterations){
      cat(paste("       Finished. Iterating search -- Round", iteration.number, sep=" "), "\n")
      if(edge.length == "optimize"){
        cat("              Optimizing edge lengths", "\n")
        #opts.edge <- opts
        opts.edge$ftol_rel <- opts$ftol_rel * (max(1,tol.step^((max.iterations+1)-iteration.number)))
        results.edge.final <- nloptr(x0=log(phy$edge.length), eval_f = OptimizeEdgeLengths, ub=upper.edge, lb=lower.edge, opts=opts.edge, par.mat=mle.pars.mat, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix, phy=phy, aa.optim_array=NULL, root.p_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=estimate.importance, neglnl=TRUE, HMM=TRUE)
        print(results.edge.final$objective)
        print(exp(results.edge.final$solution))
        phy$edge.length <- exp(results.edge.final$solution)
      }
      cat("              Optimizing model parameters", "\n")

      ParallelizedOptimizedByGene <- function(n.partition){
        if(estimate.importance == TRUE){
          tmp.par.mat <- mle.pars.mat[n.partition,c(1, max.par.model.count-1, max.par.model.count)]
          upper.bounds.gene <- upper[c(1, max.par.model.count-1, max.par.model.count)]
          lower.bounds.gene <- lower[c(1, max.par.model.count-1, max.par.model.count)]
        }else{
          tmp.par.mat <- mle.pars.mat[,c(1, max.par.model.count)]
          upper.bounds.gene <- upper[c(1, max.par.model.count)]
          lower.bounds.gene <- lower[c(1, max.par.model.count)]
        }
        optim.by.gene <- nloptr(x0=log(tmp.par.mat), eval_f = OptimizeModelParsAlphaBetaGtrFixed, ub=upper.bounds.gene, lb=lower.bounds.gene, opts=opts, alpha.beta.gtr=alpha.beta.gtr, codon.site.data=site.pattern.data.list[[n.partition]], codon.site.counts=site.pattern.count.list[[n.partition]], data.type=data.type, codon.model=codon.model, n.partitions=1, nsites.vector=nsites.vector[n.partition], index.matrix=index.matrix.red[1,], phy=phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=codon.freq.by.gene.list[[n.partition]], numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=estimate.importance, neglnl=TRUE, HMM=TRUE)
        tmp.pars <- c(optim.by.gene$objective, optim.by.gene$solution)
        return(tmp.pars)
      }
      results.set <- mclapply(1:n.partitions, ParallelizedOptimizedByGene, mc.cores=n.cores.by.gene)
      if(estimate.importance == TRUE){
        #The number of columns is 3: [1] log-likelihood, [2] C.q.phi.Ne, [3] aa transition [4] importance:
        parallelized.parameters <- t(matrix(unlist(results.set), 4, n.partitions))
      }else{
        #The number of columns is 2: [1] log-likelihood, [2] C.q.phi.Ne [3] aa transition:
        parallelized.parameters <- t(matrix(unlist(results.set), 3, n.partitions))
      }

      results.final <- NULL
      results.final$objective <- sum(parallelized.parameters[,1])
      results.final$solution <- c(t(parallelized.parameters[,-1]))
      mle.pars.mat.red <- index.matrix.red
      mle.pars.mat.red[] <- c(exp(results.final$solution), 0)[index.matrix.red]

      print(mle.pars.mat.red)
      optim.alpha.beta.gtr.all.genes <- nloptr(x0=log(alpha.beta.gtr), eval_f = OptimizeAlphaBetaGtrOnly, ub=upper.bounds.shared, lb=lower.bounds.shared, opts=opts, fixed.pars=mle.pars.mat.red, codon.site.data=site.pattern.data.list, codon.site.counts=site.pattern.count.list, data.type=data.type, codon.model=codon.model, n.partitions=n.partitions, nsites.vector=nsites.vector, index.matrix=index.matrix.red, phy=phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=codon.freq.by.gene.list, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, edge.length=edge.length, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=TRUE, verbose=verbose, n.cores.by.gene=n.cores.by.gene, n.cores.by.gene.by.site=n.cores.by.gene.by.site, estimate.importance=estimate.importance, neglnl=TRUE, HMM=TRUE)
      results.final$objective <- optim.alpha.beta.gtr.all.genes$objective
      alpha.beta.gtr <- exp(optim.alpha.beta.gtr.all.genes$solution)
      if(estimate.importance == TRUE){
        mle.pars.mat <- c()
        for(row.index in 1:dim(mle.pars.mat.red)[1]){
          mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2], mle.pars.mat.red[row.index,3]))
        }
      }else{
        mle.pars.mat <- c()
        for(row.index in 1:dim(mle.pars.mat.red)[1]){
          mle.pars.mat <- rbind(mle.pars.mat, c(mle.pars.mat.red[row.index,1], alpha.beta.gtr, mle.pars.mat.red[row.index,2]))
        }
      }
      print(results.final$objective)
      print(mle.pars.mat)
      lik.diff <- round(abs(current.likelihood-results.final$objective), 8)
      current.likelihood <- results.final$objective
      cat(paste("       Current likelihood", current.likelihood, sep=" "), paste("difference from previous round", lik.diff, sep=" "), "\n")
      iteration.number <- iteration.number + 1
    }
    #Output for use in sims#
    if(output.by.restart == TRUE){
      obj.tmp = list(np=max(index.matrix) + length(phy$edge.length) + sum(nsites.vector), loglik = -results.final$objective, AIC = -2*(-results.final$objective)+2*(max(index.matrix) + length(phy$edge.length) + sum(nsites.vector)), mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=NULL, aa.optim.type=NULL, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], codon.freq.by.aa=NULL, codon.freq.by.gene=codon.freq.by.gene.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=ip.vector)
      class(obj.tmp) = "selac"
      save(obj.tmp, file=paste(paste(codon.data.path, output.restart.filename, sep=""), number.of.current.restarts, "Rsave", sep="."))
    }
    ########################
    if(results.final$objective < best.lik){
      best.ip <- ip.vector
      best.lik <- results.final$objective
      best.solution <- mle.pars.mat
      best.edge.lengths <- phy$edge.length
      best.codon.freq.by.aa <- codon.freq.by.aa.list
      best.codon.freq.by.gene <- codon.freq.by.gene.list
    }
    number.of.current.restarts <- number.of.current.restarts + 1
    print(ip.vector)
    ip.vector[c(index.matrix[,1])] <- selac.starting.vals[number.of.current.restarts, 1]
    ip.vector[2:3] <- selac.starting.vals[number.of.current.restarts, 2:3]
    print(ip.vector)
  }
  selac.starting.vals <- best.ip
  loglik <- -(best.lik) #to go from neglnl to lnl
  mle.pars.mat <- best.solution
  codon.freq.by.aa.list <- best.codon.freq.by.aa
  codon.freq.by.gene.list <- best.codon.freq.by.gene

  if(edge.length == "optimize"){
    phy$edge.length <- best.edge.lengths
  }

  cat("Finished. Summarizing results...", "\n")
  colnames(mle.pars.mat) <- parameter.column.names

  if(edge.length == "optimize"){
    np <- max(index.matrix) + length(phy$edge.length)
  }else{
    np <- max(index.matrix)
  }

  #Counting parameters: Do we count the nsites too? Yup.
  obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, mle.pars=mle.pars.mat, index.matrix=index.matrix, partitions=partitions[1:n.partitions], opts=opts, phy=phy, nsites=nsites.vector, data.type=data.type, codon.model=codon.model, aa.optim=NULL, aa.optim.type=NULL, nuc.model=nuc.model, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=cpv.starting.parameters[3], codon.freq.by.aa=NULL, codon.freq.by.gene=codon.freq.by.gene.list, max.tol=max.tol, max.evals=max.evals, selac.starting.vals=selac.starting.vals)
  class(obj) = "selac"
  return(obj)
}



######################################################################################################################################
######################################################################################################################################
### Utility function getting raw likelihoods not just across genes, but across sites. Also allows for different optimal AA to try
######################################################################################################################################
######################################################################################################################################

#' @title Get data partiion order
#'
#' @description
#'  Provides the order of the partitions after the data is read into SELAC.
#'
#' @param codon.data.path Provides the path to the directory containing the gene specific fasta files of coding data. Must have a ".fasta" line ending.
#'
#' @details
#' Provides the order of the partitions when the data is read into SELAC. This function is mainly useful for when users want to supply their own optimal amino acid list into SELAC.
GetPartitionOrder <- function(codon.data.path){
  partitions <- system(paste("ls -1 ", codon.data.path, "*.fasta", sep=""), intern=TRUE)
  return(partitions)
}


#' @title Calculate functionality
#'
#' @description
#' Calculates the functionality of a single gene
#'
#' @param gene.length Indicates the length of the gene used to calculate functionality.
#' @param aa.data A matrix of amino acids
#' @param optimal.aa A vector of inferred optimal amino acids.
#' @param alpha The inferred Grantham composition paramter
#' @param beta The inferred Grantham polarity parameter
#' @param gamma The inferred Grantham molecular volume parameter
#' @param gp A vector of gamma rates for calculating among site hetergeneity in functionality.
#' @param aa.properties User-supplied amino acid distance properties. By default we assume Grantham (1974) properties.
#'
#' @details
#' The purpose of this function is to provide the functionality of a gene based on the inferred parameters from SelAC. The functionality is often used to scale phi.
GetFunctionality <- function(gene.length, aa.data, optimal.aa, alpha, beta, gamma, gp=NULL, aa.properties=NULL){
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
    if(aa.data[,site.index]!="NA"){
      #broke this up to make debugging easier:
      distance <- ((alpha*(aa.properties[aa.data[,site.index],1] - aa.properties[optimal.aa[site.index],1])^2 + beta*(aa.properties[aa.data[,site.index],2]-aa.properties[optimal.aa[site.index],2])^2+gamma*(aa.properties[aa.data[,site.index],3]-aa.properties[optimal.aa[site.index],3])^2)^(1/2))
      aa.distances <- c(aa.distances, (1+gp[site.index]*distance))
    }else{
      aa.distances <- c(aa.distances, 0)
      gene.length <- gene.length - 1
    }
  }
  functionality = 1/((1/gene.length) * sum(aa.distances))
  return(functionality)
}


#' @title Calculate site likelihoods under SelAC
#'
#' @description
#' Calculates the likelihoods across sites and across genes under SELAC
#'
#' @param selac.obj An object of class SELAC.
#' @param codon.data.path Provides the path to the directory containing the gene specific fasta files of coding data.
#' @param aa.optim.input A list of optimal amino acids with each list element designating a character vector for each gene. The optimal amino acids be the MLE from a selac run (default) or a list of user defined optimal A.A.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.

#'
#' @details
#' The purpose of this function is to provide the site likelihoods across genes. It is also flexible in that it allows different hypotheses about optimal acids across genes and/or site. The output is a list object, with each list entry designating 1) the tot.likelihood for that gene, and 2) the site likelihoods for that gene.
GetSelacSiteLikelihoods <- function(selac.obj, codon.data.path, aa.optim.input=NULL, fasta.rows.to.keep=NULL) {

  codon.index.matrix = CreateCodonMutationMatrixIndex()
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
  n.cores <- NULL
  obj.final <- as.list(1:length(partitions))
  n.cores.by.gene.by.site <- selac.obj$n.cores.by.gene.by.site

  for(partition.index in 1:length(partitions)){
    x <- c(selac.obj$mle.pars[partition.index,])
    gene.tmp <- read.dna(partitions[partition.index], format='fasta')
    if(!is.null(fasta.rows.to.keep)){
      gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
    }else{
      gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
    }
    codon.data.tmp <- DNAbinToCodonNumeric(gene.tmp)
    codon.data.tmp <- codon.data.tmp[phy$tip.label,]
    codon.data = NULL
    codon.data$unique.site.patterns = codon.data.tmp
    nsites <- dim(codon.data$unique.site.patterns)[2]-1
    codon.data$site.pattern.counts = rep(1, nsites)

    codon.freq.by.aa=selac.obj$codon.freq.by.aa[[partition.index]]
    codon.freq.by.gene=selac.obj$codon.freq.by.gene[[partition.index]]

    if(is.null(aa.optim.input)){
      aa.optim_array <- selac.obj$aa.optim[[partition.index]]
    }

    if(include.gamma == TRUE){
      shape = x[length(x)]
      x = x[-length(x)]
    }

    C.Phi.q.Ne <- x[1]
    C <- 4
    q <- 4e-7
    Ne <- 5e6
    Phi.q.Ne <- C.Phi.q.Ne / C
    Phi.Ne <- Phi.q.Ne / q
    Phi <- Phi.Ne / Ne
    alpha <- x[2]
    beta <- x[3]

    if(k.levels > 0){
      if(nuc.model == "JC") {
        base.freqs=c(x[4:6], 1-sum(x[4:6]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
      }
      if(nuc.model == "GTR") {
        base.freqs=c(x[4:6], 1-sum(x[4:6]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[9:length(x)], model=nuc.model, base.freqs=base.freqs)
      }
      if(nuc.model == "UNREST") {
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[6:length(x)], model=nuc.model)
      }
    }else{
      if(nuc.model == "JC") {
        base.freqs=c(x[4:6], 1-sum(x[4:6]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
      }
      if(nuc.model == "GTR") {
        base.freqs=c(x[4:6], 1-sum(x[4:6]))
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[7:length(x)], model=nuc.model, base.freqs=base.freqs)
      }
      if(nuc.model == "UNREST") {
        nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[4:length(x)], model=nuc.model)
      }
    }

    codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
    codon_mutation_matrix[is.na(codon_mutation_matrix)]=0

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
      final.likelihood.mat = matrix(0, nrow=ncats, ncol=nsites)
      for(k in sequence(ncats)){
        if(k.levels > 0){
          aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=x[7:8], k=k.levels)
        }else{
          aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
        }
        Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[k], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
        final.likelihood.mat[k,] = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
      }
      likelihood <- sum(log(colSums(exp(final.likelihood.mat)*weights.k)) * codon.data$site.pattern.counts)
    }else{
      if(k.levels > 0){
        aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=x[7:8], k=k.levels)
      }else{
        aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
      }
      Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
      final.likelihood = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
      likelihood <- sum(final.likelihood * codon.data$site.pattern.counts)
    }
    obj.gene <- NULL
    obj.gene$tot.likelihood <- likelihood
    obj.gene$partial.likelihoods <- final.likelihood
    obj.final[[partition.index]] <- obj.gene
  }
  return(obj.final)
}


#' @title Phi rate category information under SELAC+gamma
#'
#' @description
#' Provides likelihood information and best rates across sites and across genes under SELAC+gamma
#'
#' @param selac.obj An object of class SELAC.
#' @param codon.data.path Provides the path to the directory containing the gene specific fasta files of coding data.
#' @param aa.optim.input A list of optimal amino acids with each list element designating a character vector for each gene. The optimal amino acids be the MLE from a selac run (default) or a list of user defined optimal A.A.
#' @param fasta.rows.to.keep Indicates which rows to remove in the input fasta files.
#' @param n.cores.by.gene.by.site The number of cores to decidate to parallelize analyses by site WITHIN a gene. Note n.cores.by.gene*n.cores.by.gene.by.site is the total number of cores dedicated to the analysis.
#'
#' @details
#' The purpose of this function is to determine which rate category best fits each site across genes. The output is a list object, with each list entry designating the optimal rate category across sites for that gene.
GetSelacPhiCat <- function(selac.obj, codon.data.path, aa.optim.input=NULL, fasta.rows.to.keep=NULL, n.cores.by.gene.by.site=1) {

  codon.index.matrix = CreateCodonMutationMatrixIndex()
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
  n.cores.by.gene.by.site <- selac.obj$n.cores.by.gene.by.site
  volume.fixed.value <- selac.obj$volume.fixed.value
  n.cores <- NULL
  obj.final <- as.list(1:length(partitions))

  for(partition.index in 1:length(partitions)){
    x <- c(selac.obj$mle.pars[partition.index,])
    if(include.gamma == TRUE){
      shape = x[length(x)]
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
    }
    gene.tmp <- read.dna(partitions[partition.index], format='fasta')
    if(!is.null(fasta.rows.to.keep)){
      gene.tmp <- as.list(as.matrix(cbind(gene.tmp))[fasta.rows.to.keep,])
    }else{
      gene.tmp <- as.list(as.matrix(cbind(gene.tmp)))
    }
    codon.data.tmp <- DNAbinToCodonNumeric(gene.tmp)
    codon.data.tmp <- codon.data.tmp[phy$tip.label,]
    codon.data = NULL
    codon.data$unique.site.patterns = codon.data.tmp
    nsites <- dim(codon.data$unique.site.patterns)[2]-1
    codon.data$site.pattern.counts = rep(1, nsites)

    codon.freq.by.aa=selac.obj$codon.freq.by.aa[[partition.index]]
    codon.freq.by.gene=selac.obj$codon.freq.by.gene[[partition.index]]

    if(is.null(aa.optim.input)){
      aa.optim_array <- selac.obj$aa.optim[[partition.index]]
    }
    phi.likelihoods.per.site <- GetPhiLikelihoodPerSite(x, codon.data=codon.data, phy=phy, aa.optim_array=aa.optim_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=numcode, diploid=diploid, aa.properties=aa.properties, volume.fixed.value=volume.fixed.value, nuc.model=nuc.model, codon.index.matrix=codon.index.matrix, include.gamma=include.gamma, gamma.type=gamma.type, ncats=ncats, k.levels=k.levels, logspace=FALSE, verbose=FALSE, neglnl=FALSE, n.cores.by.gene.by.site=1)
    C.Phi.q.Ne <- x[1]
    C <- 4
    q <- 4e-7
    Ne <- 5e6
    Phi.q.Ne <- C.Phi.q.Ne / C
    Phi.Ne <- Phi.q.Ne / q
    Phi <- Phi.Ne / Ne
    obj <- NULL
    obj$likelihood <- phi.likelihoods.per.site
    rate.vector <- c()
    indicator.raw.vector <- c()
    indicator.weighted.vector <- c()
    model.ave.phi.vector <- c()
    model.weight.ave.phi.vector <- c()
    for(i in 1:dim(phi.likelihoods.per.site)[2]){
      #model-average based solely on the raw likelihoods:
      tmp.rate.class <- which.max(phi.likelihoods.per.site[,i])
      rate.vector <- c(rate.vector, rates.k[tmp.rate.class])
      indicator.raw.vector <- c(indicator.raw.vector, tmp.rate.class)
      aic <- -2 * phi.likelihoods.per.site[,i]
      daic <- aic - aic[tmp.rate.class]
      w.aic <- exp(-.5 * daic) / sum(exp(-.5 * daic))
      model.ave.phi.vector <- c(model.ave.phi.vector, sum(Phi * rates.k * w.aic))
      #model-average based on rescaling of the likelihoods according to laguerre weights:
      tmp.rate.class <- which.max(log(exp(phi.likelihoods.per.site[,i]) * weights.k))
      indicator.weighted.vector <- c(indicator.weighted.vector, tmp.rate.class)
      aic <- -2 * log(exp(phi.likelihoods.per.site[,i]) * weights.k)
      daic <- aic - aic[tmp.rate.class]
      w.aic <- exp(-.5 * daic) / sum(exp(-.5 * daic))
      model.weight.ave.phi.vector <- c(model.weight.ave.phi.vector, sum(Phi * rates.k * w.aic))
    }
    obj$best.rate.by.site <- rate.vector
    obj$model.ave.rawLik.phi <- model.ave.phi.vector
    obj$model.ave.weightedLik.phi <- model.weight.ave.phi.vector
    obj$indicator.by.site.rawLik <- indicator.raw.vector
    obj$indicator.by.site.weightedLik <- indicator.weighted.vector
    obj.final[[partition.index]] <- obj
  }
  return(obj.final)
}



GetAALikelihoodPerSite <- function(x, codon.data, phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix, include.gamma=FALSE, gamma.type="quadrature", ncats=4, k.levels=0, logspace=FALSE, verbose=TRUE, neglnl=FALSE, n.cores.by.gene.by.site=n.cores.by.gene.by.site) {
  if(logspace) {
    x = exp(x)
  }
  if(include.gamma == TRUE){
    shape = x[length(x)]
    x = x[-length(x)]
  }

  C.Phi.q.Ne <- x[1]
  C <- 4
  q <- 4e-7
  Ne <- 5e6
  Phi.q.Ne <- C.Phi.q.Ne / C
  Phi.Ne <- Phi.q.Ne / q
  Phi <- Phi.Ne / Ne
  alpha <- x[2]
  beta <- x[3]
  gamma <- volume.fixed.value

  if(k.levels > 0){
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[9:length(x)], model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[6:length(x)], model=nuc.model, base.freqs=NULL)
      poly.params <- x[4:5]
    }
  }else{
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[7:length(x)], model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[4:length(x)], model=nuc.model, base.freqs=NULL)
    }
  }

  if(!is.null(codon.data$unique.site.patterns)){
    codon.data.list <- codon.data
    nsites.unique <- dim(codon.data$unique.site.patterns)[2]-1
  }else{
    nsites.unique <- dim(codon.data)[2]-1
    codon.data.list <- NULL
    codon.data.list$unique.site.patterns <- codon.data
    codon.data.list$site.pattern.counts <- rep(1, nsites.unique)
  }
  nsites <- sum(codon.data$site.pattern.counts)

  #codon_mutation_matrix = c(as.vector(nuc.mutation.rates), 0)[codon.index.matrix]
  codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
  codon_mutation_matrix[is.na(codon_mutation_matrix)]=0

  optimal.vector.by.site <- rep(NA, nsites.unique)
  #unique.aa <- GetMatrixAANames(numcode)
  optimal.aa.likelihood.mat <- matrix(0, nrow=length(.unique.aa), ncol=nsites.unique)

  for(i in 1:length(.unique.aa)){
    if(.unique.aa[i]=="*"){
      optimal.aa.likelihood.mat[i,] <- rep(-1000000, nsites.unique)
    }else{
      aa.optim_array = rep(.unique.aa[i], nsites.unique)
      if(include.gamma==TRUE){
        gene_site_array <- array(1, dim=c(21, nsites.unique, ncats=ncats))
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
        final.likelihood.mat = matrix(0, nrow=ncats, ncol=nsites.unique)
        for(k.cat in sequence(ncats)){
          if(k.levels > 0){
            aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
          }else{
            aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
          }
          Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[k.cat], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999)
          tmp = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data.list, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
          tmp[is.na(tmp)] = -1000000
          gene_site_array[,,k.cat] = tmp
        }
      }else{
        gene_site_array <- c()
        if(k.levels > 0){
          aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
        }else{
          aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
        }
        Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999)
        tmp = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data.list, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
        tmp[is.na(tmp)] = -1000000
        gene_site_array = tmp
      }
    }
  }
  return(gene_site_array)
}


GetPhiLikelihoodPerSite <- function(x, codon.data, phy, aa.optim_array=NULL, codon.freq.by.aa=NULL, codon.freq.by.gene=NULL, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model, codon.index.matrix, include.gamma, gamma.type, ncats, k.levels=0, logspace=FALSE, verbose=TRUE, neglnl=FALSE, n.cores.by.gene.by.site=n.cores.by.gene.by.site) {
  if(logspace) {
    x = exp(x)
  }
  if(include.gamma == TRUE){
    shape = x[length(x)]
    x = x[-length(x)]
  }

  C.Phi.q.Ne <- x[1]
  C <- 4
  q <- 4e-7
  Ne <- 5e6
  Phi.q.Ne <- C.Phi.q.Ne / C
  Phi.Ne <- Phi.q.Ne / q
  Phi <- Phi.Ne / Ne
  alpha <- x[2]
  beta <- x[3]
  gamma <- volume.fixed.value

  if(k.levels > 0){
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[9:length(x)], model=nuc.model, base.freqs=base.freqs)
      poly.params <- x[7:8]
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[6:length(x)], model=nuc.model, base.freqs=NULL)
      poly.params <- x[4:5]
    }
  }else{
    if(nuc.model == "JC") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "GTR") {
      base.freqs=c(x[4:6], 1-sum(x[4:6]))
      #During the early stages of the optimization process it will try weird values for the base frequencies.
      if(any(base.freqs < 0)){
        return(1000000)
      }
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[7:length(x)], model=nuc.model, base.freqs=base.freqs)
    }
    if(nuc.model == "UNREST") {
      nuc.mutation.rates <- CreateNucleotideMutationMatrix(x[4:length(x)], model=nuc.model, base.freqs=NULL)
    }
  }
  #codon_mutation_matrix = CreateCodonMutationMatrix(nuc.mutation.rates) #We now make an index matrix first then just place the nucleotide rates into it:
  #codon_mutation_matrix = c(as.vector(nuc.mutation.rates), 0)[codon.index.matrix]
  codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
  codon_mutation_matrix[is.na(codon_mutation_matrix)]=0
  nsites.unique <- dim(codon.data$unique.site.patterns)[2]-1
  nsites <- sum(codon.data$site.pattern.counts)

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
    final.likelihood = matrix(0, nrow=ncats, ncol=nsites.unique)
    for(k.cat in sequence(ncats)){
      if(k.levels > 0){
        aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
      }else{
        aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
      }
      Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi*rates.k[k.cat], q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
      final.likelihood[k.cat,] = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
    }
    #likelihood <- sum(log(colSums(exp(final.likelihood.mat)*weights.k)) * codon.data$site.pattern.counts)
  }else{
    if(k.levels > 0){
      aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=poly.params, k=k.levels)
    }else{
      aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=aa.properties, normalize=FALSE, poly.params=NULL, k=k.levels)
    }
    Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
    final.likelihood = GetLikelihoodSAC_CodonForManyCharVaryingBySite(codon.data, phy, Q_codon_array, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, aa.optim_array=aa.optim_array, codon_mutation_matrix=codon_mutation_matrix, Ne=Ne, rates=NULL, numcode=numcode, diploid=diploid, n.cores.by.gene.by.site=n.cores.by.gene.by.site)
    #likelihood <- sum(final.likelihood * codon.data$site.pattern.counts)
  }

  return(final.likelihood)
}


GetGeneSiteInfo <- function(selac.obj, partition.number=1){
  phy <- selac.obj$phy
  yeast.gene <- read.dna(selac.obj$partitions[partition.number], format="fasta")
  yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
  chars <- DNAbinToCodonNumeric(yeast.gene)
  codon.data <- chars[phy$tip.label,]
  aa.data <- ConvertCodonNumericDataToAAData(codon.data, numcode=1)
  aa.optim <- selac.obj$aa.optim[[1]] #starting values for all, final values for majrule
  aa.optim.full.list <- aa.optim
  codon.freq.by.aa <- GetCodonFreqsByAA(codon.data[,-1], aa.optim, numcode=1)
  codon.freq.by.gene <- GetCodonFreqsByGene(codon.data[,-1])
  aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
  colnames(aa.optim.frame.to.add) <- colnames(codon.data)
  codon.data <- rbind(codon.data, aa.optim.frame.to.add)
  codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
  aa.optim <- codon.data$optimal.aa
  codon.index.matrix <- CreateCodonMutationMatrixIndex()
  volume.fixed.value <- 0.0003990333
  n.cores.by.gene.by.site <- selac.obj$n.cores.by.gene.by.site
  pars.to.do <- selac.obj$mle.pars[partition.number,]
  selac.all.sites.amino.acid <- GetAALikelihoodPerSite(x=log(pars.to.do), codon.data, phy, aa.optim_array=aa.optim, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="UNREST", codon.index.matrix, include.gamma=TRUE, gamma.type="quadrature", ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, n.cores.by.gene.by.site=n.cores.by.gene.by.site)

  selac.all.sites.phi <- GetPhiLikelihoodPerSite(x=log(pars.to.do), codon.data, phy, aa.optim_array=aa.optim, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="UNREST", codon.index.matrix, include.gamma=TRUE, gamma.type="quadrature", ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, n.cores.by.gene.by.site=n.cores.by.gene.by.site)

  obj <- NULL
  quadl <- LaguerreQuad(selac.obj$mle.pars[1,length(selac.obj$mle.pars[1,])],4)
  obj$phi <- selac.obj$mle.pars[partition.number,1] * quadl[1:4]
  obj$phi.weights <- quadl[5:8]
  obj$site.aa.information <- selac.all.sites.amino.acid
  obj$site.phi.information <- selac.all.sites.phi
  return(obj)
}

#library(selac)
#load("yeastRokasSelacUNRESTgamma.Rdata")
#pp <- GetGeneSiteInfo(result, partition.number=1)


######################################################################################################################################
######################################################################################################################################
### Print function for the selac class:
######################################################################################################################################
######################################################################################################################################

print.selac <- function(x,...){
  ntips=Ntip(x$phy)
  output<-data.frame(x$loglik,x$AIC,ntips,sum(x$nsites), x$k.levels, row.names="")
  names(output)<-c("-lnL","AIC", "ntax", "nsites", "k.levels")
  cat("\nFit\n")
  print(output)
  cat("\n")
  cat("\nModel options\n")
  output.part.deux <- data.frame(x$nuc.model, x$data.type, x$aa.optim.type, x$include.gamma, x$ncats, row.names="")
  names(output.part.deux) <- c("model","data", "opt.aa?", "disc.gamma", "n.cats")
  print(output.part.deux)
  cat("\n")
  cpv.starting.parameters <- GetAADistanceStartingParameters(aa.properties=x$aa.properties)
  if(x$aa.optim.type=="majrule" | x$aa.optim.type=="optimize"){
    if(x$nuc.model == "JC"){
      cat("\nSELAC Parameters\n")
      if(x$include.gamma==TRUE){
        if(x$k.levels > 1){
          output<-data.frame(x$mle.pars[1,2], x$mle.pars[1,3], x$volume.fixed.value, row.names="")
        }else{
          output<-data.frame(x$mle.pars[1,2], x$mle.pars[1,3], x$volume.fixed.value, row.names="")
        }
        names(output)<-c("c","p","v")
      }else{
        output<-data.frame(x$mle.pars[1,2], x$mle.pars[1,3], x$volume.fixed.value, row.names="")
        names(output)<-c("c","p","v")
      }
      print(output)
      cat("\n")
    }
    if(x$nuc.model == "GTR"){
      cat("\nSELAC parameters\n")
      #if(x$include.gamma==TRUE){
      #    if(x$k.levels > 1){
      #        output <- data.frame(x$mle.pars[1,2], x$mle.pars[1,3], x$volume.fixed.value, x$mle.pars[1,14], row.names="")
      #    }else{
      #        output <- data.frame(x$mle.pars[1,2], x$mle.pars[1,3], x$volume.fixed.value, x$mle.pars[1,12], row.names="")
      #    }
      #    names(output)<-c("c","p","v","disc.gamma")
      #}else{
      output<-data.frame(x$mle.pars[1,2], x$mle.pars[1,3], x$volume.fixed.value, row.names="")
      names(output)<-c("c","p","v")
      #}
      print(output)
      cat("\n")

    }
  }
}
