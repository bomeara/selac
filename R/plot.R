#' Function to plot frequency of distribution of different Wi given selac parameters

ComputeEquilibriumFrequencies <- function(nuc.model="JC", base.freqs=rep(0.25, 4), nsites=1, C=4, Phi=0.5, q=4e-7, Ne=5e6, alpha=1.83, beta=0.10, gamma=0.0003990333, volume.fixed.value=0003990333, include.stop.codon=TRUE, numcode=1, diploid=TRUE, flee.stop.codon.rate=0.9999999) {
#To test: nuc.model="JC"; base.freqs=rep(0.25, 4); nsites=1; C=4; Phi=0.5; q=4e-7; Ne=5e6; alpha=1.83; beta=0.10; gamma=0.0003990333; volume.fixed.value=0003990333; include.stop.codon=TRUE; numcode=1; diploid=TRUE; flee.stop.codon.rate=0.9999999
  nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
  codon.index.matrix = selac:::CreateCodonMutationMatrixIndex()
  codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
  codon_mutation_matrix[is.na(codon_mutation_matrix)]=0
  aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=NULL, normalize=FALSE, poly.params=NULL, k=0)
  Q_codon_array <- selac:::FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
  diag(codon_mutation_matrix) = 0
  diag(codon_mutation_matrix) = -rowSums(codon_mutation_matrix)
  .unique.aa <- c("K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E", "D", "A", "G", "V", "*", "Y", "C", "W", "F")

  #Finish the Q_array codon mutation matrix multiplication here:
  for(k in 1:21){
      if(diploid == TRUE){
          Q_codon_array[,,.unique.aa[k]] = (2 * Ne) * codon_mutation_matrix * Q_codon_array[,,.unique.aa[k]]
      }else{
          Q_codon_array[,,.unique.aa[k]] = Ne * codon_mutation_matrix * Q_codon_array[,,.unique.aa[k]]
      }
      diag(Q_codon_array[,,.unique.aa[k]]) = 0
      diag(Q_codon_array[,,.unique.aa[k]]) = -rowSums(Q_codon_array[,,.unique.aa[k]])
  }

  starting.state <- 1
  liks <- matrix(1/64, 1, 64)
  #liks[,starting.state] <- 1
  .nonstop.aa <- .unique.aa[-which(.unique.aa=="*")]
  #eq.freq <- expm(Q_codon_array[,,optimal.aa] * 1000000, method=c("Ward77")) %*% liks[1,]
  eq.freq.matrix <- matrix(nrow=64, ncol=length(.nonstop.aa))
  rownames(eq.freq.matrix) <- rownames(Q_codon_array)
  colnames(eq.freq.matrix) <- .nonstop.aa
  for (i in sequence(length(.nonstop.aa))) {
    eq.freq.matrix[,i] <- expm(Q_codon_array[,,.nonstop.aa[i]] * 1000000, method=c("Ward77")) %*% liks[1,]
  }
  return(eq.freq.matrix)
  #fitnesses <- GetFitness(focal.protein, optimal.protein, aa.distances, nsites, C, Phi, q, Ne, diploid)
}



# CODE FROM JEREMY


#
# library(selac)
# library(expm)
#
# rates <- c()
# base.freqs <- rep(.25,4)
# nuc.model <- "JC"
# nsites <- 400
# diploid <- FALSE
# numcode <- 1
#
# #Assumes JC rates:
# input.pars <- c(2, 1.83, .10)
#
# ###############
# C.Phi.q.Ne <- input.pars[1]
# C <- 4
# q <- 4e-7
# Ne <- 5e6
# Phi.q.Ne <- C.Phi.q.Ne / C
# Phi.Ne <- Phi.q.Ne / q
# Phi <- Phi.Ne / Ne
# alpha <- input.pars[2]
# beta <- input.pars[3]
# gamma <- volume.fixed.value <- 0.0003990333
#
#
# #Step 1: Get nucleotide rates:
#
# nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
#
# codon.index.matrix = selac:::CreateCodonMutationMatrixIndex()
# codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
# codon_mutation_matrix[is.na(codon_mutation_matrix)]=0
#
# aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=NULL, normalize=FALSE, poly.params=NULL, k=0)
# Q_codon_array <- selac:::FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=TRUE, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
#
# diag(codon_mutation_matrix) = 0
# diag(codon_mutation_matrix) = -rowSums(codon_mutation_matrix)
#
# ###### Not sure these lines are necessary ########
# #scale.factor <- -sum(diag(codon_mutation_matrix) * codon.freq.by.gene, na.rm=TRUE)
# #codon_mutation_matrix_scaled = codon_mutation_matrix * (1/scale.factor)
# codon_mutation_matrix_scaled <- codon_mutation_matrix
# ##################################################
#
# .unique.aa <- c("K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E", "D", "A", "G", "V", "*", "Y", "C", "W", "F")
#
# #Finish the Q_array codon mutation matrix multiplication here:
# for(k in 1:21){
#     if(diploid == TRUE){
#         Q_codon_array[,,.unique.aa[k]] = (2 * Ne) * codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[k]]
#     }else{
#         Q_codon_array[,,.unique.aa[k]] = Ne * codon_mutation_matrix_scaled * Q_codon_array[,,.unique.aa[k]]
#     }
#     diag(Q_codon_array[,,.unique.aa[k]]) = 0
#     diag(Q_codon_array[,,.unique.aa[k]]) = -rowSums(Q_codon_array[,,.unique.aa[k]])
# }
#
# starting.state <- 1
# liks <- matrix(0, 1, 64)
# liks[,starting.state] <- 1
# eq.freq <- expm(Q_codon_array[,,optimal.aa] * 1000000, method=c("Ward77")) %*% liks[1,]
#
#
#
#
#
#
#
