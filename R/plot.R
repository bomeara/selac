#' Function to plot frequency of distribution of different Wi given selac parameters

ComputeEquilibriumFrequencies <- function(nuc.model="JC", base.freqs=rep(0.25, 4), nsites=1, C=4, Phi=0.5, q=4e-7, Ne=5e6, alpha=1.83, beta=0.10, gamma=0.0003990333, include.stop.codon=TRUE, numcode=1, diploid=TRUE, flee.stop.codon.rate=0.9999999) {
#To test: nuc.model="JC"; base.freqs=rep(0.25, 4); nsites=1; C=4; Phi=0.5; q=4e-7; Ne=5e6; alpha=1.83; beta=0.10; gamma=0.0003990333; include.stop.codon=TRUE; numcode=1; diploid=TRUE; flee.stop.codon.rate=0.9999999
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
    eq.freq.matrix[,i] <- expm::expm(t(Q_codon_array[,,.nonstop.aa[i]]) * 1000000, method=c("Ward77")) %*% liks[1,]
  }
  return(eq.freq.matrix)
}

#' Function to plot a distribution of frequencies of codons given a 3d array of equilibrium frequency matrices
#'
#' @param eq.freq.matrices A 3d array of eq.freq.matrix returned from ComputeEquilibriumFrequencies
#' @param values The vector of labels for each matrix (i.e., different Phi values)
#' @param palette Color palette to use from RColorBrewer
#' @param lwd Line width
#' @param ... Other paramters to pass to plot()
#' @example
#' phi.vector <- c(0.01, .1, 0.5, 2)
#' eq.freq.matrices <- array(dim=c(64, 20, length(phi.vector)))
#' for (i in sequence(length(phi.vector))) {
#'   eq.freq.matrices[,,i] <- ComputeEquilibriumFrequencies(Phi=phi.vector[i])
#' }
#' values = paste("Phi = ", phi.vector, sep="")
#' PlotEquilbriumDistribution(eq.freq.matrices, values)
PlotEquilbriumDistribution <- function(eq.freq.matrices, values, palette="Set1", lwd=2, ...) {
  library(RColorBrewer) #MOVE TO NAMESPACE
  colors <- brewer.pal(dim(eq.freq.matrices)[3],palette)
  distributions <- list()
  total.range <- c(NA)
  for (i in sequence(dim(eq.freq.matrices)[3])) {
    sorted <- apply(eq.freq.matrices[,,i], 2, sort, decreasing=TRUE)
    distributions[[i]] <- apply(sorted, 1, mean)
    total.range <- range(c(distributions[i], total.range), na.rm=TRUE)
  }
  plot(x=c(1,64), y=total.range, type="n", bty="n", xlab="index", ylab="Average frequency", ...)
  for (i in sequence(dim(eq.freq.matrices)[3])) {
    lines(sequence(length(distributions[[i]])), distributions[[i]], lwd=lwd, col=colors[i])
  }
  legend(x="topright", legend=values, fill=colors)
}
