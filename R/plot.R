#' Function to plot frequency of distribution of different Wi given selac parameters

ComputeEquilibriumCodonFrequencies <- function(nuc.model="JC", base.freqs=rep(0.25, 4), nsites=1, C=4, Phi=0.5, q=4e-7, Ne=5e6, alpha=1.83, beta=0.10, gamma=0.0003990333, include.stop.codon=TRUE, numcode=1, diploid=TRUE, flee.stop.codon.rate=0.9999999) {
#To test: nuc.model="JC"; base.freqs=rep(0.25, 4); nsites=1; C=4; Phi=0.5; q=4e-7; Ne=5e6; alpha=1.83; beta=0.10; gamma=0.0003990333; include.stop.codon=TRUE; numcode=1; diploid=TRUE; flee.stop.codon.rate=0.9999999
  nuc.mutation.rates <- selac:::CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
  codon.index.matrix = selac:::CreateCodonMutationMatrixIndex()
  codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
  codon_mutation_matrix[is.na(codon_mutation_matrix)]=0
  aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=NULL, normalize=FALSE, poly.params=NULL, k=0)
  Q_codon_array <- selac:::FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=include.stop.codon, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
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
#' @examples
#'
#' phi.vector <- c(0.01, .1, 0.5, 2)
#' eq.freq.matrices <- array(dim=c(64, 20, length(phi.vector)))
#' for (i in sequence(length(phi.vector))) {
#'   eq.freq.matrices[,,i] <- ComputeEquilibriumCodonFrequencies(Phi=phi.vector[i])
#' }
#' values = paste("Phi = ", phi.vector, sep="")
#' PlotEquilbriumCodonDistribution(eq.freq.matrices, values)
PlotEquilbriumCodonDistribution <- function(eq.freq.matrices, values, palette="Set1", lwd=2, ...) {
  colors <- RColorBrewer::brewer.pal(dim(eq.freq.matrices)[3],palette)
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

ComputeEquilibriumAAFitness <- function(nuc.model="JC", base.freqs=rep(0.25, 4), nsites=1, C=4, Phi=0.5, q=4e-7, Ne=5e6, alpha=1.83, beta=0.10, gamma=0.0003990333, include.stop.codon=TRUE, numcode=1, diploid=TRUE, flee.stop.codon.rate=0.9999999) {
  eq.freq.matrix <- ComputeEquilibriumCodonFrequencies(nuc.model, base.freqs, nsites, C, Phi, q, Ne, alpha, beta, gamma, include.stop.codon, numcode, diploid, flee.stop.codon.rate)
  codon.fitness.matrix <- eq.freq.matrix*NA
  aa.distances <- selac:::CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=NULL, normalize=FALSE, poly.params=NULL, k=0)
  aa.names <- unique(sapply(rownames(eq.freq.matrix), selac:::TranslateCodon, numcode=numcode))
  aa.fitnesses <- matrix(nrow=length(aa.names), ncol=dim(eq.freq.matrix)[2])
  rownames(aa.fitnesses) <- aa.names
  colnames(aa.fitnesses) <- colnames(eq.freq.matrix)
  for (col.index in sequence(dim(codon.fitness.matrix)[2])) {
    for (row.index in sequence(dim(codon.fitness.matrix)[1])) {
      codon.fitness.matrix[row.index, col.index] <- selac:::GetFitness(focal.protein=selac:::TranslateCodon(rownames(codon.fitness.matrix)[row.index], numcode), optimal.protein=colnames(codon.fitness.matrix)[col.index], aa.distances, nsites=nsites, C=1, Phi=Phi, q=q)
      aa.fitnesses[selac:::TranslateCodon(rownames(codon.fitness.matrix)[row.index], numcode), col.index] <- codon.fitness.matrix[row.index, col.index]
    }
  }
  return(list(aa.fitness.matrix=aa.fitnesses, codon.fitnesses=codon.fitness.matrix, equilibrium.codon.frequency = eq.freq.matrix))
}


#' Function to plot a distribution of fitnesses W or selection coefficients S for a given optimal aa and other terms.
#'
#' @param aa.fitness.matrices, A 3d array of aa.fitness.matrix returned from ComputeEquilibriumAAFitness (first element in return)
#' @param values The vector of labels for each matrix (i.e., different Phi values)
#' @param optimal.aa Single letter code for the optimal aa. If NULL, integrates across aa.
#' @param palette Color palette to use from RColorBrewer
#' @param lwd Line width
#' @param include.stop.codon Include stop codons
#' @param type If "histogram", do a histogram plot; if "density", do a density plot
#' @param fitness If TRUE, plot fitness W; if FALSE, plot selection coefficient S (= W- 1)
#' @param scale.x.axis.by.Ne if TRUE, x axis is transformed from S to S*Ne; if FALSE no scaling is done
#' @param Ne used to scale x axis when scale.x.axis.by.Ne is TRUE
#' @param ... Other paramters to pass to plot()
#' @examples
#' phi.vector <- c(0.0000000001, 0.01, .1, 0.5, 2)
#' aa.fitness.matrices <- array(dim=c(21, 20, length(phi.vector)))
#' for (i in sequence(length(phi.vector))) {
#'  local.matrix <- ComputeEquilibriumAAFitness(Phi=phi.vector[i])$aa.fitness.matrix
#'  aa.fitness.matrices[,,i] <- local.matrix
#'  dimnames(aa.fitness.matrices) <- list(rownames(local.matrix), colnames(local.matrix), NULL)
#' }
#' values = paste("Phi = ", phi.vector, sep="")
#' PlotPerAAFitness(aa.fitness.matrices, values, optimal.aa="L")
PlotPerAAFitness <- function(aa.fitness.matrices, values, optimal.aa=NULL, palette="Set1", lwd=2, include.stop.codon=FALSE, type="histogram", fitness=TRUE, scale.x.axis.by.Ne=FALSE, legend.title=NULL,Ne=10^6, ...) {
  colors <- RColorBrewer::brewer.pal(dim(aa.fitness.matrices)[3],palette)
  distributions <- list()
  y.range <- c()

  ##scale values based on 'fitness' and set legend label if necessary
  if(fitness) {
      if(is.null(legend.title)) {
          legend.title="W"
      }
  }else{
      aa.fitness.matrices <- aa.fitness.matrices -1 #S_i = W_i - W_*
      if(is.null(legend.title)) {
          legend.title="s"
      }
      if(scale.x.axis.by.Ne){
          aa.fitness.matrices <- aa.fitness.matrices * Ne
          if(is.null(legend.title)) {
              legend.title="s Ne"
          }
      }
  }
  x.range <- c(NA)
  for (i in sequence(dim(aa.fitness.matrices)[3])) {
    distribution <- NA
    local.matrix <- aa.fitness.matrices[,,i]
    if(!include.stop.codon) {
      local.matrix <- local.matrix[which(rownames(local.matrix) != "*"),]
    }
    input.values <- NA
    if (is.null(optimal.aa)) {
      input.values <- c(local.matrix)
    } else {
      input.values <- c(local.matrix[,optimal.aa])
    }
    distribution <- list()
    if(type=="density") {
      distribution <- stats::density(input.values, to=1)
      distribution$y <- distribution$y / sum(distribution$y)
      distributions[[i]] <- distribution
    }
    if(type=="histogram") {
    #  table.of.dist <- table(input.values)

     distribution <- graphics::hist(input.values, plot=FALSE)
      distribution$x <- distribution$mids
      distribution$y <- distribution$counts / length(input.values)
      if(length(distribution$x)==1) {
        distribution$x <- median(input.values)
        distribution$y <- 1
      }
    #  distribution$x <- as.numeric(names(table.of.dist))
    #  distribution$y <- unname(table.of.dist)/length(input.values)
      distributions[[i]] <- distribution
    }
    y.range <- range(c(y.range, distribution$y), na.rm=TRUE)
    x.range <- range(c(x.range, distribution$x), na.rm=TRUE)
  }
  plot(x=x.range, y=y.range, type="n", bty="n", xlab=ifelse(fitness, "W", ifelse(scale.x.axis.by.Ne, "s Ne", "s")), ylab="Frequency", ...)
  if(type="histogram") {
    for (i in sequence(length(distributions))) {
      points(distributions[[i]]$x, distributions[[i]]$y, col=colors[i], pch=20)
      for (j in sequence(length(distributions[[i]]$x))) {
        lines(rep(distributions[[i]]$x[j],2), c(0, distributions[[i]]$y[j]), col=add.alpha(colors[i],0.2), lwd=lwd)
      }
    }
  } else {
    lines(distributions[[i]]$x, distributions[[i]]$y, col=colors[i], lwd=lwd)
  }
  legend(x="topleft", legend=values, fill=colors)
}

# from http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
add.alpha <- function(col, alpha=1){
if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
function(x)
     rgb(x[1], x[2], x[3], alpha=alpha))
}

#' Function to plot a distribution of fitnesses based on codon equilibrium freqs
#'
#' @param aa.fitness.matrices, A 3d array of aa.fitness.matrix returned from ComputeEquilibriumAAFitness (first element in return)
#' @param values The vector of labels for each matrix (i.e., different Phi values)
#' @param optimal.aa Single letter code for the optimal aa. If NULL, integrates across aa.
#' @param palette Color palette to use from RColorBrewer
#' @param lwd Line width
#' @param include.stop.codon Include stop codons
#' @param type If "histogram", do a histogram plot; if "density", do a density plot
#' @param fitness If TRUE, plot W; if FALSE, plot S (= 1 - W)
#' @param numcode The genetic code
#' @param ... Other paramters to pass to plot()
#' @examples
#' phi.vector <- c(.1, 0.5, 2)
#' codon.fitnesses.matrices <- array(dim=c(64, 20, length(phi.vector)))
#' codon.eq.matrices <- array(dim=c(64, 20, length(phi.vector)))
#' for (i in sequence(length(phi.vector))) {
#'    local.matrix <- ComputeEquilibriumAAFitness(Phi=phi.vector[i])
#'    codon.fitnesses.matrices[,,i] <-  local.matrix$codon.fitnesses
#'    codon.eq.matrices[,,i] <-  local.matrix$equilibrium.codon.frequency
#'    dimnames(codon.fitnesses.matrices) <- list(rownames(local.matrix$codon.fitnesses), colnames(local.matrix$codon.fitnesses), NULL)
#'    dimnames(codon.eq.matrices) <- list(rownames(local.matrix$equilibrium.codon.frequency), colnames(local.matrix$equilibrium.codon.frequency), NULL)
#' }
#' values = paste("Phi = ", phi.vector, sep="")
#' PlotExpectedFitness(codon.fitnesses.matrices, codon.eq.matrices, values)
PlotExpectedFitness <- function(codon.fitnesses.matrices, codon.eq.matrices, values, optimal.aa=NULL, palette="Set1", lwd=2, include.stop.codon=FALSE, type="histogram", fitness=TRUE, numcode=1, ...) {
  colors <- RColorBrewer::brewer.pal(dim(codon.fitnesses.matrices)[3],palette)
  distributions <- list()
  y.range <- c()
  if(!fitness) {
    codon.fitnesses.matrices <- 1 - codon.fitnesses.matrices
  }
  x.range <- c(NA)
  for (i in sequence(dim(codon.fitnesses.matrices)[3])) {
    distribution <- NA

    local.codon.fitness.matrix <- codon.fitnesses.matrices[,,i]
    local.codon.equilibrium.matrix <- codon.eq.matrices[,,i]
    if(!include.stop.codon) {
      aa.vector <- unname(sapply(rownames(local.codon.fitness.matrix), TranslateCodon, numcode=numcode))
      local.codon.fitness.matrix <- local.codon.fitness.matrix[which(aa.vector != "*"),]
      local.codon.equilibrium.matrix <- local.codon.equilibrium.matrix[which(aa.vector != "*"),]
    }

    input.values <- NA
    if (is.null(optimal.aa)) {
      input.values <- rep(c(local.codon.fitness.matrix), round(10000*c(local.codon.equilibrium.matrix)))
    } else {
      input.values <- rep(c(local.codon.fitness.matrix[,optimal.aa]), round(10000*c(local.codon.equilibrium.matrix[,optimal.aa])))
    }
    distribution <- list()
    if(type=="density") {
      distribution <- stats::density(input.values, to=1)
      distribution$y <- distribution$y / sum(distribution$y)
      distributions[[i]] <- distribution
    }
    if(type=="histogram") {
    #  table.of.dist <- table(input.values)

     distribution <- graphics::hist(input.values, plot=FALSE)
      distribution$x <- distribution$mids
      distribution$y <- distribution$counts / length(input.values)
      if(length(distribution$x)==1) {
        distribution$x <- median(input.values)
        distribution$y <- 1
      }
    #  distribution$x <- as.numeric(names(table.of.dist))
    #  distribution$y <- unname(table.of.dist)/length(input.values)
      distributions[[i]] <- distribution
    }
    y.range <- range(c(y.range, distribution$y), na.rm=TRUE)
    x.range <- range(c(x.range, distribution$x), na.rm=TRUE)
  }
  plot(x=x.range, y=y.range, type="n", bty="n", xlab=ifelse(fitness, "W", "S"), ylab="Frequency", ...)
  if(type=="histogram") {
    for (i in sequence(length(distributions))) {
      points(distributions[[i]]$x, distributions[[i]]$y, col=colors[i], pch=20)
      for (j in sequence(length(distributions[[i]]$x))) {
        lines(rep(distributions[[i]]$x[j],2), c(0, distributions[[i]]$y[j]), col=add.alpha(colors[i],0.2), lwd=lwd)
      }
    }
  } else {
    lines(distributions[[i]]$x, distributions[[i]]$y, col=colors[i], lwd=lwd)
  }
  legend(x="topleft", legend=values, fill=colors)
}

#' Function to plot info by site in a gene
#'
#' @param info.by.site The output of GetGeneSiteInfo
PlotGeneSiteInfo <- function(info.by.site, aa.properties=NULL) {
  if(is.null(aa.properties)) {
    aa.properties <- structure(c(0, 2.75, 1.38, 0.92, 0, 0.74, 0.58, 0, 0.33, 0, 0,
    1.33, 0.39, 0.89, 0.65, 1.42, 0.71, 0, 0.13, 0.2, 8.1, 5.5, 13,
    12.3, 5.2, 9, 10.4, 5.2, 11.3, 4.9, 5.7, 11.6, 8, 10.5, 10.5,
    9.2, 8.6, 5.9, 5.4, 6.2, 31, 55, 54, 83, 132, 3, 96, 111, 119,
    111, 105, 56, 32.5, 85, 124, 32, 61, 84, 170, 136), .Dim = c(20L,
    3L), .Dimnames = list(c("A", "C", "D", "E", "F", "G",
    "H", "I", "K", "L", "M", "N", "P", "Q", "R",
    "S", "T", "V", "W", "Y"), c("c", "p", "v"))) #properties from Grantham paper
  }
  #TODO: filter out stop aa weights
  stop("you have to finish writing code to pull out stop aa")
  #END TODO
  AIC.site.information <- -2*info.by.site$site.information
  get.delta <- function(x) {
    return(x-min(x))
  }
  normalize.rel <- function(x) {
    return(x/sum(x))
  }
  delta.AIC.site.information <- apply(AIC.site.information, c(2,3), get.delta)
  rel.likelihood.site.information <- exp(-0.5* delta.AIC.site.information)
  weight.site.information <- apply(rel.likelihood.site.information, c(2,3), normalize.rel)
  dimnames(weight.site.information)[1] <- .unique.aa
  #TODO: get weighted estimate of aa properties for given phi
  #TODO: Get weighted esetimate (weighted by phi weight) across phi
  #TODO: Plot these three properties
  #TODO: Plot sliding window of them
  #TODO: Plot average weight of Phi
}
