# Function to plot frequency of distribution of different Wi given selac parameters

ComputeEquilibriumCodonFrequencies <- function(nuc.model="JC", base.freqs=rep(0.25, 4), nsites=1, C=4, Phi=0.5, q=4e-7, Ne=5e6, alpha=1.83, beta=0.10, gamma=0.0003990333, include.stop.codon=TRUE, numcode=1, diploid=TRUE, flee.stop.codon.rate=0.9999999) {
#To test: nuc.model="JC"; base.freqs=rep(0.25, 4); nsites=1; C=4; Phi=0.5; q=4e-7; Ne=5e6; alpha=1.83; beta=0.10; gamma=0.0003990333; include.stop.codon=TRUE; numcode=1; diploid=TRUE; flee.stop.codon.rate=0.9999999
  nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
  codon.index.matrix = CreateCodonMutationMatrixIndex()
  codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
  codon_mutation_matrix[is.na(codon_mutation_matrix)]=0
  aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=NULL, normalize=FALSE, poly.params=NULL, k=0)
  Q_codon_array <- FastCreateAllCodonFixationProbabilityMatrices(aa.distances=aa.distances, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, include.stop.codon=include.stop.codon, numcode=numcode, diploid=diploid, flee.stop.codon.rate=0.9999999)
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
  aa.distances <- CreateAADistanceMatrix(alpha=alpha, beta=beta, gamma=gamma, aa.properties=NULL, normalize=FALSE, poly.params=NULL, k=0)
  aa.names <- unique(sapply(rownames(eq.freq.matrix), TranslateCodon, numcode=numcode))
  aa.fitnesses <- matrix(nrow=length(aa.names), ncol=dim(eq.freq.matrix)[2])
  rownames(aa.fitnesses) <- aa.names
  colnames(aa.fitnesses) <- colnames(eq.freq.matrix)
  for (col.index in sequence(dim(codon.fitness.matrix)[2])) {
    for (row.index in sequence(dim(codon.fitness.matrix)[1])) {
      codon.fitness.matrix[row.index, col.index] <- GetFitness(focal.protein=TranslateCodon(rownames(codon.fitness.matrix)[row.index], numcode), optimal.protein=colnames(codon.fitness.matrix)[col.index], aa.distances, nsites=nsites, C=1, Phi=Phi, q=q)
      aa.fitnesses[TranslateCodon(rownames(codon.fitness.matrix)[row.index], numcode), col.index] <- codon.fitness.matrix[row.index, col.index]
    }
  }
  return(list(aa.fitness.matrix=aa.fitnesses, codon.fitnesses=codon.fitness.matrix, equilibrium.codon.frequency = eq.freq.matrix))
}


# Computes the distribution of fitness differences between new and original mutations as well as the frequency with which those are attempted
# Returns a list; most of the return objects are the same as from ComputeEquilibriumAAFitness() but the new things are
# codon.mutation.matrix: instantaneous rate matrix for codon mutations
# codon.relative.rate.matrix: the above, but scaled so that each row sums to 1
# delta.fitness.array: 3d array: dimensions are starting codon, optimal aa, and finishing codon; entries are new codon fitness - original codon fitness
# frequency.array: 3d array: dimensions are starting codon, optimal aa, and finishing codon; entries are frequencies that that mutation is attempted given the optimal aa (column)
#
ComputeMutationFitnesses <- function(nuc.model="JC", base.freqs=rep(0.25, 4), nsites=1, C=4, Phi=0.5, q=4e-7, Ne=5e6, alpha=1.83, beta=0.10, gamma=0.0003990333, include.stop.codon=TRUE, numcode=1, diploid=TRUE, flee.stop.codon.rate=0.9999999) {
  equilibrium.values <- ComputeEquilibriumAAFitness(nuc.model=nuc.model, base.freqs=base.freqs, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, alpha=alpha, beta=beta, gamma=gamma, include.stop.codon=include.stop.codon, numcode=numcode, diploid=diploid, flee.stop.codon.rate=flee.stop.codon.rate)
  codon.fitness.matrix <- equilibrium.values$codon.fitnesses
  equilibrium.codon.frequency <- equilibrium.values$equilibrium.codon.frequency
  nuc.mutation.rates <- CreateNucleotideMutationMatrix(1, model=nuc.model, base.freqs=base.freqs)
  codon.index.matrix = CreateCodonMutationMatrixIndex()
  codon_mutation_matrix <- matrix(nuc.mutation.rates[codon.index.matrix], dim(codon.index.matrix))
  codon_mutation_matrix[is.na(codon_mutation_matrix)]=0
  diag(codon_mutation_matrix) = 0
  scaling.factors <- rowSums(codon_mutation_matrix)
  diag(codon_mutation_matrix) = -rowSums(codon_mutation_matrix)
  codon.relative.rate.matrix <- codon_mutation_matrix
  diag(codon.relative.rate.matrix) <- 0
  codon.relative.rate.matrix <- codon.relative.rate.matrix / scaling.factors #so we know the freq of each move
  delta.fitness.array <- array(dim=c(64,20,64)) #row is which codon coming from, col is which aa is optimal, depth is which codon going to
  dimnames(delta.fitness.array) <- list(rownames(equilibrium.codon.frequency), colnames(equilibrium.codon.frequency), rownames(equilibrium.codon.frequency))
  frequency.array <- delta.fitness.array
  for (optimal.aa.index in sequence(20)) {
    for (from.codon.index in sequence(64)) {
      for (to.codon.index in sequence(64)) {
        frequency.array[from.codon.index, optimal.aa.index, to.codon.index] <- codon.relative.rate.matrix[from.codon.index, to.codon.index] * equilibrium.codon.frequency[from.codon.index, optimal.aa.index]
        delta.fitness.array[from.codon.index, optimal.aa.index, to.codon.index] <- codon.fitness.matrix[to.codon.index, optimal.aa.index] - codon.fitness.matrix[from.codon.index, optimal.aa.index] #new codon - old codon
      }
    }
  }
  return(list(aa.fitness.matrix=equilibrium.values$aa.fitnesses, codon.fitnesses=equilibrium.values$codon.fitnesses, equilibrium.codon.frequency = equilibrium.values$equilibrium.codon.frequency, codon.mutation.matrix=codon_mutation_matrix, codon.relative.rate.matrix=codon.relative.rate.matrix, delta.fitness.array=delta.fitness.array, frequency.array=frequency.array))
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
#' @param legend.title Sets the title of the figure legend.
#' @param Ne used to scale x axis when scale.x.axis.by.Ne is TRUE
#' @param ... Other paramters to pass to plot()
PlotPerAAFitness <- function(aa.fitness.matrices, values, optimal.aa=NULL, palette="Set1", lwd=2, include.stop.codon=FALSE, type="histogram", fitness=TRUE, scale.x.axis.by.Ne=FALSE, legend.title=NULL, Ne=10^6, ...) {
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
#' @param codon.fitnesses.matrices A 3d array of aa.fitness.matrix returned from ComputeEquilibriumAAFitness (first element in return)
#' @param codon.eq.matrices A 3d array of codon equilibrium frequencies
#' @param values The vector of labels for each matrix (i.e., different Phi values)
#' @param optimal.aa Single letter code for the optimal aa. If NULL, integrates across aa.
#' @param palette Color palette to use from RColorBrewer
#' @param lwd Line width
#' @param include.stop.codon Include stop codons
#' @param type If "histogram", do a histogram plot; if "density", do a density plot
#' @param fitness If TRUE, plot W; if FALSE, plot S (= 1 - W)
#' @param numcode The genetic code
#' @param ... Other paramters to pass to plot()
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

# Generates a line for mutation fitness spectra
LineMutationFitnessSpectra <- function(mutation.fitness.object, optimal.aa=NULL) {
  delta.fitness.array <- mutation.fitness.object$delta.fitness.array
  frequency.array <- mutation.fitness.object$frequency.array
  if(!is.null(optimal.aa)) {
    delta.fitness.array <- delta.fitness.array[,which(dimnames(delta.fitness.array)[[2]]==optimal.aa),]
    frequency.array <- frequency.array[,which(dimnames(frequency.array)[[2]]==optimal.aa),]
  }
  result <- stats::density(x=c(delta.fitness.array), weights=c(frequency.array)/sum(frequency.array))
  return(result)
}

#' Plot fitness of mutations, weighted by frequency of those mutations
#'
#' @param mutation.fitness.object.list List that contains multiple objects from ComputeMutationFitnesses() calls
#' @param values The vector of labels for each matrix (i.e., different Phi values)
#' @param optimal.aa Single letter code for the optimal aa. If NULL, integrates across aa.
#' @param palette Color palette to use from RColorBrewer
#' @param lwd Line width
#' @param ... other arguments to pass to plot()
#'
PlotMutationFitnessSpectra <- function(mutation.fitness.object.list, values, optimal.aa=NULL, palette="Set1", lwd=2, ...) {
  colors <- add.alpha(RColorBrewer::brewer.pal(length(mutation.fitness.object.list),palette),0.5)
  results.to.plot <- lapply(mutation.fitness.object.list, LineMutationFitnessSpectra, optimal.aa=optimal.aa)
  x.range <- c(NA)
  y.range <- c(NA)
  for (i in sequence(length(results.to.plot))) {
    x.range <- range(results.to.plot[[i]]$x, na.rm=TRUE)
    y.range <- range(results.to.plot[[i]]$y, na.rm=TRUE)
  }
  plot(x=x.range, y=y.range, bty="n", xlab="W", ylab="density", type="n",...)
  for (i in sequence(length(results.to.plot))) {
    lines(results.to.plot[[i]], lwd=lwd,col=colors[i])
  }
  legend(x="topleft", legend=values, fill=colors)
}


#' Function to plot info by site in a gene
#'
#' @param all.info The output of GetGeneSiteInfo
#' @param aa.properties The aa.properties you want to use; if NULL, uses Grantham
#' @param mean.width Sliding window width
PlotGeneSiteInfo <- function(all.info, aa.properties=NULL, mean.width=10) {
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

  aa.properties.reordered <- aa.properties[order(match(rownames(aa.properties), .unique.aa)),]

  info.by.site <- all.info$site.aa.information
  dimnames(info.by.site)[[1]] <- .unique.aa
  info.by.site <- info.by.site[which(.unique.aa!="*"),,]
  AIC.site.information <- -2*info.by.site
  get.delta <- function(x) {
    return(x-min(x))
  }
  weights.integrating.phi <- matrix(nrow=dim(AIC.site.information)[1], ncol=dim(AIC.site.information)[2])
  for (site.number in sequence(dim(AIC.site.information)[2])) { #fine, I give up on the apply across three dimensions. BCO.
    local.matrix <- exp(-0.5*get.delta(AIC.site.information[,site.number,]))
    local.matrix <- local.matrix / sum(local.matrix)
    weights.integrating.phi[,site.number] <- rowSums(local.matrix)
  }
  rownames(weights.integrating.phi) <- dimnames(AIC.site.information)[[1]]

  reverse.weighted.mean <- function(w, x) {
    return(stats::weighted.mean(x, w))
  }

  average.c <- apply(weights.integrating.phi, 2, reverse.weighted.mean, x=aa.properties.reordered[,"c"])
  average.p <- apply(weights.integrating.phi, 2, reverse.weighted.mean, x=aa.properties.reordered[,"p"])
  average.v <- apply(weights.integrating.phi, 2, reverse.weighted.mean, x=aa.properties.reordered[,"v"])
  sliding.c <- zoo::rollmean(average.c, k=mean.width)
  sliding.p <- zoo::rollmean(average.p, k=mean.width)
  sliding.v <- zoo::rollmean(average.v, k=mean.width)
  #average.phi <- apply()
  par(mfcol=c(1,4))
  plot(x=sequence(length(average.c)), y=average.c, main="Composition", xlab="Site", pch=20,ylab="", bty="n", col=rgb(0,0,0,.5))
  lines(x=sequence(length(sliding.c)), y=sliding.c, lwd=2)
  plot(x=sequence(length(average.p)), y=average.p, main="Polarity", xlab="Site", pch=20, ylab="", bty="n", col=rgb(0,0,0,.5))
  lines(x=sequence(length(sliding.p)), y=sliding.p, lwd=2)
  plot(x=sequence(length(average.v)), y=average.v, main="Molecular volume", xlab="Site", pch=20, ylab="", bty="n", col=rgb(0,0,0,.5))
  lines(x=sequence(length(sliding.v)), y=sliding.v, lwd=2)
}

ComputeMutationFitnessesUnderGammaRates <- function(nuc.model="JC", base.freqs=rep(0.25, 4), nsites=1, C=4, Phi=0.5, q=4e-7, Ne=5e6, alpha=1.83, beta=0.10, gamma=0.0003990333, include.stop.codon=TRUE, numcode=1, diploid=TRUE, flee.stop.codon.rate=0.9999999, shape.gamma=1, n.pulls=1000) {
  Phi.vector <- Phi*rgamma(n.pulls, shape=shape.gamma, rate=shape.gamma)
  ComputeMutationFitnessesPhiFirst <- function(Phi=0.5, nuc.model="JC", base.freqs=rep(0.25, 4), nsites=1, C=4, q=4e-7, Ne=5e6, alpha=1.83, beta=0.10, gamma=0.0003990333, include.stop.codon=TRUE, numcode=1, diploid=TRUE, flee.stop.codon.rate=0.9999999) {
    return(ComputeMutationFitnesses(nuc.model=nuc.model, base.freqs=base.freqs, nsites=nsites, C=C, Phi=Phi, q=q, Ne=Ne, alpha=alpha, beta=beta, gamma=gamma, include.stop.codon=include.stop.codon, numcode=numcode, diploid=diploid, flee.stop.codon.rate=flee.stop.codon.rate))
  }
  results <- lapply(Phi.vector, ComputeMutationFitnessesPhiFirst, nuc.model=nuc.model, base.freqs=base.freqs, nsites=nsites, C=C, q=q, Ne=Ne, alpha=alpha, beta=beta, gamma=gamma, include.stop.codon=include.stop.codon, numcode=numcode, diploid=diploid, flee.stop.codon.rate=flee.stop.codon.rate)
  return(results)
}

# PlotCDFOfMutations <- function(mutation.fitness.object.list, values, optimal.aa=NULL, palette="Set1", lwd=2, ...) {
#   colors <- add.alpha(RColorBrewer::brewer.pal(length(mutation.fitness.object.list),palette),0.5)
#   results.to.plot <- lapply(mutation.fitness.object.list, LineMutationFitnessSpectra, optimal.aa=optimal.aa)
# }

# TODO
#
# Y axis: fixation probability relative to neutral
#
# Do for different amino acids
#
# Phi*g
#
# X axis is log(W)*Ne
#
# Lines are cdfs for diff amino acids
#
# Could have theta curve also on the plot
#
# When have gamma, pull from distribution, not just the category midpoints
