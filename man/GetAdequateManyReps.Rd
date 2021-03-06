% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelAdequacy.R
\name{GetAdequateManyReps}
\alias{GetAdequateManyReps}
\title{Parallel model adequacy test}
\usage{
GetAdequateManyReps(nreps, n.cores, model.to.reconstruct.under = "selac",
  model.to.simulate.under = "gtr", selac.obj.to.reconstruct,
  selac.obj.to.simulate, aa.optim.input = NULL,
  fasta.rows.to.keep = NULL, taxon.to.drop = 2,
  partition.number = 17, numcode = 1, for.gtr.only = NULL)
}
\arguments{
\item{nreps}{Specifies the number of repeated model adequact simulations.}

\item{n.cores}{Specifies the number of cores you want to use.}

\item{model.to.reconstruct.under}{Specifies the model that the internal nodes are to be reconstructed under assuming a single tip is pruned from the tree.}

\item{model.to.simulate.under}{Specifies the model that the simulation will be conducted along the pruned tip.}

\item{selac.obj.to.reconstruct}{The selac output object that contains the model parameters to be used in the reconstruction.}

\item{selac.obj.to.simulate}{The selac output object that contains the model parameters to be used in the simulation.}

\item{aa.optim.input}{A list of optimal amino acids with each list element designating a character vector for each gene. The optimal amino acids be the MLE from a selac run (default) or a list of user defined optimal A.A.}

\item{fasta.rows.to.keep}{Indicates which rows to remove in the input fasta files.}

\item{taxon.to.drop}{Specifies the tip based on the number in the phy object to be removed and simulated.}

\item{partition.number}{Specifies the partition number to conduct the model adequacy test.}

\item{numcode}{The ncbi genetic code number for translation. By default the standard (numcode=1) genetic code is used.}

\item{for.gtr.only}{A selac object that can be used as the reference optimal AA for when the adequacy of a GTR+G model is tested only.}
}
\description{
Performs model adequacy test using multiple cores
}
\details{
Performs a parallelized analysis of the model adequacy test. The test prunes out a user-specified taxon from the tree, performs site data reconstruction for all nodes in the tree under a user-specified model, then simulates the expected data of the pruned taxon according to a user-specified model along uniformly sampled points along the branch. The functionality of the reconstructed sequence is also calculated along the way to see how functionality changes as the simulation reaches the end of the known branch length. The output is a list with elements equally the number of repititions. Each element contains the functionality of the simulated points along equally spaced sampling points along the known branch length (i.e., edge.length * seq(0, 1, by=0.05))
}
