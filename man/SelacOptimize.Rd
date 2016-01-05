% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selac.R
\name{SelacOptimize}
\alias{SelacOptimize}
\title{Optimize parameters under the SELAC model}
\usage{
SelacOptimize(codon.data.path, n.partitions = NULL, phy,
  edge.length = "optimize", edge.linked = TRUE, optimal.aa = "optimize",
  nuc.model = "GTR", include.gamma = FALSE, ncats = 4, numcode = 1,
  diploid = TRUE, k.levels = 0, aa.properties = NULL, verbose = FALSE,
  parallel.type = "by.gene", n.cores = NULL,
  max.tol = .Machine$double.eps^0.25, fasta.rows.to.keep = NULL)
}
\arguments{
\item{codon.data.path}{Provides the path to the directory containing the gene specific fasta files of coding data.}

\item{n.partitions}{The number of partitions to analyze. The order is based on the Unix order of the fasta files in the directory.}

\item{phy}{The phylogenetic tree to optimize the model parameters.}

\item{edge.length}{A logical indicating whether or not edge lengths should be optimized.}

\item{edge.linked}{A logical indicating whether or not edge lengths should be optimized separately for each gene. By default, a single set of each lengths is optimized for all genes.}

\item{optimal.aa}{Indicates what type of optimal.aa should be used. There are three options: "none", "majrule", or "optimize".}

\item{nuc.model}{Indicates what type nucleotide model to use. There are three options: "JC", "GTR", or "UNREST".}

\item{include.gamma}{A logical indicating whether or not to include a discrete gamma model.}

\item{ncats}{The number of discrete categories.}

\item{numcode}{The The ncbi genetic code number for translation. By default the standard (numcode=1) genetic code is used.}

\item{diploid}{A logical indicating whether or not the organism is diploid or not.}

\item{k.levels}{Provides how many levels in the polynomial. By default we assume a single level (i.e., linear).}

\item{aa.properties}{User-supplied amino acid distance properties. By default we assume Grantham (1974) properties.}

\item{verbose}{Logical indicating whether each iteration be printed to the screen.}

\item{parallel.type}{Designates whether a parallell run should occur by gene ("by.gene") or by site ("by.site").}

\item{n.cores}{The number of cores to run the analyses over.}

\item{max.tol}{Supplies the relative optimization tolerance.}

\item{fasta.rows.to.keep}{Indicates which rows to remove in the input fasta files.}
}
\description{
Optimizes model parameters under the SELAC model, or the more traditional nucleotide based models
}
\details{
SELAC stands for SELection on Amino acids and/or Codons. This function takes a user supplied topology and a set of fasta formatted sequences and optimizes the parameters in the SELAC model. Selection is based on selection towards an optimal amino acid at each site. The optimal amino acid at a side could be assumed to be based on a majority rule (\code{optimal.aa="majrule"}), or actually optimized as part of the optimization routine (\code{optimal.aa="optimize"}. Note that by setting \code{optimal.aa="none"} reverts to the traditional nucleotide based model. Also of note, is that the presence of stop codons produces bad behavior. Be sure that these are removed from the data prior to running the model.
}
