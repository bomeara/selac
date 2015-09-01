test_that("GetLikelihoodSAC_CodonForSingleCharGivenOptimum", {
	ref <- c("U15717", "U15718", "U15719", "U15720",
			 "U15721", "U15722", "U15723", "U15724")
	library(ape)
	Rampho <- read.GenBank(ref)
	Rampho <- as.list(as.matrix(cbind(Rampho))[1:8, -1]) #drop first site so start at codon pos 1
	chars <- DNAbinToCodonNumeric(Rampho)
	tree <- rcoal(n=dim(chars)[1], tip.label=chars[,1])
	Q_codon <- CreateCodonFixationRateMatrix(aa_op="M", s=2, aa.distances=CreateAADistanceMatrix(), numcode=2) #vertebrate mt
	stop("Need to add codon substitution rate")
	likelihood <- GetLikelihoodSAC_CodonForSingleCharGivenOptimum(charnum=4, codon.data=chars, phy=tree, Q_codon=Q_codon)
	expect_is(likelihood, "numeric")
})

test_that("GetLikelihoodSAC_CodonForManyCharGivenFixedOptimumAndQAndRoot", {
	ref <- c("U15717", "U15718", "U15719", "U15720",
			 "U15721", "U15722", "U15723", "U15724")
	library(ape)
	Rampho <- read.GenBank(ref)
	Rampho <- as.list(as.matrix(cbind(Rampho))[1:8, -1]) #drop first site so start at codon pos 1
	chars <- DNAbinToCodonNumeric(Rampho)
	tree <- rcoal(n=dim(chars)[1], tip.label=chars[,1])
	Q_codon <- CreateCodonFixationRateMatrix(aa_op="M", s=2, aa.distances=CreateAADistanceMatrix(), numcode=2) #vertebrate mt
    stop("Need to add codon substitution rate")
	
	likelihood <- GetLikelihoodSAC_CodonForManyCharGivenFixedOptimumAndQAndRoot(codon.data=chars, phy=tree, Q_codon=Q_codon)
	expect_is(likelihood, "numeric")
})

test_that("EstimateParametersCodon", {
	ref <- c("U15717", "U15718", "U15719", "U15720",
			 "U15721", "U15722", "U15723", "U15724")
	library(ape)
	Rampho <- read.GenBank(ref)
	Rampho <- as.list(as.matrix(cbind(Rampho))[1:8, -1]) #drop first site so start at codon pos 1
	Rampho <- as.list(as.matrix(cbind(Rampho))[1:8, 1:36]) #for speed, look at first few aa only
	chars <- DNAbinToCodonNumeric(Rampho)
	tree <- rcoal(n=dim(chars)[1], tip.label=chars[,1])
	results <- EstimateParametersCodon(codon.data=chars, phy=tree, root="majrule", optimal.aa="majrule", numcode=2)
	expect_is(results$objective, "numeric")
})

test_that("GTR_likelihood", {
	tree <- read.tree("/Users/jbeaulieu/selac/data/rokasYeast.tre")
	phy <- drop.tip(tree, "Calb")
	yeast.gene <- read.dna("/Users/jbeaulieu/selac/data/gene1Yeast.fasta", format="fasta")
	yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
	chars <- DNAbinToNucleotideNumeric(yeast.gene)
	codon.data <- chars[phy$tip.label,]
	nuc.sites = SitePattern(codon.data)
	gtr_only = GetLikelihoodNucleotideForManyCharGivenAllParams(log(c(1,1,1,1,1)), nuc.sites, phy, nuc.model="GTR", include.gamma=FALSE, logspace=TRUE, ncats=4, verbose=FALSE)
	comparison <- identical(round(gtr_only,3), -8647.239)
	expect_true(comparison)
})

test_that("GTR+GAMMA_likelihood", {
	tree <- read.tree("../data/rokasYeast.tre")
	phy <- drop.tip(tree, "Calb")
	yeast.gene <- read.dna("../data/gene1Yeast.fasta", format="fasta")
	yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
	chars <- DNAbinToNucleotideNumeric(yeast.gene)
	codon.data <- chars[phy$tip.label,]
	nuc.sites = SitePattern(codon.data)
	gtr_gamma = GetLikelihoodNucleotideForManyCharGivenAllParams(log(c(.5,1,1,1,1,1)), nuc.sites, phy, nuc.model="GTR", include.gamma=TRUE, logspace=TRUE, ncats=4, verbose=FALSE)
	comparison <- identical(round(gtr_gamma,3), -8192.526)
	expect_true(comparison)
})

test_that("selac_likelihood", {
	set.seed(4)
	tree <- read.tree("../data/rokasYeast.tre")
	phy <- drop.tip(tree, "Calb")
	yeast.gene <- read.dna("../data/gene1Yeast.fasta", format="fasta")
	yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
	chars <- DNAbinToCodonNumeric(yeast.gene)
	codon.data <- chars[phy$tip.label,]
	aa.data <- ConvertCodonNumericDataToAAData(codon.data, numcode=1)
	aa.optim <- apply(aa.data[, -1], 2, GetMaxName) #starting values for all, final values for majrule
	aa.optim.full.list <- aa.optim
	empirical.codon.freq <- CodonEquilibriumFrequencies(codon.data[,-1], aa.optim, numcode=1)
	aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
	colnames(aa.optim.frame.to.add) <- colnames(codon.data)
	codon.data <- rbind(codon.data, aa.optim.frame.to.add)
	codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
	aa.optim = codon.data$optimal.aa
	codon.index.matrix = CreateCodonMutationMatrixIndex()
	selac <- GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(1.829272, 0.101799, 0.0003990333, 5e6, .25, .25, .25, rep(1,5))), codon.data, phy, aa.optim_array=aa.optim, root.p_array=empirical.codon.freq, numcode=1, aa.properties=NULL, nuc.model="GTR", codon.index.matrix, include.gamma=FALSE, ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE)
    comparison <- identical(round(selac, 3), -6102.21)
	expect_true(comparison)
})

test_that("selac+GAMMA_likelihood", {
	set.seed(4)
	tree <- read.tree("../data/rokasYeast.tre")
	phy <- drop.tip(tree, "Calb")
	yeast.gene <- read.dna("../data/gene1Yeast.fasta", format="fasta")
	yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
	chars <- DNAbinToCodonNumeric(yeast.gene)
	codon.data <- chars[phy$tip.label,]
	aa.data <- ConvertCodonNumericDataToAAData(codon.data, numcode=1)
	aa.optim <- apply(aa.data[, -1], 2, GetMaxName) #starting values for all, final values for majrule
	aa.optim.full.list <- aa.optim
	empirical.codon.freq <- CodonEquilibriumFrequencies(codon.data[,-1], aa.optim, numcode=1)
	aa.optim.frame.to.add <- matrix(c("optimal", aa.optim), 1, dim(codon.data)[2])
	colnames(aa.optim.frame.to.add) <- colnames(codon.data)
	codon.data <- rbind(codon.data, aa.optim.frame.to.add)
	codon.data <- SitePattern(codon.data, includes.optimal.aa=TRUE)
	aa.optim = codon.data$optimal.aa
	codon.index.matrix = CreateCodonMutationMatrixIndex()
	selac_gamma <- GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(1.829272, 0.101799, 0.0003990333, 5e6, .25, .25, .25, rep(1,5), 5)), codon.data, phy, aa.optim_array=aa.optim, root.p_array=empirical.codon.freq, numcode=1, aa.properties=NULL, nuc.model="GTR", codon.index.matrix, include.gamma=TRUE, ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE)
	comparison <- identical(round(selac_gamma, 3), -6097.791)
	expect_true(comparison)
})









