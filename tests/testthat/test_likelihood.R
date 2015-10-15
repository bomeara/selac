test_that("GTR_likelihood", {
	tree <- read.tree("../data/rokasYeast.tre")
	phy <- drop.tip(tree, "Calb")
	yeast.gene <- read.dna("../data/gene1Yeast.fasta", format="fasta")
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
	gtr_gamma <- GetLikelihoodNucleotideForManyCharGivenAllParams(log(c(.5,1,1,1,1,1)), nuc.sites, phy, nuc.model="GTR", include.gamma=TRUE, logspace=TRUE, ncats=4, verbose=FALSE)
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
	selac <- GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5, 1.829272, 0.101799, 5e6, .25, .25, .25, rep(1,5))), codon.data, phy, aa.optim_array=aa.optim, root.p_array=empirical.codon.freq, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="GTR", codon.index.matrix, include.gamma=FALSE, ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE)
    comparison <- identical(round(selac, 3), -6202.519)
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
	selac_gamma <- GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5, 1.829272, 0.101799, 5e6, .25, .25, .25, rep(1,5), 5)), codon.data, phy, aa.optim_array=aa.optim, root.p_array=empirical.codon.freq, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="GTR", codon.index.matrix, include.gamma=TRUE, ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE)
	comparison <- identical(round(selac_gamma, 3), -6131.114)
	expect_true(comparison)
})









