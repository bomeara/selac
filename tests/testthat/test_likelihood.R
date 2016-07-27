test_that("GTR_likelihood", {
	tree <- read.tree("rokasYeast.tre")
	phy <- drop.tip(tree, "Calb")
	yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
	yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
	chars <- DNAbinToNucleotideNumeric(yeast.gene)
	codon.data <- chars[phy$tip.label,]
	nuc.sites = SitePattern(codon.data)
	gtr_only = GetLikelihoodNucleotideForManyCharGivenAllParams(log(c(1,1,1,1,1)), nuc.sites, phy, nuc.model="GTR", include.gamma=FALSE, logspace=TRUE, ncats=4, verbose=FALSE, parallel.type="by.gene", n.cores=NULL)
	comparison <- identical(round(gtr_only,3), -8647.239)
	expect_true(comparison)
})


test_that("GTR+GAMMA_likelihood", {
	tree <- read.tree("rokasYeast.tre")
	phy <- drop.tip(tree, "Calb")
	yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
	yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
	chars <- DNAbinToNucleotideNumeric(yeast.gene)
	codon.data <- chars[phy$tip.label,]
	nuc.sites = SitePattern(codon.data)
	gtr_gamma <- GetLikelihoodNucleotideForManyCharGivenAllParams(log(c(.5,1,1,1,1,1)), nuc.sites, phy, nuc.model="GTR", include.gamma=TRUE, logspace=TRUE, ncats=4, verbose=FALSE, parallel.type="by.gene", n.cores=NULL)
	comparison <- identical(round(gtr_gamma,3), -8192.526)
	expect_true(comparison)
})


test_that("GY94_likelihood", {
    tree <- read.tree("rokasYeast.tre")
    phy <- drop.tip(tree, "Calb")
    yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
    yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
    chars <- DNAbinToCodonNumeric(yeast.gene)
    codon.data <- chars[phy$tip.label,]
    codon.data = SitePattern(codon.data)
    gy94 <- GetLikelihoodGY94_CodonForManyCharGivenAllParams(log(c(1,1)), codon.data, phy, numcode=1, logspace=TRUE, verbose=FALSE, parallel.type="by.gene", n.cores=NULL)
    comparison <- identical(round(gy94,3), -7826.123)
    expect_true(comparison)
})


test_that("FMutSel0_likelihood", {
    tree <- read.tree("rokasYeast.tre")
    phy <- drop.tip(tree, "Calb")
    yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
    yeast.gene <- as.list(as.matrix(cbind(yeast.gene))[1:7,])
    chars <- selac:::DNAbinToCodonNumeric(yeast.gene)
    codon.data <- chars[phy$tip.label,]
    codon.data = selac:::SitePattern(codon.data)
    fmutsel0 <- selac:::GetLikelihoodMutSel_CodonForManyCharGivenAllParams(log(c(.25,.25,.25, rep(1,5), .4, rep(1,19))), codon.data, phy, numcode=1,  nuc.model="GTR", logspace=TRUE, verbose=FALSE, parallel.type="by.gene", n.cores=NULL)
    comparison <- identical(round(fmutsel0,3), -8641.113)
    expect_true(comparison)
})


test_that("selac_likelihood", {
	set.seed(4)
	tree <- read.tree("rokasYeast.tre")
	phy <- drop.tip(tree, "Calb")
	yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
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
    selac <- GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, .25, .25, .25, rep(1,5))), codon.data, phy, aa.optim_array=aa.optim, root.p_array=empirical.codon.freq, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="GTR", codon.index.matrix, include.gamma=FALSE, ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, parallel.type="by.gene", n.cores=NULL)
    comparison <- identical(round(selac, 3), -6202.519)
	expect_true(comparison)
})

test_that("selac+GAMMA_likelihood", {
	set.seed(4)
	tree <- read.tree("rokasYeast.tre")
	phy <- drop.tip(tree, "Calb")
	yeast.gene <- read.dna("gene1Yeast.fasta", format="fasta")
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
	selac_gamma <- GetLikelihoodSAC_CodonForManyCharGivenAllParams(log(c(4*4e-7*.5*5e6, 1.829272, 0.101799, .25, .25, .25, rep(1,5), 5)), codon.data, phy, aa.optim_array=aa.optim, root.p_array=empirical.codon.freq, numcode=1, diploid=TRUE, aa.properties=NULL, volume.fixed.value=0.0003990333, nuc.model="GTR", codon.index.matrix, include.gamma=TRUE, ncats=4, k.levels=0, logspace=TRUE, verbose=FALSE, parallel.type="by.gene", n.cores=NULL)
	comparison <- identical(round(selac_gamma, 3), -6131.114)
	expect_true(comparison)
})

test_that("dealing_with_missing_data", {
    set.seed(4)
    phy <- rcoal(20)
    
    #part 1 -- pruning the taxa straightup:
    phy.pruned <- drop.tip(phy, c("t16", "t13"))
    phy.sort <- reorder(phy.pruned, "pruningwise")
    anc.indices <- unique(phy.sort$edge[,1])

    traits <- data.frame(taxon=phy$tip.label, trait=rep(1, length(phy$tip.label)))
    traits[1:7,2] = 2
    Q <- matrix(1,4,4)
    diag(Q) = 0
    diag(Q) = -rowSums(Q)
    scale.factor=1
    expList <- GetExpQt(phy=phy.pruned, Q=Q, scale.factor=scale.factor)
    nb.tip<-length(phy.pruned$tip.label)
    nb.node <- phy.pruned$Nnode
    nl <- nrow(Q)
    liks <- matrix(0, nb.tip + nb.node, nl)
    for(i in 1:nb.tip){
        if(!is.na(traits[i,1+1])){
            liks[i,traits[i,1+1]] <- 1
        }else{
            liks[i,] <- 1
        }
    }
    pruned.ll <- FinishLikelihoodCalculation(phy=phy.sort, liks=liks, Q=expList, root.p=rep(.25,4), anc=anc.indices)
    
    #part 2 -- Making the taxa uncertain in their scoring:
    traits <- data.frame(taxon=phy$tip.label, trait=rep(1, length(phy$tip.label)))
    traits[1:7,2] = 2
    traits[c(11,20), 2] = NA
    Q <- matrix(1,4,4)
    diag(Q) = 0
    diag(Q) = -rowSums(Q)
    scale.factor=1
    expList <- GetExpQt(phy=phy, Q=Q, scale.factor=scale.factor)
    nb.tip<-length(phy$tip.label)
    nb.node <- phy$Nnode
    nl <- nrow(Q)
    liks <- matrix(0, nb.tip + nb.node, nl)
    for(i in 1:nb.tip){
        if(!is.na(traits[i,1+1])){
            liks[i,traits[i,1+1]] <- 1
        }else{
            liks[i,] <- 1
        }
    }
    phy.sort <- reorder(phy, "pruningwise")
    anc.indices <- unique(phy.sort$edge[,1])
    indicator.ll <- FinishLikelihoodCalculation(phy=phy, liks=liks, Q=expList, root.p=rep(.25,4), anc=anc.indices)
    comparison <- identical(pruned.ll, indicator.ll)
    expect_true(comparison)
})



