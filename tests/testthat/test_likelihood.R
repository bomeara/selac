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

