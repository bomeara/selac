test_that("GetLikelihoodSAC_CodonForSingleCharGivenOptimum", {
  ref <- c("U15717", "U15718", "U15719", "U15720",
         "U15721", "U15722", "U15723", "U15724")
  library(ape)
  Rampho <- read.GenBank(ref)
  chars <- DNAbinToCodonNumeric(Rampho)
  tree <- rcoal(n=dim(chars)[1], tip.label=chars[,1])
  Q_codon <- CreateCodonFixationRateMatrix(aa_op="M", s=2, aa.distances=CreateAADistanceMatrix())
  likelihood <- GetLikelihoodSAC_CodonForSingleCharGivenOptimum(codon.data=chars, phy=tree, charnum=4, Q_codon=Q_codon)
  expect_is(likelihood, "numeric")
})
