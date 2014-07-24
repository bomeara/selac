test_that("CompareVectors", {
a<-c(1, 2, 3, 4)
b<-c(1, 3, 3, 5)
result<-CompareVectors(a, b)
expect_equal(result$num, 2)
expect_equal(result$pos, c(2,4))
})

test_that("TranslateCodon", {
	codon.string="atg"
	numcode=1
	expect_equal("M", TranslateCodon(codon.string, numcode))
})

test_that("CreateAADistanceMatrix", {
	aa.distances <- CreateAADistanceMatrix()
	expect_equal(20, dim(aa.distances)[2])
	expect_equal(aa.distances[17,4], aa.distances[4,17]) #symmetry for distances
})

test_that("CreateNucleotideMutationMatrix", {
	nuc.rates.jc <- CreateNucleotideMutationMatrix()
	expect_equal(4, dim(nuc.rates.jc)[2])
	expect_equal(2, length(unique(as.vector(nuc.rates.jc))))
	expect_error(CreateNucleotideMutationMatrix(model="GTR"))
})

test_that("CreateCodonMutationMatrix", {
	codon.rates <- CreateCodonMutationMatrix(CreateNucleotideMutationMatrix())
	expect_equal(64, dim(codon.rates)[2])
	expect_equal(0, rowSums(codon.rates[7,]))
})

test_that("CreateCodonFixationRateMatrix", {
	aa.distances <- CreateAADistanceMatrix()
	aa_op <- "M"
	s <- 0.02
})