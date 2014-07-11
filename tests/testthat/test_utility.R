test_that("CompareVectors", {
a<-c(1, 2, 3, 4)
b<-c(1, 3, 3, 5)
result<-CompareVectors(a, b)
expect_equal(result$num, 2)
expect_equal(result$pos, c(2,4))
})

