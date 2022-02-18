library(gmp)

test_that("One evaluation", {
  pol <- gmpoly("5/2 x^(2,3) + 3 x^(1,1)")
  expect_equal(gmpolyEval(pol, as.bigq(c(1, 1))), as.bigq(11, 2))
})

test_that("Two evaluations", {
  pol <- gmpoly("5/2 x^(2,3) + 3 x^(1,1)")
  X <- rbind(
    t(as.bigq(c(1, 1))),
    t(as.bigq(c(3, 4), c(4, 3)))
  )
  expect_equal(gmpolyEval(pol, X), c(as.bigq(11, 2), as.bigq(19, 3)))
})