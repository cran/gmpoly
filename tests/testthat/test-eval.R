library(gmp, warn.conflicts = FALSE)

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

test_that("Evaluation addition", {
  set.seed(666L)
  pol1 <- rgmpol(nvars = 3L, nterms = 50L)
  pol2 <- rgmpol(nvars = 3L, nterms = 50L)
  X <- rbind(
    t(as.bigq(c(1L, 1L, 1L))),
    t(as.bigq(c(2L, -3L, 4L), c(5L, 4L, 3L)))
  )
  expect_equal(
    gmpolyEval(pol1 + pol2, X), 
    gmpolyEval(pol1, X) + gmpolyEval(pol2, X)
  )
  expect_equal(
    gmpolyEval(pol1 - pol2, X), 
    gmpolyEval(pol1, X) - gmpolyEval(pol2, X)
  )
})

test_that("Evaluation multiplication", {
  set.seed(666L)
  pol1 <- rgmpol(nvars = 3L, nterms = 10L)
  pol2 <- rgmpol(nvars = 3L, nterms = 10L)
  X <- rbind(
    t(as.bigq(c(1L, 1L, 1L))),
    t(as.bigq(c(2L, -3L, 4L), c(5L, 4L, 3L)))
  )
  expect_equal(
    gmpolyEval(pol1 * pol2, X), 
    gmpolyEval(pol1, X) * gmpolyEval(pol2, X)
  )
})
