test_that("Growing works", {
  set.seed(666L)
  pol <- rgmpol(nvars = 10L, nterms = 20L)
  gpol <- gmpolyGrow(pol)
  nums <- sample.int(7L, size = 10L, replace = TRUE) - 1L
  dens <- sample.int(6L, size = 10L, replace = TRUE)
  x <- gmp::as.bigq(nums, dens)
  xplus <- c(x, gmp::as.bigq("1"))
  expect_identical(gmpolyEval(pol, x), gmpolyEval(gpol, xplus))
})

test_that("Growing preserves the sorting", {
  set.seed(666L)
  pol <- rgmpol(nvars = 3L, nterms = 50L)
  gpol <- gmpolyGrow(pol)
  expect_identical(polynomialSort(gpol), unclass(gpol))
})