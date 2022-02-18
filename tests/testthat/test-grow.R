library(gmp)

test_that("Growing preserves the sorting", {
  set.seed(666L)
  pol <- gmpoly(
    coeffs = as.bigq(rep(1L, 50L)),
    powers = matrix(rpois(50L*3L, 8), nrow = 50L, ncol = 3L)
  )
  gpol <- gmpolyGrow(pol)
  expect_false(is.unsorted(gpol[["exponents"]]))
})