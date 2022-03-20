rgmpol <- function(nvars, nterms){
  powers <- matrix(
    rpois(nterms*nvars, 4), 
    nrow = nterms, ncol = nvars
  )
  nums <- sample.int(7L, size = nterms, replace = TRUE) - 1L
  dens <- sample.int(6L, size = nterms, replace = TRUE)
  coeffs <- gmp::as.bigq(nums, dens)
  gmpoly(coeffs = coeffs, powers = powers)
}