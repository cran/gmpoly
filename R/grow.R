#' @title Grow polynomial
#' @description Grow a multivariate polynomial by including a new variable.
#'
#' @param pol a \code{\link{gmpoly}} object
#'
#' @return The "same" multivariate polynomial as \code{pol}, except that it 
#'   has an additional variable.
#' @export
#'
#' @examples library(gmpoly)
#' pol <- gmpoly("3 x^(1,2) - 1/7 x^(5,3)")
#' gmpolyGrow(pol)
gmpolyGrow <- function(pol){
  stopifnot(inherits(pol, "gmpoly"))
  m <- pol[["m"]]
  powers <- attr(pol, "powers")
  if(is.null(powers)){
    powers <- t(vapply(pol[["exponents"]], function(e){
      grlexUnrank(m, e)
    }, integer(m)))
  }
  powers <- cbind(powers, 0L)
  newpol <- list(
    "coeffs" = pol[["coeffs"]], 
    "exponents" = apply(powers, 1L, grlexRank),
    "m" = m + 1L
  )
  attr(newpol, "powers") <- powers
  class(newpol) <- "gmpoly"
  newpol
}