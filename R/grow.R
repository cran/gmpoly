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
  newpol <- list(
    "coeffs" = pol[["coeffs"]], 
    "powers" = cbind(pol[["powers"]], 0L)
  )
  class(newpol) <- "gmpoly"
  newpol
}