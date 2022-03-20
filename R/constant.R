#' @title Constant multivariate polynomial
#' @description Constructs a constant multivariate polynomial.
#'
#' @param m number of variables, a strictly positive integer
#' @param value the constant value of the polynomial; the 
#'   \code{\link[gmp]{as.bigq}} function is applied to this argument, so it can 
#'   be e.g. an integer or a character string such as \code{"2/3"} (avoid 
#'   decimal numbers)
#'
#' @return A \code{\link{gmpoly}} object.
#' @export
#' 
#' @importFrom gmp as.bigq
#'
#' @examples library(gmpoly)
#' gmpolyConstant(3, "2/3")
gmpolyConstant <- function(m, value){
  stopifnot(isStrictlyPositiveInteger(m))
  stopifnot(isScalar(value))
  value <- as.bigq(value)
  if(value == 0L){
    return(zeroPol(m))
  }
  gmpoly(coeffs = value, powers = rbind(rep(0L, m)))
}
