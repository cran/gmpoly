#' @title Arithmetic operators for multivariate polynomials
#'
#' @param e1,e2 for an unary operator, only \code{e1} must be given, a 
#'   \code{\link{gmpoly}} object; for a binary operator, at least one of 
#'   \code{e1} and \code{e2} must be a \code{\link{gmpoly}} object, and the 
#'   other must a \code{\link{gmpoly}} object as well or a scalar; the power
#'   operator (\code{^}) is an exception: one can only raise a 
#'   \code{\link{gmpoly}} object to a positive integer power
#'
#' @return A \code{\link{gmpoly}} object.
#' @export
#'
#' @examples library(gmpoly)
#' pol <- gmpoly("4 x^(2, 1, 1) + 1/2 x^(0,1,0)")
#' +pol
#' -pol
#' 2 * pol
#' pol / 2
#' pol + 5
#' pol - 5
#' pol^2
#' pol1 <- gmpoly("2 x^(1,1) - 5/3 x^(0,1)")
#' pol2 <- gmpoly("-2 x^(1,1) + 3 x^(2,1)")
#' pol1 + pol2
#' pol1 * pol2
#' pol1 == pol2
#' pol1 != pol2
Ops.gmpoly <- function(e1, e2 = NULL) {
  unary <- nargs() == 1L
  lclass <- .Method[1L] == "Ops.gmpoly" 
  rclass <- !unary && (.Method[2L] == "Ops.gmpoly")
  
  if(unary){
    if(.Generic == "+"){
      return(e1)
    }else if(.Generic == "-") {
      return(gmpoly_negate(e1))
    }else{
      stop("Unary operator '", .Generic, "' is not implemented for gmpolys.")
    }
  }
  
  if(!is.element(.Generic, c("+", "-", "*", "/", "^", "==", "!="))){
    stop("Operator '", .Generic, "' is not implemented for gmpolys.")
  }
  
  if(.Generic == "*"){
    if(lclass && rclass) {
      return(polynomialMul(e1, e2))
    }else if(lclass){
      return(gmpoly_times_scalar(e1, e2))
    }else if(rclass){
      return(gmpoly_times_scalar(e2, e1))
    }else{
      stop()
    }
  }else if(.Generic == "+"){
    if(lclass && rclass){
      if(gmpoly_eq_gmpoly(e1, gmpoly_negate(e2))){
        return(zeroPol(e1[["m"]]))
      }
      return(polynomialAdd(e1, e2))
    }else if(lclass){
      return(gmpoly_plus_scalar(e1, e2))
    }else if(rclass){
      return(gmpoly_plus_scalar(e2, e1))
    } else {
      stop()
    }
  }else if(.Generic == "-"){
    if(lclass && rclass){
      if(gmpoly_eq_gmpoly(e1, e2)){
        return(zeroPol(e1[["m"]]))
      }
      return(polynomialAdd(e1, gmpoly_negate(e2)))
    }else if(lclass){
      return(gmpoly_plus_scalar(e1, -e2))
    }else if(rclass){
      return(gmpoly_plus_scalar(gmpoly_negate(e2), e1))
    } else {
      stop()
    }
  }else if(.Generic == "^"){
    if(lclass && !rclass){
      return(gmpoly_power(e1, e2))
    }else{
      stop("Generic '^' not implemented in this case.")
    }
  } else if(.Generic == "==") {
    if(lclass && rclass){
      return(gmpoly_eq_gmpoly(e1, e2))
    }else{
      stop("Generic '==' only compares two `gmpoly` objects with one another.")
    }
  }else if(.Generic == "!="){
    if(lclass && rclass){
      return(!gmpoly_eq_gmpoly(e1, e2))
    }else{
      stop("Generic '!=' only compares two `gmpoly` objects with one another.")
    }
  }else if(.Generic == "/"){
    if(lclass && !rclass){
      if(e2 == 0){
        stop("Division by zero.")
      }
      return(gmpoly_times_scalar(e1, 1/e2))
    }else{
      stop("Generic '/' is only used to divide a `gmpoly` by a scalar.")
    }
  }
}

gmpoly_negate <- function(pol){
  if(isZeroPol(pol)){
    pol
  }else{
    pol[["coeffs"]] <- -pol[["coeffs"]]
    pol
  }
}

gmpoly_times_scalar <- function(pol, lambda){
  if(lambda == 0){
    return(zeroPol(pol[["m"]]))
  }
  if(lambda == 1 || isZeroPol(pol)){
    return(pol)
  }
  pol[["coeffs"]] <- lambda * pol[["coeffs"]]
  pol
}

# `mvp_plus_mvp` <- function(S1,S2){
#   if(is.zero(S1)){
#     return(S2)
#   } else if(is.zero(S2)){
#     return(S1)
#   } else {
#     jj <- mvp_add(
#       allnames1=S1[[1]],allpowers1=S1[[2]],coefficients1=S1[[3]],
#       allnames2=S2[[1]],allpowers2=S2[[2]],coefficients2=S2[[3]]
#     )
#     return(mvp(jj[[1]],jj[[2]],jj[[3]]))
#   }
# }

#' @importFrom gmp as.bigq
#' @noRd
gmpoly_plus_scalar <- function(pol, x){
  if(x == 0){
    return(pol)
  }
  m <- pol[["m"]]
  scalarPol <- gmpoly(coeffs = as.bigq(x), powers = rbind(rep(0L, m)))
  polynomialAdd(pol, scalarPol)
}

gmpoly_power <- function(pol, n){
  stopifnot(isPositiveInteger(n))
  if(n == 0){
    gmpoly(sprintf("x^(%s)", toString(rep("0", pol[["m"]]))))
  }else{
    Reduce(polynomialMul, rep(list(pol), n))
  }
}

gmpoly_eq_gmpoly <- function(pol1, pol2){
  pol1[["m"]] == pol2[["m"]] && all(pol1[["coeffs"]] == pol2[["coeffs"]]) &&
    all(pol1[["exponents"]] == pol2[["exponents"]])
}
