#' @importFrom gmp as.bigq
#' @noRd
zeroPol <- function(m){
  pol <- list(
    "coeffs" = as.bigq(0L),
    "powers" = t(rep(0L, m))
  )
  attr(pol, "zero") <- TRUE
  class(pol) <- "gmpoly"
  pol
}

isZeroPol <- function(pol){
  isTRUE(attr(pol, "zero"))
}

#' @title Define a multivariate polynomial
#' @description Define a multivariate polynomial from a string or from its 
#'   coefficients and powers. 
#'
#' @param string a string such as \code{"x^(1,2,1) + 4 x^(0,2,3)"}, or you can 
#'   define the polynomial with \code{coeffs} and \code{powers}; in this case
#'   set \code{string} to \code{NULL} or to nothing (i.e. missing)
#' @param coeffs the coefficients of the polynomial, in case you don't define 
#'   it with the \code{string} argument; this must be a vector of 
#'   \code{\link[gmp]{bigq}} numbers
#' @param powers the powers of the polynomial, in case you don't define it 
#'   with the \code{string} argument; this must be a matrix of integers, one 
#'   row for each term (hence \code{nrow(powers)} must equal 
#'   \code{length(coeffs)}), and the number of columns is the number of 
#'   variables
#'
#' @return A \code{gmpoly} object.
#' @export
#' 
#' @importFrom gmp is.bigq
#'
#' @examples library(gmpoly)
#' gmpoly("5/2 x^(2,3) + 3 x^(1,1)")
#' gmpoly("5/2 x^(2,3) - 3 x^(1,1)")
#' gmpoly("-x^(1,2,3) + 4/7 x^(3,1,1)")
#' library(gmp)
#' gmpoly(
#'   coeffs = as.bigq(c(5, 7), c(8, 9)), 
#'   powers = rbind(c(1, 1, 1), c(2, 2, 2))
#' )
gmpoly <- function(string, coeffs = NULL, powers = NULL){
  if(missing(string) || is.null(string)){
    stopifnot(is.bigq(coeffs))
    stopifnot(is.matrix(powers))
    stopifnot(length(coeffs) == nrow(powers))
    m <- ncol(powers)
    if(all(coeffs == 0L)){
      return(zeroPol(m))
    }
    storage.mode(powers) <- "integer"
    stopifnot(all(powers >= 0L))
    pol <- polynomialSort(list(
      "coeffs" = coeffs,
      "powers" = powers
    ))
    if(anyDuplicated(powers) || any(coeffs == 0)){
      pol <- polynomialCompress(pol)
    }
  }else{
    stopifnot(isString(string))
    pol <- stringToPol(string)
    if(all(pol[["coeffs"]] == 0L)){
      return(zeroPol(ncol(pol[["powers"]])))
    }
  }
  class(pol) <- "gmpoly"
  pol
}

#' @title Print a multivariate polynomial
#' @description Print a multivariate polynomial of class \code{gmpoly}.
#'
#' @param x a \code{\link{gmpoly}} object
#' @param ... ignored
#'
#' @return No value, just prints the polynomial.
#' @export
print.gmpoly <- function(x, ...){
  cat("gmpoly object algebraically equal to\n")
  cat(polAsString(x, x[["powers"]]))
  cat("\n")
}


#' @title Conversion to 'mvp' polynomial
#' @description Convert a \code{\link{gmpoly}} polynomial to a 
#'   \code{\link[mvp:mvp-package]{mvp}} polynomial.
#'
#' @param pol a \code{\link{gmpoly}} object
#'
#' @return A \code{\link[mvp:mvp-package]{mvp}} object.
#' @export
#' 
#' @importFrom gmp asNumeric
#' @importFrom mvp mvp
#'
#' @examples library(gmpoly)
#' pol <- gmpoly("5/2 x^(2,2,3) + 3 x^(1,0,1)")
#' gmpoly2mvp(pol)
gmpoly2mvp <- function(pol){
  powers <- pol[["powers"]]
  nterms <- nrow(powers)
  mvp(
    vars = rep(list(paste0("x_", 1L:ncol(powers))), nterms),
    powers = lapply(seq_len(nterms), function(i) powers[i, ]),
    coeffs = asNumeric(pol[["coeffs"]])
  )
}
