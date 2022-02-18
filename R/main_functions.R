polynomialSort <- function(pol){
  exponents <- pol[["exponents"]]
  idx <- order(exponents)
  spol <- list(
    "coeffs" = pol[["coeffs"]][idx],
    "exponents" = exponents[idx],
    "m" = pol[["m"]]
  )
  attr(spol, "powers") <- attr(pol, "powers")[idx, , drop = FALSE]
  spol
}


#' @importFrom gmp NA_bigq_
#' @noRd
polynomialCompress <- function(pol){
  coeffs <- pol[["coeffs"]]
  # zeros <- coeffs == 0L
  # if(any(zeros)){
  #   coeffs <- coeffs[!zeros]
  #   exponents <- pol[["exponents"]][!zeros]
  # }else{
  exponents <- pol[["exponents"]]
  # }
  o1 <- o2 <- length(exponents)
  exponents2 <- integer(o2)
  coeffs2 <- rep(NA_bigq_, o2)
  get <- put <- 0L
  while(get < o1){
    get <- get + 1L
    if(0L == put){
      put <- put + 1L
      coeffs2[put] <- coeffs[get]
      exponents2[put] <- exponents[get]
    }else{
      if(exponents2[put] == exponents[get]){
        coeffs2[put] <- coeffs2[put] + coeffs[get]
      }else{
        put <- put + 1L
        coeffs2[put] <- coeffs[get]
        exponents2[put] <- exponents[get]
      }
    }
  }
  o2 <- put
  get <- put <- 1L
  while(get <= o2){
    if(coeffs2[get] != 0L){
      coeffs2[put] <- coeffs2[get]
      exponents2[put] <- exponents2[get]
      put <- put + 1L
    }
    get <- get + 1L
  }
  o2 <- put - 1L
  list(
    "coeffs" = coeffs2[1L:o2],
    "exponents" = exponents2[1L:o2],
    "m" = pol[["m"]]
  )
}


polynomialAdd <- function(pol1, pol2){
  m <- pol1[["m"]]
  if(m != pol2[["m"]]){
    stop(
      "Adding polynomials is possible only if the number of variables is the ",
      "same for the two polynomials.",
      call. = TRUE
    )
  }
  if(isZeroPol(pol1)){
    return(pol2)
  }
  if(isZeroPol(pol2)){
    return(pol1)
  }
  coeffs <- c(pol1[["coeffs"]], pol2[["coeffs"]])
  exponents <- c(pol1[["exponents"]], pol2[["exponents"]])
  powers1 <- attr(pol1, "powers")
  powers2 <- attr(pol2, "powers")
  powers <- if(!is.null(powers1) && !is.null(powers2)){
    rbind(powers1, powers2)
  }
  pol <- list(
    "coeffs" = coeffs,
    "exponents" = exponents,
    "m" = m
  )
  attr(pol, "powers") <- powers
  spol <- polynomialSort(pol)
  notCompressed <- anyDuplicated(exponents) || any(coeffs == 0)
  if(notCompressed){
    spol <- polynomialCompress(spol)
  }
  class(spol) <- "gmpoly"
  spol
}

#' @importFrom gmp outer
#' @noRd
polynomialMul <- function(pol1, pol2){
  m <- pol1[["m"]]
  if(m != pol2[["m"]]){
    stop(
      "Multiplying polynomials is possible only if the number of variables is ",
      "the same for the two polynomials.",
      call. = TRUE
    )
  }
  if(isZeroPol(pol1) || isZeroPol(pol2)){
    return(zeroPol(m))
  }
  coeffs1 <- pol1[["coeffs"]]
  coeffs2 <- pol2[["coeffs"]]
  o1 <- length(coeffs1)
  o2 <- length(coeffs2)
  exponents1 <- pol1[["exponents"]]
  exponents2 <- pol2[["exponents"]]
  powers1 <- attr(pol1, "powers")
  powers2 <- attr(pol2, "powers")
  if(is.null(powers1)){
    powers1 <- t(vapply(exponents1, function(e){
      grlexUnrank(m, e)
    }, integer(m)))
  }
  if(is.null(powers2)){
    powers2 <- t(vapply(exponents2, function(e){
      grlexUnrank(m, e)
    }, integer(m)))
  }
  coeffs <- c(gmp::outer(coeffs1, coeffs2))
  nterms <- o1 * o2
  exponents <- integer(nterms)
  powers <- matrix(NA_integer_, nrow = nterms, ncol = m)
  o <- 1L
  for(j in 1L:o2){
    for(i in 1L:o1){
      f1 <- powers1[i, ]
      f2 <- powers2[j, ]
      f <- f1 + f2
      powers[o, ] <- f
      exponents[o] <- grlexRank(f)
      o <- o + 1L
    }
  }
  pol <- list(
    "coeffs" = coeffs,
    "exponents" = exponents,
    "m" = m
  )
  attr(pol, "powers") <- powers
  spol <- polynomialSort(pol)
  notCompressed <- anyDuplicated(exponents) || any(coeffs == 0)
  if(notCompressed){
    spol <- polynomialCompress(spol)
  }
  class(spol) <- "gmpoly"
  spol
}
