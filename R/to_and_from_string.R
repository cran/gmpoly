polAsString <- function(pol, powers){
  coeffs <- pol[["coeffs"]]
  if(ncol(powers) != 1L){
    powers <- apply(powers, 1L, paste0, collapse = ",")    
  }
  terms <- paste0(coeffs, " x^(", powers, ")")
  s <- paste0(terms, collapse = " + ")
  s <- gsub("+ -", "- ", s, fixed = TRUE)
  s <- gsub(" 1 x", " x", s, fixed = TRUE)
  sub("^1 ", "", s)
}

#' @importFrom purrr transpose
#' @importFrom gmp as.bigq
#' @noRd
stringToPol <- function(p){
  p <- gsub("\\)\\s*-\\s*(\\d*/*\\d*)\\s*", ")+-\\1", p)#, perl = TRUE)
  p <- gsub("^-\\s*x", "-1x", trimws(p, "left"))
  terms <- strsplit(p, "+", fixed = TRUE)[[1L]]
  csts <- !grepl("x", terms)
  terms[csts] <- paste0(terms[csts], "x^(0")
  ss <- transpose(strsplit(terms, "x^(", fixed = TRUE))#, .names = c("coeff", "power"))
  coeffs <- as.bigq(unlist(ss[[1L]], recursive = FALSE))
  coeffs[is.na(coeffs)] <- as.bigq(1L)
  powers <- sub(")", "", unlist(ss[[2L]], recursive = FALSE), fixed = TRUE)
  powers <- lapply(strsplit(powers, ","), as.integer)
  positive <- vapply(powers, function(a){
    all(a >= 0)
  }, logical(1L))
  if(!all(positive)){
    stop("Negative powers are not allowed.")
  }
  i <- 1L
  # m <- length(powers[[1L]])
  # nterms <- length(powers)
  # while(m == 1L && i < nterms){
  #   i <- i + 1L
  #   m <- length(powers[[i]])
  # }
  powers <- do.call(rbind, powers)
  pol <- polynomialSort(list(
    "coeffs" = coeffs, 
    "powers" = powers
  ))
  if(anyDuplicated(powers) || any(coeffs == 0)){
    pol <- polynomialCompress(pol)
  }
  if(all(pol[["coeffs"]] == 0L)){
    return(zeroPol(ncol(powers)))
  }
  class(pol) <- "gmpoly"
  pol
}
