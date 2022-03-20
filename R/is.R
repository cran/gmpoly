isPositiveInteger <- function(n){
  is.vector(n) && is.numeric(n) && length(n) == 1L && !is.na(n) && floor(n) == n
}

isStrictlyPositiveInteger <- function(n){
  isPositiveInteger(n) && n != 0
}

isScalar <- function(x){
  is.atomic(x) && length(x) == 1L && !is.na(x)
}

isString <- function(x){
  isScalar(x) && is.character(x)
}
