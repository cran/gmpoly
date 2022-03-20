lexorder <- function(M){
  do.call(order, lapply(seq_len(ncol(M)), function(i) M[, i]))
}
