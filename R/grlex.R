grlexRank <- function(a){
  a <- as.integer(a)
  m <- length(a)
  nm <- sum(a)
  ns <- nm + m - 1L
  ks <- m - 1L
  xs <- integer(ks)
  xs[1L] <- a[1L] + 1L
  for(i in seq_len(m)[-c(1L, 2L)]){
    xs[i-1L] <- xs[i-2L] + a[i-1L] + 1L
  }
  rank <- 1L
  for(i in seq_len(ks)){
    if(i == 1L){
      tim1 <- 0L
    }
    else{
      tim1 <- xs[i-1L]
    }
    if(tim1 + 1L <= xs[i] - 1L){
      for(j in (tim1 + 1L):(xs[i] - 1L)){
        rank <- rank + as.integer(choose(ns - j, ks - i))
      }
    }
  }
  for(n in seq_len(nm)-1L){
    rank <- rank + as.integer(choose(n + m - 1L, n))
  }
  rank
}


grlexUnrank <- function(m, rank){
  rank <- as.integer(rank)
  if(m == 1L) {
    return(rank - 1L)
  }
  m <- as.integer(m)
  rank1 <- 1L
  nm <- -1L
  while(TRUE){
    nm <- nm + 1L
    r <- as.integer(choose(nm + m - 1L, nm))
    if(rank < rank1 + r){
      break
    }
    rank1 <- rank1 + r
  }
  rank2 <- rank - rank1
  ks <- m - 1L
  ns <- nm + m - 1L
  xs <- integer(ks)
  j <- 1L
  for(i in seq_len(ks)){
    r <- as.integer(choose(ns - j, ks - i))
    
    while(r <= rank2 && 0L < r){
      rank2 <- rank2 - r
      j <- j + 1L
      r <- as.integer(choose(ns - j, ks - i))
    }
    xs[i] <- j
    j <- j + 1L
  }
  a <- integer(m)
  a[1L] <- xs[1L] - 1L
  for(i in seq_len(m)[-c(1L, 2L)]){
    a[i - 1L] <- xs[i - 1L] - xs[i - 2L] - 1L
  }
  a[m] <- ns - xs[ks]
  a
}
