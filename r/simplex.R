
check_simplex <- function(pzx) {
  
  if (min(pzx) < 0 || max(pzx) > 1) {
    return(TRUE)
  }
  FALSE
}

locate_simplex <- function(prx, Ainv, dls = FALSE) {
  
  b <- Ainv %*% prx
  a <- prx
  lambda <- 1
  cnt <- 0
  while(sum(b < 0)) {
    cnt <- cnt+1
    idx <- which(b < 0)
    lambda_max <- min(a[idx] / (a-b)[idx])
    idx_max <- idx[which.min(a[idx] / (a-b)[idx])]
    
    a <- a + lambda_max * (b-a)
    b[idx_max] <- 0
    lambda <- lambda - lambda_max
  }
  if (dls) cat(cnt, "Max(Lambda) steps\n")
  a + lambda * (b-a)
}

search_lambda <- function(prx, Ainv) {
  
  idx <- which.min(Ainv %*% prx)
  
  weights <- Ainv[idx, ]
  
  a <- 0
  b <- 1
  while(b-a > 0.01) {
    
    mid <- (b+a)/2
    cand <- prx * (1-mid) + (mid) * Ainv %*% prx
    if (check_simplex(cand)) {
      b <- mid
    } else {
      a <- mid
    }
  }
  
  prx * (1-a) + (a) * Ainv %*% prx
}
