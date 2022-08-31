
A_xy <- function(fi_xy, k, order = Inf) {
  
  Axy <- array(0, dim = c(2^k, 2^k))
  
  for (i in seq_len(2^k)) {
    
    for (j in seq_len(2^k)) {
      
      i_bin <- as.integer(intToBits(i-1))
      j_bin <- as.integer(intToBits(j-1))
      
      if (any(i_bin > j_bin)) next
      if (all(i_bin == j_bin)) next
      if (sum(j_bin - i_bin) > order) next
      Axy[i, j] <- fi_xy^sum(j_bin - i_bin) * (1 - fi_xy)^sum(i_bin)
      
    }
  }
  
  Axy <- Axy + diag(1 - colSums(Axy))
  
  assert_that(!is.element(0, eigen(Axy)$values),
              msg = "Singular transition matrix")
  
  Axy
}

invert_A_xy <- function(A, pz) {
  
  B <- array(0, dim = dim(A))
  k <- as.integer(log(nrow(B), 2))
  
  for (i in seq_len(2^k)) {
    
    for (j in seq_len(2^k)) {
      
      i_bin <- as.integer(intToBits(i-1))
      j_bin <- as.integer(intToBits(j-1))

      B[i, j] <- A[j, i] * pz[i] / sum(A[j, ] * pz)
      
    }
  }
  
  B
}
