
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
  
  assertthat::assert_that(!is.element(0, eigen(Axy)$values),
                          msg = "Singular transition matrix")
  
  Axy
}

invert_A_xy <- function(A, pz, pr) {
  
  B <- array(0, dim = dim(A))
  k <- as.integer(log(nrow(B), 2))
  
  for (i in seq_len(2^k)) {
    
    for (j in seq_len(2^k)) {
      
      i_bin <- as.integer(intToBits(i-1))
      j_bin <- as.integer(intToBits(j-1))

      B[i, j] <- A[j, i] * pz[i] / pr[j]
      
    }
  }
  
  B
}

infer_pz <- function(z, pr, phi, k, weights = FALSE) {
  
  zbit <- as.integer(intToBits(z-1)[seq.int(1, k, 1)])
  idx <- which(zbit == 0)
  
  if (length(idx) == 0) {
    r_geq_z <- list(integer(0L))
  } else if (length(idx) == 1) {
    r_geq_z <- list(integer(0L), idx)
  } else {
    
    r_geq_z <- unlist(lapply(0:length(idx),
                             combn, 
                             x = idx,
                             simplify = FALSE), 
                      recursive = FALSE)
  }
  # browser()
  if (length(phi) == 1) {
    
    r_wgh <- unlist(
      lapply(r_geq_z, function(r) (-1)^length(r) * phi^length(r) * 
               (1-phi)^(-length(r) - sum(zbit > 0)))
    )
  } else if (length(phi) == k) {
    
    r_wgh <- unlist(
      lapply(r_geq_z, function(r) (-1)^length(r) * 
               prod(phi[r]) / 
               prod(1 - phi[c(r, which(zbit > 0))]))
    )
  }
  
  r_idx <- unlist(lapply(r_geq_z, function(r) sum(2^((1:k) - 1)[r]) + z))
  
  if (weights) {
    wgh <- rep(0, 2^k)
    wgh[r_idx] <- r_wgh
    return(wgh)
  }
  
  sum(r_wgh * pr[r_idx])
}

# # explore order of the sup norm
# 
# k_seq <- seq.int(2, 10)
# fi <- 0.2
# res <- lapply(
#   k_seq,
#   function(k) {
#     B <- solve(A_xy(fi, k))
#     c(max(eigen(B)$values), max(B), (1-fi)^(-k))
#   }
# )
# 
# df <- data.frame(do.call(rbind, res))
# df <- cbind(df, k = k_seq)
# 
# ggplot(df, aes(x = k, y = X1)) +
#   geom_line() + geom_point() + theme_bw()

