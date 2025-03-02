
z_2d <- function(k, n_samp, seed = 2022, Sigma) {
  
  set.seed(seed)
  k <- ncol(Sigma)
  pz <- cpp_scaling_2d(Sigma)
  
  z <- sample.int(2^k, size = n_samp, prob = pz, replace = TRUE)
  cpp_idx_to_bit(z-1, k)
}

param_mod_gen <- function(n, k, seed, A, Sigma, lam, mu, icept_x, icept_y, beta) {
  
  set.seed(seed)
  
  Z <- z_2d(k, n, seed, Sigma)
  X <- rbinom(n, 1, expit(Z %*% lam + icept_x))
  u <- runif(n)
  Y <- as.integer(u < expit(Z %*% mu + X * beta + icept_y))
  
  R <- Z
  for (x in c(T, F)) {
    
    for (y in c(T, F)) {
      
      idx <- X == x & Y == y
      for (i in seq_len(k)) {
        
        # get the cpp bit to idx, idx to bit
        z_idx <- 1 + cpp_bit_to_idx_mat(Z[idx, ])
        r_idx <- vapply(
          z_idx, 
          function(z_id) {
            
            sample.int(n = 2^k, size = 1, prob = A[[1+x]][[1+y]][, z_id])
          }, numeric(1L)
        )
        R[idx, ] <- cpp_idx_to_bit(r_idx-1, k)
      }
    }
  }
  
  list(X = X, Y = Y, R = R)
}