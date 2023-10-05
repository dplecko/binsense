
infer_Sigma <- function(X, Y, R, fi) {
  
  k <- ncol(R)
  mu_hat <- array(0, dim = c(k, k))
  
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      idx <- X == x & Y == y
      Rxy <- R[idx, , drop = FALSE]
      
      fixy <- fi[[1+x]][[1+y]]
      fi_mat <- array((1-fixy)^2, dim = c(k, k))
      diag(fi_mat) <- 1-fixy
      
      mu_hat <- mu_hat + (t(Rxy) %*% Rxy / nrow(Rxy) / fi_mat) * mean(idx)
    }
  }
  
  attr(mu_hat, "cor_mat") <- TRUE
  tail(r_descent_2d(array(0, dim = c(k, k)), r = mu_hat, fi = fi), n = 1)[[1]]
}