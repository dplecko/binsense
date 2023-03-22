
likelihood_2d_fi0 <- function(Sigma, z) {
  
  pz <- cpp_scaling_2d(Sigma)
  1/ nrow(z) * sum(t(z) %*% z * Sigma) - log(sum(pz))
}

gradient_2d_fi0 <- function(Sigma, z) {
  
  pz <- cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  grad <- grad_ij_fi0(z)
  k <- ncol(Sigma)
  
  for (i in seq_len(k)) {
    for (j in seq.int(i, k, 1)) {
      
      grad[i, j] <- grad[i, j] - compute_grad_logc_ij(i, j, pz)
      
      grad[j, i] <- grad[i, j]
    }
  }
  
  grad
}

grad_ij_fi0 <- function(z) {
  
  nsamp <- nrow(z)
  grad <- 2 / nsamp * t(z) %*% z
  diag(grad) <- diag(grad) / 2
  
  grad
}