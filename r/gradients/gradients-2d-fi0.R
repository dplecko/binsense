
likelihood_2d_fi0 <- function(Sigma, z, fi = 0) {
  
  if (!is.null(attr(z, "cor_mat"))) {

    cor_mat <- z
  } else {
    
    nsamp <- nrow(z)
    cor_mat <- 1 / nsamp * t(z) %*% z
    fi_mat <- array((1-fi)^2, dim = dim(cor_mat))
    diag(fi_mat) <- 1-fi
    cor_mat <- cor_mat / fi_mat
  }
  
  pz <- cpp_scaling_2d(Sigma)
  sum(cor_mat * Sigma) - log(sum(pz))
}

gradient_2d_fi0 <- function(Sigma, z, fi = 0) {
  
  pz <- cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  grad <- grad_ij_fi0(z, fi)
  k <- ncol(Sigma)
  
  for (i in seq_len(k)) {
    for (j in seq.int(i, k, 1)) {
      
      grad[i, j] <- grad[i, j] - compute_grad_logc_ij(i, j, pz)
      
      grad[j, i] <- grad[i, j]
    }
  }
  
  grad
}

grad_ij_fi0 <- function(z, fi = 0) {
  
  if (!is.null(attr(z, "cor_mat"))) {

    cor_mat <- z
  } else {
    
    nsamp <- nrow(z)
    cor_mat <- 1 / nsamp * t(z) %*% z
    fi_mat <- array((1-fi)^2, dim = dim(cor_mat))
    diag(fi_mat) <- 1-fi
    cor_mat <- cor_mat / fi_mat
  }
  
  grad <- 2 * cor_mat
  diag(grad) <- diag(grad) / 2
  
  grad
}