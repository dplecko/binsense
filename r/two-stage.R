
#' @importFrom utils tail
infer_Sigma_IM <- function(X, Y, R, fi) {
  
  cor_mat <- XYR_to_cormat(X, Y, R, fi)
  tail(cormat_to_Sigma(cor_mat), n = 1)[[1]]
}

XYR_to_cormat <- function(X, Y, R, fi) {
  
  k <- ncol(R)
  mu_hat <- array(0, dim = c(k, k))
  
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      idx <- X == x & Y == y
      Rxy <- R[idx, , drop = FALSE]
      
      fixy <- fi[[1+x]][[1+y]]
      if (length(fixy) == 1 & ncol(R) > 1) fixy <- rep(fixy, ncol(R))
      
      fi_mat <- (1 - fixy) %*% t(1-fixy)
      diag(fi_mat) <- 1 - fixy
      
      mu_hat <- mu_hat + (t(Rxy) %*% Rxy / nrow(Rxy) / fi_mat) * mean(idx)
    }
  }
  
  mu_hat
}

cormat_to_Sigma <- function(cor_mat, n_iter = 100, multistep = FALSE, verbose = FALSE) {
  
  # get dimension k based on the cor_mat
  k <- ncol(cor_mat)
  
  # initiate the list for each step
  Sigma_t <- list()
  Sigma_t[[1]] <- array(0, dim = c(k, k))
  
  # initiate Armijo-Goldstein parameters
  beta <- 1 / 2
  alpha0 <- 1
  alpha <- alpha0
  
  # get indices for the lower-triangular part of Sigma
  idx <- !upper.tri(Sigma_t[[1]])
  
  for (i in seq.int(2, 1 + n_iter)) {
    
    # compute gradient
    curr_grad <- gradient_2d_fi0(Sigma_t[[i-1]], cor_mat)
    
    # make Sigma into a vector
    Sigma_comp <- Sigma_t[[i-1]][idx]
    
    # compute the Hessian matrix
    pz <- cpp_scaling_2d(Sigma_t[[i-1]])
    Hess <- cpp_hessian_2d(pz / sum(pz), which(idx) - 1, ncol(Sigma_t[[i-1]]))
    
    update <- tryCatch(
      solve(Hess, curr_grad[idx]),
      error = function(e) e
    )
    
    if (inherits(update, "error")) {
      
      update <- curr_grad[idx]
      message("Hessian Singular\n")
    }
    
    curr_obj <- likelihood_2d_fi0(Sigma_t[[i-1]], cor_mat) 
    Sigma_cand <- Sigma_pf(Sigma_t[[i-1]][idx], as.vector(alpha * update), k)
    new_obj <- likelihood_2d_fi0(Sigma_cand, cor_mat)
    improve <- new_obj - curr_obj
    
    nabla2 <- vec_norm(curr_grad)^2
    
    # if improvement < 0, need to halve until better
    arm_cnt <- 0
    while (improve < 0) {
      
      arm_cnt <- arm_cnt + 1
      alpha <- alpha * beta
      Sigma_cand <- Sigma_pf(Sigma_t[[i-1]][idx], as.vector(alpha * update), k)
      new_obj <- likelihood_2d_fi0(Sigma_cand, cor_mat)
      improve <- new_obj - curr_obj
      if (verbose) cat(arm_cnt, "\n")
    }
    
    halved <- FALSE
    h_bet <- TRUE
    while (h_bet) {
      
      h_bet <- FALSE
      
      Sigma_cand_h <- Sigma_pf(Sigma_t[[i-1]][idx], 
                               as.vector(alpha * beta * update), k)
      new_obj_h <- likelihood_2d_fi0(Sigma_cand_h, cor_mat)
      improve_h <- new_obj_h - curr_obj
      
      if (improve_h > improve) {
        
        halved <- TRUE
        improve <- improve_h
        Sigma_cand <- Sigma_cand_h
        alpha <- alpha * beta
        
        if (multistep) h_bet <- TRUE
      }
    }
    
    # if not halved, check doubling
    if (!halved) {
      
      d_bet <- TRUE
      # is doubling better?
      while(d_bet) {
        
        d_bet <- FALSE
        Sigma_cand_d <- Sigma_pf(Sigma_t[[i-1]][idx], as.vector(alpha / beta * update), k)
        new_obj_d <- likelihood_2d_fi0(Sigma_cand_d, cor_mat)
        improve_d <- new_obj_d - curr_obj
        
        if (improve_d > improve) {
          
          improve <- improve_d
          Sigma_cand <- Sigma_cand_d
          alpha <- alpha / beta
          
          if (multistep) d_bet <- TRUE
        }
      }
    }
    
    new_obj <- likelihood_2d_fi0(Sigma_cand, cor_mat)
    if (verbose) {
      
      cat("Likelihood at iteration", i, "=", new_obj, 
          "\n")
    }
    improve <- new_obj - curr_obj
    Sigma_t[[i]] <- Sigma_cand
    if (improve < 10^(-10)) break
  }
  
  Sigma_t
}

Sigma_pf <- function(Sigma_comp, update, k) {
  
  Sigma_comp <- Sigma_comp + update
  
  Sigma_pf <- array(0, dim = c(k, k))
  Sigma_pf[!upper.tri(Sigma_pf)] <- Sigma_comp
  
  Sigma_pf <- (Sigma_pf + t(Sigma_pf))
  diag(Sigma_pf) <- diag(Sigma_pf) / 2
  
  Sigma_pf
}

likelihood_2d_fi0 <- function(Sigma, cor_mat) {
  
  pz <- cpp_scaling_2d(Sigma)
  sum(cor_mat * Sigma) - log(sum(pz))
}

compute_grad_logc_ij <- function(i, j, pz) {
  
  dim <- as.integer(log(length(pz), 2))
  idx <- cpp_idx(1, i, j, dim)
  (1 + (i != j)) * 1 / sum(pz) * sum(pz[idx])
}

gradient_2d_fi0 <- function(Sigma, cor_mat) {
  
  pz <- cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  grad <- grad_ij_fi0(cor_mat)
  k <- ncol(Sigma)
  
  for (i in seq_len(k)) {
    for (j in seq.int(i, k, 1)) {
      
      grad[i, j] <- grad[i, j] - compute_grad_logc_ij(i, j, pz)
      
      grad[j, i] <- grad[i, j]
    }
  }
  
  grad
}

grad_ij_fi0 <- function(cor_mat) {
  
  grad <- 2 * cor_mat
  diag(grad) <- diag(grad) / 2
  
  grad
}
