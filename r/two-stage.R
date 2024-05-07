
#' @importFrom utils tail
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

r_descent_2d <- function(Sigma, r, fi, n_iter = 100, multistep = FALSE,
                         verbose = FALSE) {
  
  Sigma_t <- list()
  Sigma_t[[1]] <- Sigma
  beta <- 1 / 2
  alpha0 <- 1
  alpha <- alpha0
  k <- ncol(Sigma_t[[1]])
  
  idx <- !upper.tri(Sigma_t[[1]])
  
  for (i in seq.int(2, 1 + n_iter)) {
    
    # compute gradient
    curr_grad <- gradient_2d_fi0(Sigma_t[[i-1]], r, fi)
    
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
      cat("Hessian Singular\n")
    }
    
    curr_obj <- likelihood_2d_fi0(Sigma_t[[i-1]], r, fi) 
    Sigma_cand <- Sigma_pf(Sigma_t[[i-1]][idx], as.vector(alpha * update), k)
    new_obj <- likelihood_2d_fi0(Sigma_cand, r, fi)
    improve <- new_obj - curr_obj
    
    nabla2 <- vec_norm(curr_grad)^2
    
    # if improvement < 0, need to halve until better
    arm_cnt <- 0
    while (improve < 0) {
      
      arm_cnt <- arm_cnt + 1
      alpha <- alpha * beta
      Sigma_cand <- Sigma_pf(Sigma_t[[i-1]][idx], as.vector(alpha * update), k)
      new_obj <- likelihood_2d_fi0(Sigma_cand, r, fi)
      improve <- new_obj - curr_obj
      if (verbose) cat(arm_cnt, "\n")
    }
    
    halved <- FALSE
    h_bet <- TRUE
    while (h_bet) {
      
      h_bet <- FALSE
      
      Sigma_cand_h <- Sigma_pf(Sigma_t[[i-1]][idx], 
                               as.vector(alpha * beta * update), k)
      new_obj_h <- likelihood_2d_fi0(Sigma_cand_h, r, fi)
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
        new_obj_d <- likelihood_2d_fi0(Sigma_cand_d, r, fi)
        improve_d <- new_obj_d - curr_obj
        
        if (improve_d > improve) {
          
          improve <- improve_d
          Sigma_cand <- Sigma_cand_d
          alpha <- alpha / beta
          
          if (multistep) d_bet <- TRUE
        }
      }
    }
    
    new_obj <- likelihood_2d_fi0(Sigma_cand, r, fi)
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

compute_grad_logc_ij <- function(i, j, pz) {
  
  dim <- as.integer(log(length(pz), 2))
  idx <- cpp_idx(1, i, j, dim)
  (1 + (i != j)) * 1 / sum(pz) * sum(pz[idx])
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
