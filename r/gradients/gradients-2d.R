
r_idx <- function(rk, i, j, k) {
  
  rk_bit <- to_bits(rk - 1, k)
  vapply(
    seq_len(2^k),
    function(x) {
      
      x_bit <- to_bits(x-1, k)
      
      if (any(x_bit < rk_bit)) return(FALSE)
      
      if (x_bit[i] == 0 | x_bit[j] == 0) return(FALSE)
      
      return(TRUE)
    }, logical(1L)
  )
}

# computes P(r | z)
r_fi_hash <- function(rk, fi, k) {
  
  rk_bit <- to_bits(rk - 1, k)
  vapply(
    seq_len(2^k),
    function(z) {
      
      z_bit <- to_bits(z-1, k)
      
      if (any(z_bit < rk_bit)) return(0)
      
      fi ^ sum(z_bit - rk_bit) * (1 - fi) ^ sum(rk_bit == 1)
      
    }, numeric(1L)
  )
}

compute_grad_ij_old <- function(i, j, r, pz, fi, mode = "cpp") {
  
  k <- as.integer(log(length(pz), 2))
  
  rhash <- 1 + rowSums(r * t(replicate(nrow(r), 2^(seq_len(k) - 1))))
  
  k_grad_ij <- parallel::mclapply(
    rhash, function(rk) {
      
      if (mode == "cpp") {
        
        fi_wgh_hash <- cpp_fi_hash(rk, fi, k)
        idx <- cpp_idx(rk, i, j, k)
      } else {
        
        fi_wgh_hash <- r_fi_hash(rk, fi, k)
        idx <- r_idx(rk, i, j, k) 
      }
      
      sum(fi_wgh_hash[idx] * pz[idx]) / sum(fi_wgh_hash * pz)
    }, mc.cores = n_cores()
  )
  
  k_grad_ij <- unlist(k_grad_ij)
  
  if (i != j) k_grad_ij <- k_grad_ij * 2
  
  mean(k_grad_ij)
}

compute_grad_ij <- function(i, j, r, pz, fi, mc_cores = n_cores()) {
  
  k <- as.integer(log(length(pz), 2))
  
  rhash <- 1 + rowSums(r * t(replicate(nrow(r), 2^(seq_len(k) - 1))))
  
  # mc_cores <- 1 # n_cores()
  folds <- rep(seq_len(mc_cores), length(rhash) / 6 + 1)[seq_along(rhash)]
  
  if (mc_cores == 1) {
    
    k_grad_ij <- cpp_grad_ij(rhash, i, j, k, fi, pz)
  } else {
    
    k_grad_ij <- parallel::mclapply(
      seq_len(mc_cores), function(fold) {
        
        rhash_fold <- rhash[folds == fold]
        cpp_grad_ij(rhash_fold, i, j, k, fi, pz)
      }, mc.cores = mc_cores
    )
    
    k_grad_ij <- do.call(c, k_grad_ij)
  }
  
  if (i != j) k_grad_ij <- k_grad_ij * 2
  
  mean(k_grad_ij)
}

likelihood_2d <- function(Sigma, r, fi, mode = "cpp") {
  
  pz <- scaling_2d(Sigma)
  pz <- pz / sum(pz)
  
  k <- as.integer(log(length(pz), 2))
  
  rhash <- 1 + rowSums(r * t(replicate(nrow(r), 2^(seq_len(k) - 1))))
  
  ll_k <- vapply(
    rhash, function(rk) {
      
      if (mode == "cpp") {
        
        fi_wgh_hash <- cpp_fi_hash(rk, fi, k)
      } else {
        
        fi_wgh_hash <- r_fi_hash(rk, fi, k)
      }
      
      log(sum(fi_wgh_hash * pz))
      
    }, numeric(1L)
  )
  
  mean(ll_k)
  
}

compute_grad_logc_ij <- function(i, j, pz) {
  
  dim <- as.integer(log(length(pz), 2))
  idx <- cpp_idx(1, i, j, dim)
  (1 + (i != j)) * 1 / sum(pz) * sum(pz[idx])
}

gradient_2d <- function(Sigma, r, fi) {
  
  pz <- scaling_2d(Sigma)
  pz <- pz / sum(pz)
  grad <- array(0, dim = dim(Sigma))
  k <- ncol(Sigma)
  
  for (i in seq_len(k)) {
    for (j in seq.int(i, k, 1)) {
      
      grad[i, j] <- compute_grad_ij(i, j, r, pz, fi) - 
        compute_grad_logc_ij(i, j, pz)
      
      grad[j, i] <- grad[i, j]
    }
  }
  
  grad
}

descent_2d <- function(Sigma0, r, fi, gamma = 1, n_iter = 50, 
                       true_Sigma = NULL, verbose = FALSE) {
  
  Sigma_dist <- loglc <- rep(0, n_iter)
  Sigma <- list()
  Sigma[[1]] <- Sigma0
  
  loglc[1] <- likelihood_2d(Sigma[[1]], r, fi)
  if (!is.null(true_Sigma)) Sigma_dist[1] <- sum((Sigma[[1]] - true_Sigma)^2)
  
  for (i in seq.int(2, n_iter)) {
    
    step_targ <- gamma / (1 + i)
    curr_grad <- gradient_2d(Sigma[[i-1]], r, fi)
    
    if (!is.null(true_Sigma)) perf_grad <- true_Sigma - Sigma[[i-1]]
    
    
    gamma_i <- step_targ / vec_norm(curr_grad)
    
    Sigma[[i]] <- Sigma[[i-1]] + gamma_i * curr_grad
    
    loglc[i] <- likelihood_2d(Sigma[[i]], r, fi)
    
    if (!is.null(true_Sigma)) Sigma_dist[i] <- sum((Sigma[[i]] - true_Sigma)^2)
    
    if (verbose) cat("Log(L_c) =", loglc[i], "\n")
  }
  
  list(Sigma = Sigma, Sigma_dist = Sigma_dist, loglc = loglc)
}
