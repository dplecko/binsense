
r_mu <- function(Sigma) {
  
  k <- ncol(Sigma)
  pz <- cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  cpp_mu(pz, k)
}

mom_obj <- function(Sigma, mu_hat) {
  
  k <- ncol(Sigma)
  pz <- cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  mu_sig <- cpp_mu(pz, k)
  
  sum((mu_sig - mu_hat)^2)
}

mom_gradient <- function(Sigma, mu_hat, alpha0, nsamp,
                         reset_alpha = F, allow_speedup = !reset_alpha,
                         verbose = F) {
  
  beta <- 1/2
  
  k <- ncol(Sigma)
  idx <- which(!upper.tri(Sigma))
  iter <- z_cnt <- nz_cnt <- 0
  curr_obj <- Inf
  
  # choose the stopping cleverly 
  c_rad <- mean(mu_hat * (1 - mu_hat) / nsamp)
  cat("Convergence radius", 10^4 * c_rad, "times smaller than before\n")
  stopifnot(c_rad > 0)
  
  arm_hist <- alpha_hist <- NULL
  alpha <- alpha0
  while (curr_obj > 1 / 500 * c_rad * k^2) { # you can choose the stopping cleverly
    
    iter <- iter + 1
    if (reset_alpha) alpha <- alpha0
    
    
    pz <- cpp_scaling_2d(Sigma)
    pz <- pz / sum(pz)
    mu <- cpp_mu(pz, k)
    curr_obj <- mom_obj(Sigma, mu_hat)
    delta_mu <- mu - mu_hat
    hess <- cpp_hessian_2d(pz / sum(pz), idx - 1, k)
    grad <- hess %*% delta_mu[idx]
    
    # cat("\n", round(as.vector(grad)[1:5], 3), "\n")
    
    Sigma_cand <- Sigma_pf(Sigma[idx], as.vector(alpha * grad), k)
    new_obj <- mom_obj(Sigma_cand, mu_hat)
    
    improve <- curr_obj - new_obj # expecting this to be positive (minimizing!)
    nabla2 <- vec_norm(grad)^2
    arm_cnt <- 0
    while (improve < abs(nabla2 * alpha / 2)) {
      
      arm_cnt <- arm_cnt + 1
      alpha <- alpha * beta
      Sigma_cand <- Sigma_pf(Sigma[idx], as.vector(alpha * grad), k)
      new_obj <- mom_obj(Sigma_cand, mu_hat)
      improve <- curr_obj - new_obj
    }
    
    # if (iter == 4) browser()
    
    if (arm_cnt == 0 & allow_speedup) {

      while (improve >= abs(nabla2 * alpha / 2)) {

        arm_cnt <- arm_cnt - 1
        alpha <- alpha / beta
        Sigma_cand <- Sigma_pf(Sigma[idx], as.vector(alpha * grad), k)
        new_obj <- mom_obj(Sigma_cand, mu_hat)
        improve <- curr_obj - new_obj
        
      }

      alpha <- alpha * beta # return to the optimal value
      Sigma_cand <- Sigma_pf(Sigma[idx], as.vector(alpha * grad), k)
      arm_cnt <- arm_cnt + 1
    }
    
    if (verbose) {
      
      if (arm_cnt < -1) cat("\nSpeed-Up triggered", "\n")
      if (arm_cnt > 0) cat("\n Slow-Down triggered", "\n")
      if (iter %% 1000 == 0) cat(log(-1/alpha, 2))
    }
    if (arm_cnt != 0 & (iter %% 10000) == 0) {

      cat("Armijo count needed", arm_cnt, "in iter", iter,
          "with objective", new_obj,"\n")
    }
    
    arm_hist <- c(arm_hist, arm_cnt)
    alpha_hist <- c(alpha_hist, alpha)
    Sigma <- Sigma_cand
  }
  
  # return(arm_hist)
  # cat(arm_hist)
  # cat(alpha_hist)
  return(Sigma)
}
