#' * Newton-Raphson with Armijo-Goldstein *

r_descent_2d <- function(Sigma, r, fi, n_iter = 100, multistep = FALSE) {
  
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
      
      update <- curr_grad
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
      cat(arm_cnt, "\n")
    }
    
    #' * further improvement could be to make halving / doubling multi-step *
    # is halving better?
    
    halved <- FALSE
    h_bet <- TRUE
    while (h_bet) {
      
      h_bet <- FALSE
      
      Sigma_cand_h <- Sigma_pf(Sigma_t[[i-1]][idx], as.vector(alpha * beta * update), k)
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
    cat("Likelihood at iteration", i, "=", new_obj, 
        "\n")
    improve <- new_obj - curr_obj
    Sigma_t[[i]] <- Sigma_cand
    if (improve < 0) browser()
    if (improve < 10^(-10)) break
  }
  
  Sigma_t
}
