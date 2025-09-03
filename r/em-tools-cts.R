
gradient_2d_fi0_cts <- function(cor_mat, Amut, Ez, idx, Sigma0, Sig0Lam, p, d) {
  
  cor_mat[seq_len(p), seq_len(p)] <- 2 * cor_mat[seq_len(p), seq_len(p)]
  cor_mat[-seq_len(p), -seq_len(p)] <- -1 * cor_mat[-seq_len(p), -seq_len(p)]
  diag(cor_mat) <- diag(cor_mat) / 2
  cor_mat[idx] - cpp_logc_ij_cts(Amut, Ez, which(idx) - 1, Sigma0, Sig0Lam, p, d)
} 

likelihood_2d_fi0_cts <- function(theta, cor_mat, d_index) {
  
  Umat <- theta_rearrange(theta, d_index)
  Sigma0 <- solve(theta[d_index, d_index, drop=FALSE])
  
  Ez <- theta_Ez(theta, d_index)
  detSigma0 <- det(Sigma0)
  if (detSigma0 < 0) return(NaN)
  sum(cor_mat * Umat) - log(sum(Ez)) - 1 / 2 * log(det(Sigma0))
}

theta_rearrange <- function(theta, d_index) {
  
  A11 <- theta[-d_index, -d_index]
  A21 <- theta[d_index, -d_index, drop = FALSE]

  A12 <- matrix(0, nrow = nrow(theta) - length(d_index), ncol = length(d_index))
  A22 <- -0.5 * theta[d_index, d_index]
  
  # Combine:
  left <- rbind(A11, A21)
  right <- rbind(A12, A22)
  cbind(left, right)
}

theta_Ez <- function(theta, d_index) {
  
  Sigma0 <- solve(theta[d_index, d_index, drop = FALSE])
  Lambda <- theta[d_index, -d_index, drop = FALSE]
  Sigma <- theta[-d_index, -d_index]
  
  cpp_scaling_2d(Sigma + 1 / 2 * t(Lambda) %*% Sigma0 %*% Lambda)
}

cormat_to_Sigma_cts <- function(cor_mat, k, d, n_iter = 100, multistep = FALSE) {
  
  # initiate the list for each step
  theta_t <- list()
  # theta_t[[1]] <- array(0, dim = c(k+d, k+d))
  # diag(theta_t[[1]])[k+seq_len(d)] <- 1
  # theta_t[[1]][seq_len(k), seq_len(k)] <- cor_mat[seq_len(k), seq_len(k)]
  theta_t[[1]] <- cor_mat
  
  # initiate Armijo-Goldstein parameters
  beta <- 1 / 2
  alpha0 <- 1
  alpha <- alpha0
  
  d_index <- seq.int(k+1, k+d, 1)
  k_index <- seq.int(1, k, 1)
  
  # get indices for the lower-triangular part of Sigma
  idx <- !upper.tri(theta_t[[1]])
  
  for (i in seq.int(2, 1 + n_iter)) {
    
    Sigma <- theta_t[[i-1]][k_index, k_index]
    Omega <- theta_t[[i-1]][d_index, d_index]
    Sigma0 <- solve(Omega)
    Lambda <- theta_t[[i-1]][d_index, k_index, drop = FALSE]
    Sig0Lam <- Sigma0 %*% Lambda
    
    Ez <- cpp_scaling_2d(Sigma + 1 / 2 * t(Lambda) %*% Sigma0 %*% Lambda)
    Amut <- sum(Ez)
    
    # compute gradient
    curr_grad <- gradient_2d_fi0_cts(cor_mat, Amut, Ez, idx, Sigma0, Sig0Lam, k, d)
    
    # make Sigma into a vector
    theta_comp <- theta_t[[i-1]][idx]
    
    # compute the Hessian matrix
    Hess <- cpp_hess_cts_fast(Amut, Ez, which(idx) - 1, Sigma0, Sig0Lam, k, d)
    
    update <- tryCatch(
      solve(Hess, curr_grad),
      error = function(e) e
    )
    
    if (inherits(update, "error")) {
      
      update <- curr_grad
      message("Hessian Singular\n")
    }
    
    curr_obj <- likelihood_2d_fi0_cts(theta_t[[i-1]], cor_mat, d_index)
    theta_cand <- Sigma_pf(theta_t[[i-1]][idx], as.vector(alpha * update), k+d)
    new_obj <- likelihood_2d_fi0_cts(theta_cand, cor_mat, d_index)
    if (is.nan(new_obj)) new_obj <- -Inf
    improve <- new_obj - curr_obj
    
    nabla2 <- vec_norm(curr_grad)^2
    
    # if improvement < 0, need to halve until better
    arm_cnt <- 0
    halved <- FALSE
    while (improve < 0) {
      
      arm_cnt <- arm_cnt + 1
      alpha <- alpha * beta
      theta_cand <- Sigma_pf(theta_t[[i-1]][idx], as.vector(alpha * update), k+d)
      new_obj <- likelihood_2d_fi0_cts(theta_cand, cor_mat, d_index)
      if (is.nan(new_obj)) new_obj <- -Inf
      improve <- new_obj - curr_obj
      if (improve > 0) halved <- TRUE
      
      if (arm_cnt > 30) break
    }
    
    h_bet <- TRUE
    h_bet_cnt <- 0
    while (h_bet) {
      
      h_bet <- FALSE
      
      theta_cand_h <- Sigma_pf(theta_t[[i-1]][idx], 
                               as.vector(alpha * beta * update), k+d)
      new_obj_h <- likelihood_2d_fi0_cts(theta_cand_h, cor_mat, d_index)
      if (is.nan(new_obj_h)) new_obj_h <- -Inf
      improve_h <- new_obj_h - curr_obj
      
      if (improve_h > improve) {
        
        halved <- TRUE
        improve <- improve_h
        theta_cand <- theta_cand_h
        alpha <- alpha * beta
        h_bet_cnt <- h_bet_cnt + 1
        
        if (multistep) h_bet <- TRUE
      }
      
      if (h_bet_cnt > 30) break
    }
    
    # if not halved, check doubling
    if (!halved) {
      
      d_bet <- TRUE
      d_bet_cnt <- 0
      # is doubling better?
      while(d_bet) {
        
        d_bet <- FALSE
        theta_cand_d <- Sigma_pf(theta_t[[i-1]][idx], as.vector(alpha / beta * update), k+d)
        new_obj_d <- likelihood_2d_fi0_cts(theta_cand_d, cor_mat, d_index)
        if (is.nan(new_obj_d)) new_obj_d <- -Inf
        improve_d <- new_obj_d - curr_obj
        
        if (improve_d > improve) {
          
          improve <- improve_d
          theta_cand <- theta_cand_d
          alpha <- alpha / beta
          d_bet_cnt <- d_bet_cnt + 1
          
          if (multistep) d_bet <- TRUE
        }
        
        if (d_bet_cnt > 30) break
      }
    }
    
    new_obj <- likelihood_2d_fi0_cts(theta_cand, cor_mat, d_index)
    improve <- new_obj - curr_obj
    theta_t[[i]] <- theta_cand
    if (improve < 10^(-10)) break
  }
  
  theta_t
}

fit_lmb_rw <- function(X, Y, R, W) {
  
  xfit <- glm(X ~ ., data = data.frame(X = X, R = R, W = W), family = "binomial")
  xz_coeff <- grepl("^R", names(xfit$coefficients))
  xw_coeff <- grepl("^W", names(xfit$coefficients))
  lambda_z <- xfit$coefficients[xz_coeff]
  lambda_w <- xfit$coefficients[xw_coeff]
  lambda_icept <- xfit$coefficients[1]
  
  yfit <- glm(Y ~ ., data = data.frame(Y = Y, X = X, R = R, W = W), 
              family = "binomial")
  yz_coeff <- grepl("^R", names(yfit$coefficients))
  yw_coeff <- grepl("^W", names(yfit$coefficients))
  
  mu_icept <- yfit$coefficients[1]
  mu_z <- yfit$coefficients[yz_coeff]
  mu_w <- yfit$coefficients[yw_coeff]
  beta0 <- yfit$coefficients["X"]
  
  list(lambda_z = lambda_z, lambda_w = lambda_w, lambda_icept = lambda_icept,
       mu_z = mu_z, mu_w = mu_w, mu_icept = mu_icept,
       beta = beta0, beta_se = summary(yfit)$coefficients["X", "Std. Error"])
}