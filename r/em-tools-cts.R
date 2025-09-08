
gradient_2d_fi0_cts <- function(cor_mat, m_w, Amut, Ez, idx, 
                                Sigma0, Sig0Lam, Sig0b, p, d) {
  
  cor_mat[seq_len(p), seq_len(p)] <- 2 * cor_mat[seq_len(p), seq_len(p)]
  cor_mat[-seq_len(p), -seq_len(p)] <- -1 * cor_mat[-seq_len(p), -seq_len(p)]
  diag(cor_mat) <- diag(cor_mat) / 2
  c(cor_mat[idx], m_w) - 
    cpp_logc_ij_cts(Amut, Ez, which(idx) - 1, Sigma0, Sig0Lam, Sig0b, p, d)
} 

likelihood_2d_fi0_cts <- function(theta, cor_mat, m_w, d_index) {

  Umat <- SLO_for_like(theta[["SLO"]], d_index)
  Sigma0 <- solve(theta[["SLO"]][d_index, d_index, drop=FALSE])
  
  Ez <- theta_Ez(theta, d_index)
  detSigma0 <- det(Sigma0)
  if (detSigma0 < 0) return(NaN)
  sum(cor_mat * Umat) + sum(theta[["b"]] * m_w) - 
    log(sum(Ez)) - 1 / 2 * log(det(Sigma0)) # normalizing constants
}

SLO_for_like <- function(SLO, d_index) {
  
  A11 <- SLO[-d_index, -d_index]
  A21 <- SLO[d_index, -d_index, drop = FALSE]

  A12 <- matrix(0, nrow = nrow(SLO) - length(d_index), ncol = length(d_index))
  A22 <- -0.5 * SLO[d_index, d_index]
  
  # combine
  left <- rbind(A11, A21)
  right <- rbind(A12, A22)
  cbind(left, right)
}

theta_Ez <- function(theta, d_index) {
  
  SLO <- theta[["SLO"]]
  b <- theta[["b"]]
  
  Sigma0 <- solve(SLO[d_index, d_index, drop = FALSE])
  Lambda <- SLO[d_index, -d_index, drop = FALSE]
  Sigma <- SLO[-d_index, -d_index]
  
  Ez <- cpp_scaling_2d(Sigma + 1 / 2 * t(Lambda) %*% Sigma0 %*% Lambda +
                       diag(as.vector(b %*% Sigma0 %*% Lambda)))
  Ez * exp(1 / 2 * as.vector(b %*% Sigma0 %*% b))
}

cormat_to_Sigma_cts <- function(cor_mat, m_w, k, d, n_iter = 100, multistep = FALSE) {
  
  # initiate the list for each step
  theta_t <- list()
  theta_t[[1]] <- list(SLO = cor_mat, b = m_w)
  
  # initiate Armijo-Goldstein parameters
  beta <- 1 / 2
  alpha0 <- 1
  alpha <- alpha0
  
  d_index <- seq.int(k+1, k+d, 1)
  k_index <- seq.int(1, k, 1)
  
  # get indices for the lower-triangular part of SLO matrix
  idx <- !upper.tri(theta_t[[1]][["SLO"]])
  
  for (i in seq.int(2, 1 + n_iter)) {
    
    # extract Sigma, Omega, Lambda from SLO
    Sigma <- theta_t[[i-1]][["SLO"]][k_index, k_index]
    Omega <- theta_t[[i-1]][["SLO"]][d_index, d_index]
    Lambda <- theta_t[[i-1]][["SLO"]][d_index, k_index, drop = FALSE]
    b <- theta_t[[i-1]][["b"]]
    
    Sigma0 <- solve(Omega)
    Sig0Lam <- Sigma0 %*% Lambda
    Sig0b <- as.vector(Sigma0 %*% b)
    
    Ez <- theta_Ez(theta_t[[i-1]], d_index)
    Amut <- sum(Ez)
    
    # compute gradient
    curr_grad <- gradient_2d_fi0_cts(cor_mat, m_w, Amut, Ez, idx, 
                                     Sigma0, Sig0Lam, Sig0b, k, d)
    
    # compute the Hessian matrix
    Hess <- cpp_hess_cts_fast(Amut, Ez, which(idx) - 1, 
                              Sigma0, Sig0Lam, Sig0b, k, d)

    update <- tryCatch(
      solve(Hess, curr_grad),
      error = function(e) e
    )
    
    if (inherits(update, "error")) {
      
      update <- curr_grad
      message("Hessian Singular\n")
    }
    
    curr_obj <- likelihood_2d_fi0_cts(theta_t[[i-1]], cor_mat, m_w, d_index)
    
    theta_pf <- function(theta, update, k, d) {
      
      SLO <- theta[["SLO"]]
      b <- theta[["b"]]
      
      SLO_del <- update[seq_len(length(update) - d)]
      
      SLO_up <- array(0, dim = c(k+d, k+d))
      SLO_up[!upper.tri(SLO_up)] <- SLO_del
      
      SLO_up <- SLO_up + t(SLO_up)
      diag(SLO_up) <- diag(SLO_up) / 2
      
      SLO_pf <- SLO + SLO_up
      b_pf <- b + tail(update, n = d)
      
      list(SLO = SLO_pf, b = b_pf)
    }
    
    theta_cand <- theta_pf(theta_t[[i-1]], as.vector(alpha * update), k, d)
    new_obj <- likelihood_2d_fi0_cts(theta_cand, cor_mat, m_w, d_index)
    if (is.nan(new_obj)) new_obj <- -Inf
    improve <- new_obj - curr_obj
    
    nabla2 <- vec_norm(curr_grad)^2
    
    # if improvement < 0, need to halve until better
    arm_cnt <- 0
    halved <- FALSE
    while (improve < 0) {
      
      arm_cnt <- arm_cnt + 1
      alpha <- alpha * beta
      theta_cand <- theta_pf(theta_t[[i-1]], as.vector(alpha * update), k, d)
      new_obj <- likelihood_2d_fi0_cts(theta_cand, cor_mat, m_w, d_index)
      if (is.nan(new_obj)) new_obj <- -Inf
      improve <- new_obj - curr_obj
      if (improve > 0) halved <- TRUE
      
      if (arm_cnt > 30) break
    }
    
    h_bet <- TRUE
    h_bet_cnt <- 0
    while (h_bet) {
      
      h_bet <- FALSE
      
      theta_cand_h <- theta_pf(theta_t[[i-1]], as.vector(alpha * beta * update), 
                               k, d)
      new_obj_h <- likelihood_2d_fi0_cts(theta_cand_h, cor_mat, m_w, d_index)
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
        theta_cand_d <- theta_pf(theta_t[[i-1]], as.vector(alpha / beta * update), 
                                 k, d)
        new_obj_d <- likelihood_2d_fi0_cts(theta_cand_d, cor_mat, m_w, d_index)
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
    
    new_obj <- likelihood_2d_fi0_cts(theta_cand, cor_mat, m_w, d_index)
    improve <- new_obj - curr_obj
    theta_t[[i]] <- theta_cand
    if (improve < 10^(-10)) break
  }
  
  theta_t
}

lmb_infer <- function(lmb, data = NULL, type = c("logitx_z", "logity_z", 
                                                 "logitx_w", "logity_w", 
                                                 "logitx_zw", "logity_zw")) {
  
  type <- match.arg(type, c("logitx_z", "logity_z", 
                            "logitx_w", "logity_w", 
                            "logitx_zw", "logity_zw"))
  
  # data can only be NULL for x_z, y_z
  if (is.null(data)) {
    
    if (type == "logitx_z") return(cpp_p_z(lmb$lambda_z) + lmb$lambda_icept)
    if (type == "logity_z") return(cpp_p_z(lmb$mu_z) + lmb$mu_icept)
  } else {
    
    w <- data[["w"]]
    x <- data[["x"]]
    z <- data[["z"]]
    
    if (type == "logitx_w") {
      
      w_lin <- sum(lmb$lambda_w * w)
      if (!is.null(lmb$tau_xw)) 
        w_lin <- w_lin + sum(lmb$lambda_tauw * lmb$tau_xw(w))
      return(w_lin)
    } else if (type == "logity_w") {
      
      w_lin <- sum(lmb$mu_w * w) + lmb$beta * x
      if (!is.null(lmb$tau_yw)) 
        w_lin <- w_lin + sum(lmb$mu_tauw * lmb$tau_yw(w))
      return(w_lin)
    } else if (type == "logity_zw") { # this case is vectorized
      
      z_lin <- z %*% lmb$mu_z + lmb$mu_icept
      w_lin <- w %*% lmb$mu_w
      if (!is.null(lmb$tau_yw))
        w_lin <- w_lin + lmb$tau_yw(w) %*% lmb$mu_tauw
      return(z_lin + w_lin)
    }
  }
}

fit_lmb_rw <- function(X, Y, R, W, tau_xw = NULL, tau_yw = NULL) {
  
  colnames(R) <- colnames(W) <- NULL
  xdat <- data.frame(X = X, R = R, W = W)
  if (!is.null(tau_xw)) xdat <- cbind(xdat, tauW = tau_xw(W))
  xfit <- glm(X ~ ., data = xdat, family = "binomial")
  xz_coeff <- grepl("^R", names(xfit$coefficients))
  xw_coeff <- grepl("^W", names(xfit$coefficients))
  lambda_z <- xfit$coefficients[xz_coeff]
  lambda_w <- xfit$coefficients[xw_coeff]
  lambda_icept <- xfit$coefficients[1]
  
  if (!is.null(tau_xw)) {
    
    xtauw_coeff <- grepl("^tauW", names(xfit$coefficients))
    lambda_tauw <- xfit$coefficients[xtauw_coeff]
  }
  
  ydat <- data.frame(Y = Y, X = X, R = R, W = W)
  if (!is.null(tau_yw)) ydat <- cbind(ydat, tauW = tau_yw(W))
  yfit <- glm(Y ~ ., data = ydat, family = "binomial")
  yz_coeff <- grepl("^R", names(yfit$coefficients))
  yw_coeff <- grepl("^W", names(yfit$coefficients))
  
  mu_icept <- yfit$coefficients[1]
  mu_z <- yfit$coefficients[yz_coeff]
  mu_w <- yfit$coefficients[yw_coeff]
  beta0 <- yfit$coefficients["X"]
  
  if (!is.null(tau_yw)) {
    
    ytauw_coeff <- grepl("^tauW", names(yfit$coefficients))
    mu_tauw <- yfit$coefficients[ytauw_coeff]
  }
  
  lmb <- list(lambda_z = lambda_z, lambda_w = lambda_w, lambda_icept = lambda_icept,
              mu_z = mu_z, mu_w = mu_w, mu_icept = mu_icept,
              beta = beta0, beta_se = summary(yfit)$coefficients["X", "Std. Error"])
  
  if (!is.null(tau_xw)) {
    lmb$lambda_tauw <- lambda_tauw
    lmb$tau_xw <- tau_xw
  }
  
  if (!is.null(tau_yw)) {
    lmb$mu_tauw <- mu_tauw
    lmb$tau_yw <- tau_yw
  }
  
  lmb
}