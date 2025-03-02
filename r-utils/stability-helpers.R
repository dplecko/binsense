
px_z <- function(th) {
  
  lambda <- th$lambda 
  l_icept <- th$l_icept
  
  k <- length(lambda)
  all_z <- lapply(
    seq_len(2^k) - 1,
    function(i) binsensate:::cpp_idx_to_bit(i, k)
  )
  all_z <- do.call(rbind, all_z)
  
  expit(as.vector(all_z %*% lambda) + l_icept)
}

py_xz <- function(th, x) {
  
  mu <- th$mu
  beta <- th$beta
  m_icept <- th$m_icept
  
  k <- length(mu)
  all_z <- lapply(
    seq_len(2^k) - 1,
    function(i) binsensate:::cpp_idx_to_bit(i, k)
  )
  all_z <- do.call(rbind, all_z)
  expit(as.vector(all_z %*% mu) + m_icept + beta * x)
}

p_v <- function(th, fi, x, y) {
  
  Sigma <- th$Sigma
  k <- ncol(Sigma)
  
  pz <- binsensate:::cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  
  pxz <- px_z(th)
  if (x == 0) pxz <- 1 - pxz
  
  pyxz <- py_xz(th, x)
  if (y == 0) pyxz <- 1 - pyxz
  
  phid <- pz * pxz * pyxz
  
  prxyz <- binsensate:::cpp_A_xy(fi, k)
  
  t(t(prxyz) * phid) # r varies across rows, z across columns (!)
}

nabla_Q <- function(th_p, th, th_s, fi, full = FALSE) {
  
  k <- ncol(th_s$Sigma)
  
  all_z <- lapply(
    seq_len(2^k) - 1,
    function(i) binsensate:::cpp_idx_to_bit(i, k)
  )
  all_z <- do.call(rbind, all_z)
  
  # get P_{\theta^*}(r, x, y)
  grad_tot <- rep(0, 2 * k + 3)
  grad_termsx <- grad_termsy <- list(list(NULL, NULL), list(NULL, NULL))
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      # get P_{\theta^*}(r, x, y)
      pv_s <- p_v(th_s, fi, x, y)
      p_rxy_s <- rowSums(pv_s)
      
      # get P_{\theta}(z | r, x, y)
      pv <- p_v(th, fi, x, y)
      pz_rxy <- pv / rowSums(pv)
      
      # get x - P_{\theta'}(x_1 | z)
      term_x <- x - px_z(th_p)
      # term_x[TRUE] <- 1
      
      pz_rxy_tx <- t(t(pz_rxy) * term_x)
      grad_x <- colSums((pz_rxy_tx * p_rxy_s) %*% cbind(1, all_z))
      
      grad_termsx[[1+x]][[1+y]] <- pz_rxy_tx %*% cbind(1, all_z)
      
      
      # cat("Check that", pz_rxy_tx * p_rxy, "=", colSums(p_v(th, fi, x, y)) * term_x, "\n")
      
      # get y - P_{\theta'}(y_1 | z, x)
      term_y <- y - py_xz(th_p, x)
      # term_y[TRUE] <- 1
      
      pz_rxy_ty <- t(t(pz_rxy) * term_y)
      grad_y <- colSums((pz_rxy_ty * p_rxy_s) %*% cbind(1, all_z, x))
      
      grad_termsy[[1+x]][[1+y]] <- pz_rxy_ty %*% cbind(1, all_z, x)
      
      grad_tot <- grad_tot + c(grad_x, grad_y)
    }
  }
  
  if (!full) return(grad_tot)
  list(grad_tot = grad_tot, grad_termsx = grad_termsx, 
       grad_termsy = grad_termsy)
}

# find approximate maximizer at some theta given the true theta_star
M_theta <- function(th, th_s, fi, n_mc = 10^5) {
  
  k <- ncol(th_s$Sigma)
  all_z <- lapply(
    seq_len(2^k) - 1,
    function(i) binsensate:::cpp_idx_to_bit(i, k)
  )
  all_z <- do.call(rbind, all_z)
  
  # compute the maximizer based on Monte Carlo sampling
  fit_lmb <- function(X, Y, R) {
    
    xfit <- glm(X ~ ., data = data.frame(X = X, R = R), family = "binomial")
    lambda0 <- xfit$coefficients[-1]
    lambda_icept <- xfit$coefficients[1]
    
    yfit <- glm(Y ~ ., data = data.frame(Y = Y, X = X, R = R), 
                family = "binomial")
    r_coeff <- grepl("^R", names(yfit$coefficients))
    
    mu_icept <- yfit$coefficients[1]
    mu0 <- yfit$coefficients[r_coeff]
    beta0 <- yfit$coefficients["X"]
    
    list(lambda = lambda0, mu = mu0, beta = beta0,
         l_icept = lambda_icept, m_icept = mu_icept)
  }
  
  X <- Y <- Z <- NULL
  # get P(x, y)
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      p_v <- p_v(th_s, fi, x, y)
      
      n_mc_xy <- round(sum(p_v) * n_mc)
      
      pz_rxy <- p_v(th, fi, x, y)
      pz_rxy <- pz_rxy / rowSums(pz_rxy)
      
      p_MC <- colSums(rowSums(p_v) * pz_rxy)
      samp <- sample.int(2^k, n_mc_xy, prob = p_MC, replace = TRUE)
      
      X <- c(X, rep(x, n_mc_xy))
      Y <- c(Y, rep(y, n_mc_xy))
      Z <- rbind(Z, all_z[samp, ])
    }
  }
  
  fit <- fit_lmb(X, Y, Z)
  fit$Sigma <- th_s$Sigma
  fit
}

# theta to vector
em_params <- function(th) {
  
  c(th$lambda, th$mu, th$beta, th$l_icept, th$m_icept)
}

# generate theta randomly (multiple options)
gen_params <- function(k, rand = FALSE, eps = 0) {
  
  if (!rand) {
    
    Sigma <- matrix(0, ncol = k, nrow = k)
    lambda <- mu <- rep(eps, k)
    m_icept <- l_icept <- eps
    beta <- eps
  } else {
    
    Sigma <- gen_Sigma(k, real = FALSE)
    lambda <- runif(k, -1, 1)
    l_icept <- runif(1, -1, 1)
    mu <- runif(k, -1, 1)
    m_icept <- runif(1, -1, 1)
    beta <- runif(1, -1, 1)
  }
  
  list(Sigma = Sigma, lambda = lambda, mu = mu,
       l_icept = l_icept, m_icept = m_icept, beta = beta)
}

# eigenvalues of the "covariance" matrix
evals_Sigma <- function(Sigma) {
  
  z_mc <- z_2d(k, 10^5, Sigma = th_s$Sigma)
  eigen(t(z_mc) %*% z_mc / nrow(z_mc))$values
}

perturb_th <- function(th, size = 1) {
  
  for (cmp in c("lambda", "mu", "l_icept", "m_icept", "beta")) {
    
    th[[cmp]] <- th[[cmp]] + size * runif(length(th[[cmp]]), -1, 1)
  }
  
  th
}
