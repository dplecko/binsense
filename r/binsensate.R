
#' Binary Sensitivity Analysis 
#'
#' 
#'
#' @param X The treatment variable (vector with 0/1 entries).
#' @param Y The outcome variable (vector with 0/1 entries).
#' @param R The set of proxies for the confounding variables Z (matrix with 0/1 
#' entries).
#' @param fi Sensitivity parameters of the analysis. List with two entries, both 
#' of which are lists. The entry `fi[[1+x]][[1+y]]` contains the sensitivity 
#' parameter corresponding to the probability that a confounder is unobserved for
#' the specific value of X = x, Y = y.
#' @param W The set of correctly observed continuous confounders, which should be
#' mean centered and variance normalized. Defaults to `NULL`, meaning no continuous
#' confounders are included. If `W` is included, the continous EM routine is used
#' automatically (`W` cannot be used with the \code{\link{backward_solver}}).
#' @param method Character taking values `"IF"` for independent fidelity, or
#' `"ZINF"` for the zero-inflation fidelity pattern.
#' @param A Sensitivity parameters that  List with two entries, both 
#' of which are lists. The entry `A[[1+x]][[1+y]]` contains the entries of the 
#' fidelity matrix A_xy for X = x, Y = y. Default to NULL (in which case) `fi` is
#' used to infer the value of A_xy.
#'@param solver Character scalar describing the solution method, 
#' either `"backward"` (non-parametric approach) or `"em"` (parametric EM).
#' @param ... Further arguments passed to subroutines \code{\link{em_solver}} or
#' \code{\link{em_solver_cts}}.
#' @importFrom Rcpp sourceCpp
#' @useDynLib binsensate
#' @export
binsensate <- function(X, Y, R, fi, W = NULL, method = c("IF", "ZINF"), A = NULL,
                       solver = c("backward", "em"), ...) {
  
  solver <- match.arg(solver, c("backward", "em"))
  method <- match.arg(method, c("IF", "ZINF"))
  
  if (!is.null(A)) {
    
    fi <- NULL
    method <- "matrix"
    if (is.matrix(A)) A <- list(list(A, A), list(A, A))
  } else {
    
    # if agnostic fi, create a list
    if (!is.list(fi)) fi <- list(list(fi, fi), list(fi, fi))
    fi <- expand_fi(fi, ncol(R), method)
    A <- fi_to_A(fi, ncol(R), method)
  }
  
  if (method == "ZINF") zinf_verify(X, Y, R, fi)
  
  if (!is.null(W)) {
    
    if (solver == "backward")
      message("Switching to `em` solver since `W` argument is provided.\n")
    solver <- "em-cts"
  }

  switch(
    solver,
    backward = backward_solver(X, Y, R, fi, method, A),
    em = em_solver(X, Y, R, fi, method, A, ...),
    `em-cts` = em_solver_cts(X, Y, R, W, fi, method, A, ...)
  )
}

#' Non-Parametric Binary Sensitivity Analysis
#' 
#' @inheritParams binsensate
#' @export
backward_solver <- function(X, Y, R, fi, method, A) {
  
  k <- ncol(R)
  enc_mat <- t(replicate(nrow(R), 2^(seq_len(k) - 1L)))
  fi_change <- attr(fi, "fi_change")
  if (is.null(fi_change)) fi_change <- FALSE
  
  pr <- pz <- pxy <- list(list(0, 0), list(0, 0))
  
  x0 <- y0 <- 0
  x1 <- y1 <- 1
  
  for (x in c(T, F)) {
    
    for (y in c(T, F)) {
      
      idx <- X == x & Y == y
      
      cprxy <- tabulate(
        rowSums(
          R[idx, , drop = FALSE] * enc_mat[idx, , drop = FALSE]
        ) + 1
      )
      cprxy <- c(cprxy, rep(0L, 2^k - length(cprxy))) # zero-padding
      cprxy <- cprxy / sum(cprxy)
      
      # get P(x, y) weights
      pxy[[x + 1]][[y + 1]] <- mean(idx)
      
      # get the bin-counting P(r | x, y) estimator
      pr[[x + 1]][[y + 1]] <- cprxy
      
      # obtain the inverse of the matrix
      if (is.element(method, c("IF", "ZINF"))) {
        
        A_inv_xy <- fi_to_Ainv_xy(fi[[1+x]][[1+y]], k, method)
      } else A_inv_xy <- solve(A[[1+x]][[1+y]])
      
      # obtain P(z | x,y) estimator using matrix inversion
      pz[[x + 1]][[y + 1]] <- A_inv_xy %*% pr[[x + 1]][[y + 1]]
      
      # perform simplex check for IF family
      if (check_simplex(pz[[x + 1]][[y + 1]])) {
        
        initial <- pz[[x + 1]][[y + 1]]
        
        # check if ill-posedness is statistically significant
        if (ill_pose_sig(pz[[x + 1]][[y + 1]], pr[[x + 1]][[y + 1]], 
                         A_inv_xy, nsamp = sum(idx))) {
          
          # update parameters to more appropriate values for IF method
          if (method == "IF") { 
            
            # trade for a well-posed inverse problem
            fi_new <- trade_inv_prob_IF(pz[[x + 1]][[y + 1]], 
                                        pr[[x + 1]][[y + 1]], 
                                        A_inv_xy, fi[[x + 1]][[y + 1]])
            
            # re-run with new fi
            fi_new <- fi_trade(fi, fi_new, x, y)
            Anew <- fi_to_A(fi_new, k, method)
            attr(fi_new, "fi_change") <- TRUE
            return(
              backward_solver(X, Y, R, fi_new, method, Anew)
            )
          } else if (method == "ZINF") {
            
            fi_new <- trade_inv_prob_ZINF(pz[[x + 1]][[y + 1]], 
                                          pr[[x + 1]][[y + 1]], 
                                          fi[[x + 1]][[y + 1]])
            
            # re-run with new fi
            fi_new <- fi_trade(fi, fi_new, x, y)
            Anew <- fi_to_A(fi_new, k, method)
            attr(fi_new, "fi_change") <- TRUE
            return(
              backward_solver(X, Y, R, fi_new, method, Anew)
            )
            
          } else { # for other methods
            
            stop(
              paste(
                "Sensitivity parameter likely not compatible with observed data.\n",
                " Try adjusting the parameter to a smaller value."
              )
            )
          }
        } else {
          
          neg_idx <- pz[[x + 1]][[y + 1]] < 0
          pz[[x + 1]][[y + 1]][neg_idx] <- 0
          pz[[x + 1]][[y + 1]] <- pz[[x + 1]][[y + 1]] / sum(pz[[x + 1]][[y + 1]])
        }
      }
    }
  }
  
  if (fi_change) 
    message("fi parameter was updated to be compatible with observed data.\n")
  
  # get the marginal P(z)
  pz_inf <- rep(0, 2^k)
  for (x in c(T, F)) {
    for (y in c(T, F)) {
      pz_inf <- pz_inf + pxy[[x + 1]][[y + 1]] * pz[[x + 1]][[y + 1]]
    }
  }
  
  # compute P(y | x_1, z) (vectorized)
  py_x1z <- pz[[x1 + 1]][[y1 + 1]] * pxy[[x1 + 1]][[y1 + 1]] / (
    pz[[x1 + 1]][[y1 + 1]] * pxy[[x1 + 1]][[y1 + 1]] + 
      pz[[x1 + 1]][[y0 + 1]] * pxy[[x1 + 1]][[y0 + 1]]
  )
  
  # compute P(y | x_0, z) (vectorized)
  py_x0z <- pz[[x0 + 1]][[y1 + 1]] * pxy[[x0 + 1]][[y1 + 1]] / (
    pz[[x0 + 1]][[y1 + 1]] * pxy[[x0 + 1]][[y1 + 1]] + 
      pz[[x0 + 1]][[y0 + 1]] * pxy[[x0 + 1]][[y0 + 1]]
  )
  
  py_x1z[is.nan(py_x1z)] <- 0
  py_x0z[is.nan(py_x0z)] <- 0
  ate <- sum((py_x1z - py_x0z) * pz_inf)
  
  return(
    list(pz = pz_inf, ATE = ate, fi = fi, solver = "backward", 
         fi_change = fi_change)
  )
}

#' Binary Sensitivity Analysis using Expectation-Maximization in an Exponential
#' Family
#' 
#' @inheritParams binsensate
#' @param mc_per_samp Number of Monte Carlo samples drawn for each observed 
#' sample in the data during the expectation-maximization algorithm.
#' @param n_epoch Number of epochs for the expectation-maximization algorithm.
#' @param rand_init Whether the expectation-maximization parameters should have
#' a random initialization. Defaults to `FALSE`.
#' @param se Logical indicating whether to compute standard errors using the
#' method of Louis. Defaults to `FALSE`.
#' @param detect_viol Logical indicating whether to test for ill-posedness of the
#' inverse problem. Defaults to `FALSE`.
#' @importFrom stats glm runif
#' @export
em_solver <- function(X, Y, R, fi, method, A, mc_per_samp = 10, n_epoch = 10, 
                      rand_init = FALSE, se = FALSE, detect_viol = FALSE) {
  
  # Step (0) - infer Sigma
  lmb <- Sigma <- list()
  
  k <- ncol(R)
  
  if (method == "IF") { # for IF method, infer Sigma using method of moments
    
    Sigma[[1]] <- infer_Sigma_IF(X, Y, R, fi)
  } else { # if method not IF, infer Sigma from the R correlation matrix
    
    cor_mat <- t(R) %*% R / nrow(R)
    Sigma[[1]] <- tail(cormat_to_Sigma(cor_mat), n = 1L)[[1]]
  }
  
  # get initial pz
  pz <- cpp_scaling_2d(Sigma[[1]])
  pz <- pz / sum(pz)
  
  # Step (0+) - infer initial lambda, mu, beta
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
         lambda_icept = lambda_icept, mu_icept = mu_icept)
  }
  
  if (!rand_init) {
    
    lmb[[1]] <- fit_lmb(X, Y, R)
  } else {
    
    lmb[[1]] <- list(lambda = runif(ncol(Sigma), -1, 1), 
                     mu = runif(ncol(Sigma), -1, 1), beta = runif(1, -1, 1), 
                     lambda_icept = runif(1, -1, 1), mu_icept = runif(1, -1, 1))
  }
  
  if (k > 5 | method == "ZINF") { # possibly consider a chunking approach
    
    # re-order
    enc <- 1 + rowSums(R * t(replicate(nrow(R), 2^(seq_len(k) - 1))))
    reord <- order(enc)
    R <- R[reord, ]
    X <- X[reord]
    Y <- Y[reord]
    enc <- enc[reord]
    
    chunks <- rle(enc)
    enc_vals <- chunks$values
    starts <- c(1, 1 + cumsum(chunks$lengths))
    starts <- starts[-length(starts)]
    stops <- cumsum(chunks$lengths)
    zero_idx <- which(enc == 1) # seq_len(stops[1])
    nzero_idx <- setdiff(seq_len(length(X)), zero_idx)
    
    by_chunk <- if ((length(enc_vals) / length(X)) < 0.2) TRUE else FALSE
  } else by_chunk <- FALSE
  
  for (i in seq_len(n_epoch)) {
    
    # get px_z, py_z
    px_z <- expit(cpp_p_z(lmb[[i]]$lambda) + lmb[[i]]$lambda_icept)
    logity_z <- cpp_p_z(lmb[[i]]$mu) + lmb[[i]]$mu_icept
    beta <- lmb[[i]]$beta
    
    if (by_chunk | method == "ZINF") {
      
      if (method == "ZINF") {
        
        mc_twist <- list()
        r_idx <- 1
        X_chunk <- X[zero_idx]
        Y_chunk <- Y[zero_idx]
        ret <- c()
        for (x in c(0, 1)) for (y in c(0, 1)) {
          
          num_xy <- sum(X_chunk == x & Y_chunk == y)
          if (num_xy == 0) next
          # P(r | z, x, y) conditionals
          pr_zxy <- A[[1 + x]][[1 + y]][r_idx, ]
          
          if (x == 1) px <- px_z else px <- 1 - px_z
          
          py_z <- expit(logity_z + beta * x)
          if (y == 1) py <- py_z else py <- 1 - py_z
          
          probs <- pz * px * py * pr_zxy
          probs <- probs / sum(probs)
          
          z_mc_idx <- sample(length(probs), size = num_xy * mc_per_samp, 
                             prob = probs, replace = TRUE)
          z_mc <- cpp_idx_to_bit(z_mc_idx - 1, k)
          ret <- rbind(ret, cbind(X = x, Y = y, Z = z_mc))
        }
        
        mc_twist[[1]] <- ret
        nzero_mc <- rep(nzero_idx, mc_per_samp)
        mc_twist[[2]] <- cbind(X = X[nzero_mc], Y = Y[nzero_mc], Z = R[nzero_mc, ])
        
      } else {
        
        mc_twist <- Map(
          function(r_idx, start, stop) {
            
            X_chunk <- X[start:stop]
            Y_chunk <- Y[start:stop]
            ret <- c()
            for (x in c(0, 1)) for (y in c(0, 1)) {
              
              num_xy <- sum(X_chunk == x & Y_chunk == y)
              if (num_xy == 0) next
              # P(r | z, x, y) conditionals
              pr_zxy <- A[[1 + x]][[1 + y]][r_idx, ]
              
              if (x == 1) px <- px_z else px <- 1 - px_z
              
              py_z <- expit(logity_z + beta * x)
              if (y == 1) py <- py_z else py <- 1 - py_z
              
              probs <- pz * px * py * pr_zxy
              probs <- probs / sum(probs)
              
              z_mc_idx <- sample(length(probs), size = num_xy * mc_per_samp, 
                                 prob = probs, replace = TRUE)
              
              z_mc <- cpp_idx_to_bit(z_mc_idx - 1, k)
              ret <- rbind(ret, cbind(X = x, Y = y, Z = z_mc))
            }

            ret
          },
          enc_vals, starts, stops
        )
      }
    } else {
      
      mc_twist <- lapply(
        seq_len(nrow(R)),
        function(j) {
          
          x <- X[j]
          y <- Y[j]
          r <- R[j, ]
          
          r_idx <- 1 + sum(r * 2^(seq_len(k)-1))
          
          # P(r | z, x, y) conditionals
          pr_zxy <- A[[1 + x]][[1 + y]][r_idx, ]
          
          if (x == 1) px <- px_z else px <- 1 - px_z
          
          # if (j %% 100 == 0) pb$tick()
          
          py_z <- expit(logity_z + beta * x)
          if (y == 1) py <- py_z else py <- 1 - py_z
          
          probs <- pz * px * py * pr_zxy
          probs <- probs / sum(probs)
          
          z_mc_idx <- sample(length(probs), size = mc_per_samp, prob = probs,
                             replace = TRUE)
          z_mc <- cpp_idx_to_bit(z_mc_idx - 1, length(r)) # returns a matrix
          
          cbind(X = x, Y = y, Z = z_mc)
        }
      )
    }
    mc_twist <- as.data.frame(do.call(rbind, mc_twist))
    
    # re-fit for Sigma
    Z_mc <- as.matrix(mc_twist[, -c(1, 2)])
    cor_mat <- t(Z_mc) %*% Z_mc / nrow(Z_mc)
    Sigma[[i+1]] <- tail(cormat_to_Sigma(cor_mat), n = 1L)[[1]]
    
    pz <- cpp_scaling_2d(Sigma[[i+1]])
    pz <- pz / sum(pz)
    
    # re-fit for lambda, mu, beta
    lmb[[i+1]] <- fit_lmb(mc_twist[, 1], mc_twist[, 2], mc_twist[, -c(1, 2)])
    
    # early stopping criterion
    scale <- 1
    converged <- vec_norm(lmb[[i+1]]$lambda - lmb[[i]]$lambda)^2 < 
      scale * ncol(R) / nrow(R) &
      vec_norm(lmb[[i+1]]$mu - lmb[[i]]$mu)^2 < scale * ncol(R) / nrow(R) &
      lmb[[i+1]]$beta - lmb[[i]]$beta < scale / nrow(R)
    if (converged) break
  }
  
  # Louis' method for Standard Errors (SEs)
  if (se) {
    
    # get the latest Monte Carlo sample
    Z_mc <- as.matrix(mc_twist[, -c(1, 2)])
    X_mc <- mc_twist[, 1]
    Y_mc <- mc_twist[, 2]
    
    Sigma_idx <- which(!upper.tri(Sigma[[i+1]]))
    
    # compute via Exponential families
    tX <- do.call(rbind, lapply(
      seq_len(nrow(Z_mc)),
      function(i) {
        
        unit_cor <- 2 * (Z_mc[i, ] %*% t(Z_mc[i, ]))
        diag(unit_cor) <- diag(unit_cor) / 2
        c(unit_cor[Sigma_idx], X_mc[i], Z_mc[i, ] * X_mc[i], Y_mc[i],
          Z_mc[i, ] * Y_mc[i], Y_mc[i] * X_mc[i])
      }
    ))
    
    covtX <- t(tX) %*% tX / nrow(tX) - colMeans(tX) %*% t(colMeans(tX))
    
    covtX_y <- array(0, dim = dim(covtX))
    for (i in seq_len(length(X))) {
      
      mc_idx <- seq.int((i-1) * mc_per_samp + 1, i * mc_per_samp)
      
      # Monte Carlo sample
      tXi <- t(tX[mc_idx, ])
      tXi_c <- tXi - rowMeans(tXi)
      
      covtX_y <- covtX_y + tXi_c %*% t(tXi_c) / mc_per_samp
    }
    
    covtX_y <- covtX_y / length(X)
    se_mat <- tryCatch(
      solve(covtX - covtX_y),
      error = function(e) e
    )
    if (inherits(se_mat, "error")) {
      
      message(paste(
        "Observed Fisher information matrix nearly singular.",
        "Louis' method cannot be used."
      ))
      louis_se <- NA
    } else {
      
      louis_se <- sqrt(se_mat[length(covtX - covtX_y)] / length(X))
      if (louis_se < 0) louis_se <- NA
    }
  } else louis_se <- NULL
  
  # detecting incompatibility of fi with the observed data
  if (detect_viol) {
    
    if (incompat_detect(X, Y, R, fi, A, lmb, Sigma)) {
      
      message("Inverse problem may be ill-posed with current fi parameter.\n")
      fi_change <- TRUE 
    } else fi_change <- FALSE
  } else fi_change <- FALSE
  
  # ATE computation
  tl <- length(lmb)
  px_z <- expit(cpp_p_z(lmb[[tl]]$lambda) + lmb[[tl]]$lambda_icept)
  logity_z <- cpp_p_z(lmb[[tl]]$mu) + lmb[[tl]]$mu_icept
  beta <- lmb[[tl]]$beta
  ATE_hat <- sum( (expit(logity_z + beta) - expit(logity_z)) * pz)
  if (!is.null(louis_se)) {
    
    ATE_lwr <- sum( 
      (expit(logity_z + (beta - 1.96 * louis_se)) - expit(logity_z)) * pz
    )
    ATE_upr <- sum( 
      (expit(logity_z + (beta + 1.96 * louis_se)) - expit(logity_z)) * pz
    )
  } else ATE_lwr <- ATE_upr <- NULL
  
  list(
    pz = pz,
    ATE = ATE_hat,
    ATE_lwr = ATE_lwr,
    ATE_upr = ATE_upr,
    solver = "two-stage-ecm",
    fi_change = fi_change,
    Sigma = Sigma,
    beta = beta,
    beta_se = louis_se,
    theta = lmb
  )
}

#' Binary Sensitivity Analysis using Expectation-Maximization in an Exponential
#' Family with Correctly Observed Continuous Confounders
#' 
#' @inheritParams binsensate
#' @inheritParams em_solver
#' @param mc_cores Number of cores used for parallelization during Monte Carlo
#' sampling in the M-step.
#' @importFrom parallel mclapply
#' @export
em_solver_cts <- function(X, Y, R, W, fi, method, A, 
                          mc_per_samp = 10, n_epoch = 10, 
                          rand_init = FALSE, se = FALSE, detect_viol = FALSE,
                          mc_cores = NULL) {
  
  # Step (0) - infer Sigma
  lmb <- theta <- Sigma <- list()
  k <- ncol(R)
  d <- ncol(W)
  k_index <- seq_len(k)
  d_index <- k + seq_len(d)
  
  cor_mat <- t(cbind(R, W)) %*% cbind(R, W) / nrow(R)
  theta[[1]] <- tail(cormat_to_Sigma_cts(cor_mat, k, d), n = 1L)[[1]]
  
  # extract Sigma, Lambda, Sigma0
  Sigma <- theta[[1]][seq_len(k), seq_len(k)]
  Lambda <- theta[[1]][-seq_len(k), seq_len(k), drop=FALSE]
  Omega <- theta[[1]][-seq_len(k), -seq_len(k)]
  
  # Step (0+) - infer initial lambda, mu, beta
  if (!rand_init) {
    
    lmb[[1]] <- fit_lmb_rw(X, Y, R, W)
  } else {
    
    lmb[[1]] <- list(lambda_z = runif(ncol(R), -1, 1), 
                     lambda_w = runif(ncol(W), -1, 1), 
                     mu_z = runif(ncol(R), -1, 1), mu_w = runif(ncol(W), -1, 1), 
                     beta = runif(1, -1, 1), 
                     lambda_icept = runif(1, -1, 1), mu_icept = runif(1, -1, 1))
  }
  
  # check if fi == 0; MC steps not needed
  A_id <- fi_to_A_xy(rep(0, ncol(R)), ncol(R), "IF")
  mfi <- 0
  for (x in c(1, 2)) for (y in c(1, 2)) mfi <- max(mfi, max(abs(A[[x]][[y]] - A_id)))
  if (mfi == 0) {
    n_epoch <- 0
    ZW_mc <- cbind(Z = R, W = W)
    se <- FALSE
  }
  
  final_run <- FALSE # flag for last iteration (for Z,W Monte Carlo)
  for (i in seq_len(n_epoch)) {
    
    # get px_z, py_z
    logitx_z <- cpp_p_z(lmb[[i]]$lambda_z) + lmb[[i]]$lambda_icept
    logity_z <- cpp_p_z(lmb[[i]]$mu_z) + lmb[[i]]$mu_icept
    beta <- lmb[[i]]$beta
  
    mc_fn <- function(j) {
      
      x <- X[j]
      y <- Y[j]
      r <- R[j, ]
      w <- W[j, ]
      
      # compute P(z, w)
      pzw <- cpp_pzw_prop(Sigma, t(Lambda) %*% w)
      
      r_idx <- 1 + sum(r * 2^(seq_len(k)-1))
      
      # P(r | z, x, y) conditionals
      pr_zxy <- A[[1 + x]][[1 + y]][r_idx, ]
      
      px_zw <- expit(logitx_z + sum(lmb[[i]]$lambda_w * w))
      if (x == 1) px <- px_zw else px <- 1 - px_zw
      
      # if (j %% 100 == 0) pb$tick()
      py_zw <- expit(logity_z + sum(lmb[[i]]$mu_w * w) + beta * x)
      if (y == 1) py <- py_zw else py <- 1 - py_zw
      
      probs <- pzw * px * py * pr_zxy
      probs <- probs / sum(probs)
      
      z_mc_idx <- sample(length(probs), size = mc_per_samp, prob = probs,
                         replace = TRUE)
      z_mc <- cpp_idx_to_bit(z_mc_idx - 1, length(r)) # returns a matrix
      
      cbind(X = x, Y = y, Z = z_mc, 
            W = matrix(w, ncol = d, nrow = mc_per_samp, byrow = TRUE))
    }
    
    if (!is.null(mc_cores)) {
      
      mc_twist <- mclapply(seq_len(nrow(R)), mc_fn, mc.cores = mc_cores)
    } else {
      
      mc_twist <- lapply(seq_len(nrow(R)), mc_fn)
    }
    mc_twist <- as.data.frame(do.call(rbind, mc_twist))
    
    # re-fit for Sigma
    ZW_mc <- as.matrix(mc_twist[, -c(1, 2)])
    cor_mat <- t(ZW_mc) %*% ZW_mc / nrow(ZW_mc)
    
    if (final_run) break
    
    theta[[i+1]] <- tail(cormat_to_Sigma_cts(cor_mat, k, d), n = 1L)[[1]]
    
    # re-fit for lambda, mu, beta
    lmb[[i+1]] <- fit_lmb_rw(X = mc_twist[, 1], Y = mc_twist[, 2], 
                             R = mc_twist[, 2 + k_index], 
                             W = mc_twist[, 2 + d_index])
    
    # early stopping criterion
    scale <- 1
    converged <- vec_norm(lmb[[i+1]]$lambda_z - lmb[[i]]$lambda_z)^2 < 
      scale * ncol(R) / nrow(R) &
      vec_norm(lmb[[i+1]]$mu_z - lmb[[i]]$mu_z)^2 < scale * ncol(R) / nrow(R) &
      lmb[[i+1]]$beta - lmb[[i]]$beta < scale / nrow(R)
    if (converged) final_run <- TRUE
  }
  
  # Louis' method for Standard Errors (SEs)
  if (se) {

    # get the latest Monte Carlo sample
    ZW_mc <- as.matrix(mc_twist[, -c(1, 2)])
    X_mc <- mc_twist[, 1]
    Y_mc <- mc_twist[, 2]

    theta_idx <- which(!upper.tri(theta[[i]]))

    # compute via Exponential families
    tX <- do.call(rbind, lapply(
      seq_len(nrow(ZW_mc)),
      function(i) {

        unit_cor <- 2 * (ZW_mc[i, ] %*% t(ZW_mc[i, ]))
        unit_cor[seq_len(k), seq_len(k)] <- 2 * unit_cor[seq_len(k), seq_len(k)]
        unit_cor[-seq_len(k), -seq_len(k)] <- -1 * unit_cor[-seq_len(k), -seq_len(k)]
        diag(unit_cor) <- diag(unit_cor) / 2

        c(unit_cor[theta_idx], X_mc[i], ZW_mc[i, ] * X_mc[i], Y_mc[i],
          ZW_mc[i, ] * Y_mc[i], Y_mc[i] * X_mc[i])
      }
    ))

    covtX <- t(tX) %*% tX / nrow(tX) - colMeans(tX) %*% t(colMeans(tX))

    covtX_y <- array(0, dim = dim(covtX))
    for (i in seq_len(length(X))) {

      mc_idx <- seq.int((i-1) * mc_per_samp + 1, i * mc_per_samp)

      # Monte Carlo sample
      tXi <- t(tX[mc_idx, , drop=FALSE])
      tXi_c <- tXi - rowMeans(tXi)

      covtX_y <- covtX_y + tXi_c %*% t(tXi_c) / mc_per_samp
    }

    covtX_y <- covtX_y / length(X)
    se_mat <- tryCatch(
      solve(covtX - covtX_y),
      error = function(e) e
    )
    
    if (inherits(se_mat, "error")) {

      message(paste(
        "Observed Fisher information matrix nearly singular.",
        "Louis' method cannot be used."
      ))
      louis_se <- NA
    } else {

      louis_se <- sqrt(se_mat[length(covtX - covtX_y)] / length(X))
      if (louis_se < 0) louis_se <- NA
    }
  } else louis_se <- NULL
  
  if (length(i) == 0) louis_se <- lmb[[1]]$beta_se
  
  # detecting incompatibility of fi with the observed data
  if (detect_viol) {
    
    if (incompat_detect(X, Y, R, fi, A, lmb, Sigma)) {
      
      message("Inverse problem may be ill-posed with current fi parameter.\n")
      fi_change <- TRUE 
    } else fi_change <- FALSE
  } else fi_change <- FALSE
  
  # ATE computation
  tl <- length(lmb)
  
  beta <- lmb[[tl]]$beta
  logity_zw <- cbind(1, ZW_mc) %*% c(lmb[[tl]]$mu_icept, lmb[[tl]]$mu_z, lmb[[tl]]$mu_w)
  ATE_hat <- mean( (expit(logity_zw + beta) - expit(logity_zw)))
  if (!is.null(louis_se)) {
    
    ATE_lwr <- mean( 
      (expit(logity_zw + (beta - 1.96 * louis_se)) - expit(logity_zw))
    )
    ATE_upr <- mean( 
      (expit(logity_zw + (beta + 1.96 * louis_se)) - expit(logity_zw))
    )
  } else ATE_lwr <- ATE_upr <- NULL
  
  list(
    # pz = pz,
    ATE = ATE_hat,
    ATE_lwr = ATE_lwr,
    ATE_upr = ATE_upr,
    solver = "two-stage-ecm",
    fi_change = fi_change,
    SLO = theta,
    beta = beta,
    beta_se = louis_se,
    LMB = lmb
  )
}
