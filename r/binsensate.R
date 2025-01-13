
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
#' @param A Sensitivity parameters that  List with two entries, both 
#' of which are lists. The entry `A[[1+x]][[1+y]]` contains the entries of the 
#' fidelity matrix A_xy for X = x, Y = y. Default to NULL (in which case) `fi` is
#' used to infer the value of A_xy.
#'@param solver Character scalar describing the solution method, 
#' either `"backward"` (non-parametric approach) 
#' or `"two-stage-em"` (parametric EM).
#' @param ... Further arguments passed to subroutine \code{\link{two_stage_em}}.
#' @importFrom Rcpp sourceCpp
#' @useDynLib binsensate
#' @export
binsensate <- function(X, Y, R, fi, method = c("IM", "ZINF"), A = NULL,
                       solver = c("backward", "two-stage-em"), ...) {
  
  solver <- match.arg(solver, c("backward", "two-stage-em"))
  
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
  
  switch(
    solver,
    backward = backward_solver(X, Y, R, fi, method, A),
    `two-stage-em` = two_stage_em(X, Y, R, fi, method, A, ...)
  )
}

#' Non-Parametric Binary Sensitivity Analysis
#' 
#' @inheritParams binsensate
#' @export
backward_solver <- function(X, Y, R, fi, method, A) {
  
  k <- ncol(R)
  enc_mat <- t(replicate(nrow(R), 2^(seq_len(k) - 1L)))
  
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
      if (is.element(method, c("IM", "ZINF"))) {
        
        A_inv_xy <- fi_to_Ainv_xy(fi[[1+x]][[1+y]], k, method)
      } else A_inv_xy <- solve(A[[1+x]][[1+y]])
      
      # obtain P(z | x,y) estimator using matrix inversion
      pz[[x + 1]][[y + 1]] <- A_inv_xy %*% pr[[x + 1]][[y + 1]]
      
      # perform simplex check for IM family
      if (check_simplex(pz[[x + 1]][[y + 1]])) {
        
        initial <- pz[[x + 1]][[y + 1]]
        
        # check if ill-posedness is statistically significant
        if (ill_pose_sig(pz[[x + 1]][[y + 1]], pr[[x + 1]][[y + 1]], 
                         A_inv_xy, nsamp = sum(idx))) {
          
          # update parameters to more appropriate values for IM method
          if (method == "IM") { 
            
            # trade for a well-posed inverse problem
            fi_new <- trade_inv_prob(pz[[x + 1]][[y + 1]], 
                                     pr[[x + 1]][[y + 1]], 
                                     A_inv_xy, fi[[x + 1]][[y + 1]])
            
            # re-run with new fi
            fi_new <- fi_trade(fi, fi_new, x, y)
            Anew <- fi_to_A(fi_new, k, method)
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
    list(pz = pz_inf, ATE = ate, fi = fi, solver = "backward")
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
#' @importFrom stats glm runif
#' @export
two_stage_em <- function(X, Y, R, fi, method, A, mc_per_samp = 10, n_epoch = 10, 
                         rand_init = FALSE) {
  
  # Step (0) - infer Sigma
  lmb <- Sigma <- list()
  
  k <- ncol(R)
  
  if (method == "IM") { # for IM method, infer Sigma using method of moments
    
    Sigma[[1]] <- infer_Sigma_IM(X, Y, R, fi)
  } else { # if method not IM, infer Sigma from the R correlation matrix
    
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
  
  for (i in seq_len(n_epoch)) {
    
    # get px_z, py_z
    px_z <- expit(cpp_p_z(lmb[[i]]$lambda) + lmb[[i]]$lambda_icept)
    logity_z <- cpp_p_z(lmb[[i]]$mu) + lmb[[i]]$mu_icept
    beta <- lmb[[i]]$beta
    
    # sample mc_samples for each actual sample
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
    mc_twist <- as.data.frame(do.call(rbind, mc_twist))
    
    # re-fit for Sigma
    if (!(method == "IM")) {
      
      Z_mc <- as.matrix(mc_twist[, -c(1, 2)])
      cor_mat <- t(Z_mc) %*% Z_mc / nrow(Z_mc)
      Sigma[[i+1]] <- tail(cormat_to_Sigma(cor_mat), n = 1L)[[1]]
    } else { # in the case of independent missingess, do not need to update Sigma
      
      Sigma[[i+1]] <- Sigma[[i]]
    }
    
    pz <- cpp_scaling_2d(Sigma[[i+1]])
    pz <- pz / sum(pz)
    
    # re-fit for lambda, mu, beta
    lmb[[i+1]] <- fit_lmb(mc_twist[, 1], mc_twist[, 2], mc_twist[, -c(1, 2)])
    
    # early stopping criterion
    scale <- 1
    if (
      vec_norm(lmb[[i+1]]$lambda - lmb[[i]]$lambda)^2 < 
      scale * ncol(R) / nrow(R) &
      vec_norm(lmb[[i+1]]$mu - lmb[[i]]$mu)^2 < scale * ncol(R) / nrow(R) &
      lmb[[i+1]]$beta - lmb[[i]]$beta < scale / nrow(R)
    ) break
  }
  
  px_z <- expit(cpp_p_z(lmb[[i+1]]$lambda) + lmb[[i+1]]$lambda_icept)
  logity_z <- cpp_p_z(lmb[[i+1]]$mu) + lmb[[i+1]]$mu_icept
  beta <- lmb[[i+1]]$beta
  
  list(
    pz = pz,
    ATE = sum( (expit(logity_z + beta) - expit(logity_z)) * pz),
    solver = "two-stage-ecm",
    Sigma = Sigma,
    beta = beta,
    theta = lmb
  )
}
