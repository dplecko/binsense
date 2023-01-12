
gradient_1d <- function(lambda, r, fi) {
  
  k <- ncol(r)
  ri <- colMeans(r)
  
  ri + (1-ri) * expit(lambda + log(fi)) - expit(lambda) 
}

likelihood_1d <- function(lambda, r, fi, lambda_ext = NULL) {
  
  k <- ncol(r)
  ri <- colMeans(r)
  
  sum(
    ri * log(expit(lambda) * (1-fi)) + 
      (1-ri) * log( fi * expit(lambda) + 1 - expit(lambda))
  )
}

descent_1d <- function(lambda0, r, fi, gamma = 1, n_iter = 50,
                       true_lambda = NULL, verbose = FALSE) {
  
  lambda_dist <- loglc <- rep(0, n_iter)
  lambda <- list()
  lambda[[1]] <- lambda0
  loglc[1] <- likelihood_1d(lambda[[1]], r, fi)
  
  if (!is.null(true_lambda)) {
    
    lambda_dist[1] <- sum((lambda[[1]] - true_lambda)^2)
  }
  
  for (i in seq.int(2, n_iter)) {
    
    gamma_i <- gamma / (1 + i)
    lambda[[i]] <- lambda[[i-1]] + gamma_i * gradient_1d(lambda[[i-1]], r, fi)
    
    loglc[i] <- likelihood_1d(lambda[[i]], r, fi)
    
    if (!is.null(true_lambda)) {
      
      lambda_dist[i] <- sum((lambda[[i]] - true_lambda)^2)
    }
    
    if (verbose) cat("Log(L_c) =", loglc[i], "\n")
  }
  
  list(lambda = lambda, lambda_dist = lambda_dist, loglc = loglc)
}
