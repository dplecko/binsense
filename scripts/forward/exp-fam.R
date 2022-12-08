
to_bits <- function(i, k) as.integer(intToBits(i)[seq_len(k)])

exp_fam_gen <- function(k = 3, n_samp, seed = 2022) {
  
  set.seed(seed)
  
  lambda <- runif(k, -1, 1)
  lambda <- lambda - mean(lambda)
  
  pz <- scaling_fact(lambda)
  
  pz <- pz / sum(pz)
  
  # can speed up with binary search?
  z <- sample.int(2^k, size = n_samp, prob = pz, replace = TRUE) 
  
  z <- do.call(rbind, lapply(z, function(zi) to_bits(zi, k = k)))
  attr(z, "lambda") <- lambda
  
  z
}

scaling_fact <- function(lambda) {
  
  k <- length(lambda)
  pz <- vapply(
    seq_len(2^k), 
    function(x) exp(sum(to_bits(x-1, k) * lambda)),
    numeric(1L)
  )
}

r_dat_gen <- function(k, n_samp, fi, seed = 2022) {
  
  z <- exp_fam_gen(k, n_samp, seed)
  
  mask <- array(rbinom(k * n_samp, 1, 1-fi), dim = dim(z))
  
  r <- z * mask
  
  attr(r, "lambda") <- attr(z, "lambda")
  
  r
}

compute_gradient <- function(lambda, r, fi) {
  
  ri <- colMeans(r)
  
  grad1 <- ri + (1-ri) * (fi * exp(lambda)) / (1 + fi * exp(lambda))
  
  grad2 <- rep(0, length(lambda))
  for (i in seq_along(grad2)) {
    
    grad2[i] <- sum(scaling_fact(lambda[-i]))
  }
  
  grad2 <- grad2 / sum(scaling_fact(lambda))
  
  grad1 - grad2
}

compute_likelihood <- function(lambda, r, fi) {
  
  ri <- colMeans(r)
  sum(ri + (1-ri) * lambda * (fi * exp(lambda)) / (1 + fi * exp(lambda))) -
    log(sum(scaling_fact(lambda)))
}

gradient_descent <- function(lambda0, r, fi, gamma = 1, n_iter = 50,
                             true_lambda = NULL) {
  
  lambda_dist <- loglc <- rep(0, n_iter)
  lambda <- list()
  lambda[[1]] <- lambda0
  loglc[1] <- compute_likelihood(lambda[[1]], r, fi)
  lambda_dist[1] <- sum((lambda[[1]] - true_lambda)^2)
  for (i in seq.int(2, n_iter)) {
    
    gamma_i <- gamma / (10 + i)
    lambda[[i]] <- lambda[[i-1]] + gamma_i * compute_gradient(lambda[[i-1]], r, fi)
    
    loglc[i] <- compute_likelihood(lambda[[i]], r, fi)
    lambda_dist[i] <- sum((lambda[[i]] - true_lambda)^2)
    cat("Log(L_c) =", loglc[i], "\n")
  }
  
  list(lambda = lambda, lambda_dist = lambda_dist, loglc = loglc)
}


k <- 3
fi <- 0.3
n_samp <- 10
r <- r_dat_gen(k, n_samp, fi)
true_lambda <- attr(r, "lambda")

em <- gradient_descent(rep(0, ncol(r)), r, fi = fi, true_lambda = true_lambda,
                       gamma = 5, n_iter = 100)

df <- data.frame(lambda_dist = em$lambda_dist, loglc = em$loglc,
                 n_iter = seq_along(em$loglc))

ggplot(reshape2::melt(df, id.vars = "n_iter"), aes(x = n_iter, y = value)) +
  geom_line() + theme_bw() +
  facet_grid(rows = vars(variable), scales = "free")
