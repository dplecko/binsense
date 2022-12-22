
library(ggplot2)

z_1d <- function(k = 3, n_samp, seed = 2022) {
  
  set.seed(seed)
  
  lambda <- runif(k, -1, 1)
  
  pz <- scaling_1d(lambda)
  
  pz <- pz / sum(pz)
  
  # can speed up with binary search? -> not really a bottleneck
  z <- sample.int(2^k, size = n_samp, prob = pz, replace = TRUE) 
  
  z <- do.call(rbind, lapply(z, function(zi) to_bits(zi-1, k = k)))
  attr(z, "lambda") <- lambda
  
  z
}

scaling_1d <- function(lambda) {
  
  k <- length(lambda)
  pz <- vapply(
    seq_len(2^k), 
    function(x) exp(sum(to_bits(x-1, k) * lambda)),
    numeric(1L)
  )
  
  pz
}

r_1d <- function(k, n_samp, fi, seed = 2022) {
  
  z <- z_1d(k, n_samp, seed)
  
  mask <- array(rbinom(k * n_samp, 1, 1-fi), dim = dim(z))
  
  r <- z * mask
  
  attr(r, "lambda") <- attr(z, "lambda")
  
  r
}

# check sampler correctness
# fi <- 0.3
# r <- r_1d(7, 10^4, fi)
# colMeans(r) - expit(attr(r, "lambda")) * (1-fi)

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
                             true_lambda = NULL) {
  
  lambda_dist <- loglc <- rep(0, n_iter)
  lambda <- list()
  lambda[[1]] <- lambda0
  loglc[1] <- likelihood_1d(lambda[[1]], r, fi)
  lambda_dist[1] <- sum((lambda[[1]] - true_lambda)^2)
  for (i in seq.int(2, n_iter)) {
    
    gamma_i <- gamma / (1 + i)
    lambda[[i]] <- lambda[[i-1]] + gamma_i * gradient_1d(lambda[[i-1]], r, fi)
    
    loglc[i] <- likelihood_1d(lambda[[i]], r, fi)
    lambda_dist[i] <- sum((lambda[[i]] - true_lambda)^2)
    cat("Log(L_c) =", loglc[i], "\n")
  }
  
  list(lambda = lambda, lambda_dist = lambda_dist, loglc = loglc)
}


k <- 10
fi <- 0.3
n_samp <- 10000
r <- r_1d(k, n_samp, fi, seed = 22)
true_lambda <- attr(r, "lambda")

em <- gradient_descent(rep(0, ncol(r)), r, fi = fi, true_lambda = true_lambda,
                       gamma = 10, n_iter = 100)

df <- data.frame(lambda_dist = em$lambda_dist, loglc = em$loglc,
                 n_iter = seq_along(em$loglc))

p1 <- ggplot(df, aes(x = n_iter, y = loglc)) +
  geom_line() + theme_bw()

p2 <- ggplot(df, aes(x = n_iter, y = lambda_dist)) +
  geom_line() + theme_bw() # +
  # geom_hline(yintercept = likelihood_1d(true_lambda, r, fi), 
  #            color = "red", linetype = "dashed", linewidth = 1)

cowplot::plot_grid(p1, p2, ncol = 1L)


#' * check if likelihood is correct *

compute_marginal_grad <- function(lambda, r, fi) {
  
  ri <- colMeans(r)
  
  ri + (1-ri) * expit(log(fi) + lambda) - expit(lambda)
}

llc <- rep(0, 100)
lam_seq <- seq(-3, 3, length.out = 100)

for (i in seq_along(llc)) {
  
  llc[i] <- compute_marginal_grad(lam_seq[i], fi = 0.3, r)
}


plot(lam_seq, llc, pch = 19)
abline(v = true_lambda, col = "red")
abline(h = 0, col = "blue")

# true_lambda_ext <- c(true_lambda, -sum(true_lambda))
# 
# df$llx <- 0
# for (i in 1:100) {
#   df$llx[i] <- likelihood_1d(NULL, r, fi, lambda_ext = true_lambda_ext - i)
# }
# 
# ggplot(df, aes(x = n_iter, y = llx)) +
#   geom_line() + theme_bw() +
#   geom_hline(yintercept = likelihood_1d(NULL, r, fi, lambda_ext = true_lambda_ext), 
#              color = "red", linetype = "dashed", linewidth = 1)
# 
# likelihood_1d(NULL, r, fi, lambda_ext = true_lambda_ext)
# 
# likelihood_1d(NULL, r, fi, lambda_ext = true_lambda_ext - 1)
# 
# sum(scaling_1d(true_lambda_ext)) / sum(scaling_1d(true_lambda_ext - 1))
# 
# 
# 
# 
# log_llc <- function(r, lambda, fi) {
# 
#   t1 <- log( (1-fi) * exp(lambda )) * mean(r) - mean(r) * log(1 + exp(lambda))
# 
#   q1 <-  (1-mean(r)) * (
#     fi * exp(lambda) / (1 + fi * exp(lambda)) * 
#       log( (exp(lambda) * fi) / (1 + exp(lambda)) ) +
#       1 / (1 + fi * exp(lambda)) * log( 1 / (1 + exp(lambda)) )
#   )
# 
#   q2 <- (1-mean(r)) * (
#     fi * exp(lambda) / (1 + fi * exp(lambda)) * 
#       log( fi * exp(lambda) / (1 + fi * exp(lambda) ) ) +
#         1 / (1 + fi * exp(lambda)) * log( 1 / (1 + fi * exp(lambda)) )
#   )
# 
#   t1 + q1  - q2
# }
# 
# fi <- 0.3
# lam1seq <- seq(-2, 2, 0.01)
# r1 <- rbinom(10^4, 1, expit(true_lambda[1])) * rbinom(10^4, 1, 1 - fi)
# llc <- vapply(lam1seq, function(lam) log_llc(r1, lam, fi), numeric(1L))
# df <- data.frame(lambda = lam1seq, llc = llc)
# 
# ggplot(df, aes(x = lambda, y = llc)) +
#   geom_line() + theme_bw() +
#   geom_vline(xintercept = true_lambda[1], 
#              color = "red", linetype = "dashed", linewidth = 1)
# 

