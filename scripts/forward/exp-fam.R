

root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), source))

k <- 10
fi <- 0.3
n_samp <- 10000
r <- r_1d(k, n_samp, fi, seed = 22)
true_lambda <- attr(r, "lambda")

em <- descent_1d(rep(0, ncol(r)), r, fi = fi, true_lambda = true_lambda,
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
