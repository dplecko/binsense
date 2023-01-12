
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

experiment <- FALSE
if (experiment) {
  
  k <- 10
  fi <- 0.3
  n_samp <- 10000
  r <- r_2d(k, n_samp, fi, seed = 2022)
  true_Sigma <- attr(r, "Sigma")
  true_pz <- attr(r, "pz")
  
  # global \phi hash dictionary
  fi_wgh_hash <- list()
  
  em <- descent_2d(
    # tail(em$Sigma, n = 1L)[[1]], 
    array(0, dim = c(k, k)),
    r, fi = fi, 
    true_Sigma = true_Sigma,
    gamma = 2, n_iter = 5
  )
  
  df <- data.frame(Sigma_dist = em$Sigma_dist, loglc = em$loglc,
                   n_iter = seq_along(em$loglc))
  
  p1 <- ggplot(df, aes(x = n_iter, y = Sigma_dist)) +
    geom_line() + theme_bw() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed")
  
  p2 <- ggplot(df, aes(x = n_iter, y = loglc)) +
    geom_line() + theme_bw() +
    geom_hline(yintercept = likelihood_2d(true_Sigma, r, fi), 
               color = "red", linetype = "dashed")
  
  cowplot::plot_grid(p1, p2, ncol = 1L)
  
  # check that true_Sigma minimizes likelihood
  Sigma_hat <- tail(em$Sigma, n = 1L)[[1]]
  
  line_cmb <- function(A, B, lambda) lambda * A + (1 - lambda) * B
  
  loglc <- rep(0, 100)
  lambda <- seq(-5, 5, length.out = 100)
  for (i in seq_along(loglc)) {
    loglc[i] <- likelihood_2d(line_cmb(Sigma_hat, true_Sigma, lambda[i]), r, fi)
  }
  
  plot(lambda, loglc, pch = 19)
  abline(v = 0, col = "red")
  abline(v = 1, col = "blue")
}

