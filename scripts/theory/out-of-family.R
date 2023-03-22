

root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

n <- 5000
nrep <- 10
kseq <- 2:5
fiseq <- c(0, 0.01, 0.025, 0.05, 0.1, 0.2)

res <- expand.grid(k = kseq, fi = fiseq, rep = seq_len(nrep))
res$err_in <- res$err_out <- res$n <- 0

for (i in seq_len(nrow(res))) {
  
  k <- res$k[i]
  fi <- res$fi[i]
  seed <- res$rep[i]
  
  fixy <- list(list(0, 0), list(0, 0))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- fi
  dat <- synth_data_mid(n, k, fi = fixy, class = "expfam-2d", seed = seed)
  R <- dat$R
  
  pz <- attr(R, "pz")
  pz <- pz / sum(pz)
  
  z_all <- do.call(rbind, lapply(seq_len(2^k), function(zi) to_bits(zi-1, k = k)))
  lambda <- attr(R, "lambda")
  mu <- attr(R, "mu")
  
  x_ad <- expit(as.vector(z_all %*% lambda - (sum(lambda) / 2)))
  y_ad <- expit(as.vector(z_all %*% mu - (sum(mu) / 2)))
  
  pz_x0y0 <- pz * x_ad * y_ad
  pz_x0y0 <- pz_x0y0 / sum(pz_x0y0)
  
  Rx0y0 <- R[dat$X == 0 & dat$Y == 0, ]
  
  de_oof <- descent_2d(cov(Rx0y0), Rx0y0, fi = fi, n_iter = 20, 
                       type = "Newton-Raphson")
  
  Sigma_hat <- tail(de_oof$Sigma, n = 1L)[[1L]]
  
  pz_hat <- cpp_scaling_2d(Sigma_hat)
  pz_hat <- pz_hat / sum(pz_hat)
  err_oof <- vec_norm(pz_x0y0 - pz_hat)
  
  Rtx0y0 <- r_2d(k, nrow(Rx0y0), fi = fi, Sigma = Sigma_hat)
  
  de_itf <- descent_2d(cov(Rtx0y0), Rtx0y0, fi = fi, n_iter = 20, 
                       type = "Newton-Raphson")
  
  Sigma_that <- tail(de_itf$Sigma, n = 1L)[[1L]]
  pz_that <- cpp_scaling_2d(Sigma_that)
  pz_that <- pz_that / sum(pz_that)
  err_itf <- vec_norm(pz_hat - pz_that)
  
  res$err_in[i] <- err_itf
  res$err_out[i] <- err_oof
  res$n[i] <- nrow(Rx0y0)
}

ggplot(reshape2::melt(res, id.vars = names(res)[1:4]), 
       aes(x = fi, y = value, color = variable)) +
  geom_point() + geom_smooth() +
  facet_grid(cols = vars(k)) + theme_bw()

