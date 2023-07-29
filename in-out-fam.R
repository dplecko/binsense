
library(ggplot2)

# generate the R
k <- 5
n <- 1750
fi <- 0.1

res <- NULL
for (i in 1:71) {
  
  seed <- i
  fixy <- list(list(0, 0), list(0, 0))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- fi
  dat <- synth_data_mid(n, k, fi = fixy, class = "expfam-2d", seed = seed)
  
  idx_x1y0 <- dat$X == 1 & dat$Y == 0
  
  n_eff <- sum(idx_x1y0)
  
  # take R | x, y
  r_x1y0 <- dat$R[idx_x1y0, ]
  
  # take equally many R samples
  r_rand <- dat$R[sample.int(nrow(dat$R), n_eff), ]
  
  # intial point Sigma matrix
  Sigma_zero <- array(0, dim = c(k, k))
  
  # apply moments method to both
  fi_mat <- array((1-fi)^2, dim = c(k, k))
  diag(fi_mat) <- 1-fi
  
  mu_hat_x1y0 <- t(r_x1y0) %*% r_x1y0 / nrow(r_x1y0) / fi_mat
  Sigma_hat_x1y0 <- mom_gradient(Sigma_zero, mu_hat_x1y0, -1, nrow(r_x1y0))
  
  mu_hat_rand <- t(r_rand) %*% r_rand / nrow(r_rand) / fi_mat
  Sigma_hat_rand <- mom_gradient(Sigma_zero, mu_hat_rand, -1, nrow(r_x1y0))
  
  # compute fitted values from fits
  pz_rand <- cpp_scaling_2d(Sigma_hat_rand)
  pz_rand <- pz_rand / sum(pz_rand)
  
  pz_x1y0 <- cpp_scaling_2d(Sigma_hat_x1y0)
  pz_x1y0 <- pz_x1y0 / sum(pz_x1y0)
  
  # extract theoretical values and compute L2 error
  
  err_x1y0 <- sum((pz_x1y0 - dat$pz_all[[2]][[1]])^2)
  err_rand <- sum((pz_rand - dat$pz)^2)
  
  res <- rbind(res, c(seed, err_x1y0, err_rand))
  cat("\r", i)
}

res <- as.data.frame(res)
names(res) <- c("seed", "out", "in")

ggplot(
  melt(res, id.vars = "seed"),
  aes(x = value, fill = variable)
) +
  geom_density(alpha = 0.4) + theme_minimal() + ylim(c(0, 500))
