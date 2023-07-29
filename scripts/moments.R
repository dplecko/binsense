
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

k <- 4
fi <- 0.1
r <- r_2d(k, 10^4, fi, seed = 2029)

fi_mat <- array((1-fi)^2, dim = c(k, k))
diag(fi_mat) <- 1-fi
mu_hat <- t(r) %*% r / nrow(r) / fi_mat

Sigma_zero <- array(0, dim = c(k, k))

Sigma_hat <- mom_gradient(Sigma_zero, mu_hat, -1)

mom_obj(Sigma_hat, mu_hat) * 10^5

mom_obj(attr(r, "Sigma"), mu_hat) * 10^5

seed_seq <- 2023:2032
kseq <- 2:7
grid <- expand.grid(seed = seed_seq, k = kseq)
grid$fluct <- 0

for (i in seq_len(nrow(grid))) {
  
  k <- grid$k[i]
  seed <- grid$seed[i]
  fi <- 0.1
  r <- r_2d(k, 10^4, fi, seed = seed)
  
  fi_mat <- array((1-fi)^2, dim = c(k, k))
  diag(fi_mat) <- 1-fi
  mu_hat <- t(r) %*% r / nrow(r) / fi_mat
  
  grid$fluct[i] <- mom_obj(attr(r, "Sigma"), mu_hat) * 10^5
}

ggplot(grid, aes(x = k, y = fluct / k^2)) +
  geom_point() + geom_smooth() + theme_bw()


