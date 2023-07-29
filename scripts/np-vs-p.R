
k <- 5
nsamp <- 500
fi <- 0.1

res <- NULL
for (i in 1:100) {
  
  R <- r_2d(k, nsamp, fi = fi, seed = i)
  
  # get the ground truth
  pz_true <- attr(R, "pz")
  pz_true <- pz_true / sum(pz_true)
  pz_true
  
  # non-parametric, empirical estimator
  enc_mat <- t(replicate(nrow(R), 2^(seq_len(k) - 1L)))
  cprxy <- tabulate(
    rowSums(
      R * enc_mat
    ) + 1
  )
  cprxy <- c(cprxy, rep(0L, 2^k - length(cprxy))) # zero-padding
  cprxy <- cprxy / sum(cprxy)
  
  pz_np <- vapply(
    seq_len(2^k),
    function(z) infer_pz(z, cprxy, fi, k),
    numeric(1L)
  )
  
  # parametric estimator
  Sigma_zero <- array(0, dim = c(k, k))
  
  fi_mat <- array((1-fi)^2, dim = c(k, k))
  diag(fi_mat) <- 1-fi
  
  mu_hat <- t(R) %*% R / nrow(R) / fi_mat
  Sigma_hat <- mom_gradient(Sigma_zero, mu_hat, -1, nsamp)
  pz_par <- cpp_scaling_2d(Sigma_hat)
  pz_par <- pz_par / sum(pz_par)
  pz_par
  
  res <- rbind(``
    res,
    data.frame(method = "non-par", loss = vec_norm(pz_true - pz_np)),
    data.frame(method = "par", loss = vec_norm(pz_true - pz_par))
  )
}

ggplot(res, aes(x = loss, fill = method)) +
  geom_density(alpha = 0.3) + theme_bw()


