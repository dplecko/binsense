
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

kseq <- c(2, 3, 4, 5, 6)
nseq <- c(100, 500, 1000, 2000, 5000)
fiseq <- c(0, 0.05, 0.1, 0.15, 0.2, 0.3)
nrep <- 100

res <- expand.grid(k = kseq, n = nseq, fi = fiseq, iter = seq_len(nrep))
res$err <- 0
res <- as.data.table(res)

for (i in seq_len(nrow(res))) {
  
  k <- res$k[i]
  n <- res$n[i]
  fi <- res$fi[i]
  
  rhash_all <- seq_len(2^k)
  
  # get true P(r)
  
  # get \widehat{P}(r)
  r_samp <- r_2d(k, n, fi, seed = i)
  rhash_samp <- 1 + rowSums(r_samp * t(replicate(nrow(r_samp), 2^(seq_len(k) - 1))))
  pz <- attr(r_samp, "pz")
  pz <- pz / sum(pz)
  
  A <- do.call(rbind, lapply(rhash_all, function(ri) cpp_fi_hash(ri, fi, k)))
  pr <- as.vector(A %*% pz)
  
  # get P(z | r)P(r)
  ex_zr <- cpp_grad_ij(rhash_all, 1, 2, k, fi, pz)
  
  pop <- sum(ex_zr * pr)
  smp <- mean(ex_zr[rhash_samp])
  
  res$err[i] <- (pop - smp)^2
  
  cat("\r", i, "finished")
}

# should have O(1/n), not care about k, not care about fi.
ggplot(res, aes(x = log(n), y = err * n)) +
  geom_point() + theme_bw() +
  facet_grid(rows = vars(k), cols = vars(fi)) +
  geom_smooth()

# 2nd case (population level, Term T2)
kseq <- c(2, 3, 4, 5, 6, 7, 8)
fiseq <- c(0, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3)
nrep <- 100

res <- expand.grid(k = kseq, fi = fiseq, iter = seq_len(nrep))
res$gamma <- res$dist <- 0
res <- as.data.table(res)

for (i in seq_len(nrow(res))) {
  
  k <- res$k[i]
  fi <- res$fi[i]
  
  rhash_all <- seq_len(2^k)
  
  # get true P(r)
  
  # get \widehat{P}(r)
  set.seed(i)
  
  Sigma_s <- gen_Sigma(k)
  pz_s <- cpp_scaling_2d(Sigma_s)
  pz_s <- pz_s / sum(pz_s)
  
  Sigma_c <- gen_Sigma(k)
  pz_c <- cpp_scaling_2d(Sigma_c)
  pz_c <- pz_c / sum(pz_c)
  
  A <- do.call(rbind, lapply(rhash_all, function(ri) cpp_fi_hash(ri, fi, k)))
  pr <- as.vector(A %*% pz_s)
  
  # get P(z | r)P(r)
  ex_zrts <- cpp_grad_ij(rhash_all, 1, 2, k, fi, pz_s)
  ex_zrtc <- cpp_grad_ij(rhash_all, 1, 2, k, fi, pz_c)
  
  res$gamma[i] <- sum( (ex_zrts - ex_zrtc) * pr )^2 / 
                  vec_norm(Sigma_s - Sigma_c)^2
  
  res$dist[i] <- vec_norm(Sigma_s - Sigma_c)
  
  res$gamma[i] <- sqrt(res$gamma[i])
  cat("\r", i, "finished")
}

# should have O(1/n), not care about k, and inspect fi?
ggplot(res, aes(x = k, y = gamma)) +
  geom_point() + theme_bw() +
  facet_grid(rows = vars(fi)) +
  geom_smooth()

ggplot(res, aes(x = gamma, y = dist)) +
  geom_point() + theme_bw() + 
  facet_grid(rows = vars(k), cols = vars(fi))


#' * Case 3: * population level, Term T2, sample paths
kseq <- 2:6
fiseq <- c(0, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3)
deltaseq <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5)
nrep <- 10


Sigma_norm <- function(Sigma) vec_norm(Sigma[!lower.tri(Sigma)])

res <- expand.grid(k = kseq, fi = fiseq, delta = deltaseq, 
                   iter = seq_len(nrep))
res$gamma <- res$dist <- 0
res <- as.data.table(res)

for (i in seq_len(nrow(res))) {
  
  k <- res$k[i]
  fi <- res$fi[i]
  seed <- res$iter[i]
  delta <- res$delta[i]
  
  rhash_all <- seq_len(2^k)
  
  # get true P(r)
  
  # get \widehat{P}(r)
  set.seed(seed)
  
  # get Sigma
  Sigma_s <- gen_Sigma(k)
  pz_s <- cpp_scaling_2d(Sigma_s)
  pz_s <- pz_s / sum(pz_s)
  
  # get perturbation
  delta_Sigma <- gen_Sigma(k)
  delta_Sigma <- delta_Sigma / Sigma_norm(delta_Sigma)
  
  Sigma_c <- Sigma_s + delta * delta_Sigma
  pz_c <- cpp_scaling_2d(Sigma_c)
  pz_c <- pz_c / sum(pz_c)
  
  A <- do.call(rbind, lapply(rhash_all, function(ri) cpp_fi_hash(ri, fi, k)))
  pr <- as.vector(A %*% pz_s)
  
  # get P(z | r)P(r)
  
  grad_diff <- Sigma_s
  grad_diff[,] <- 0
  for (j in seq.int(1, k)) {
    
    for (l in seq.int(j, k)) {
      
      ex_zrts <- cpp_grad_ij(rhash_all, j, l, k, fi, pz_s)
      ex_zrtc <- cpp_grad_ij(rhash_all, j, l, k, fi, pz_c)
      grad_diff[j, l] <- sum( (ex_zrts - ex_zrtc) * pr )
    }
  }
  
  res$gamma[i] <- Sigma_norm(grad_diff) / Sigma_norm(Sigma_s - Sigma_c)
  
  res$dist[i] <- Sigma_norm(Sigma_s - Sigma_c)
  
  # res$gamma[i] <- sqrt(res$gamma[i])
  cat("\r", i, "finished")
}

# should have O(1/n), not care about k, and inspect fi?
ggplot(res[delta < 5], aes(x = delta, y = gamma, color = factor(iter))) +
  geom_point() + geom_line() + theme_bw() +
  facet_grid(rows = vars(fi), cols = vars(k))
