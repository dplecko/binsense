
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r-utils")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE),
                 source))

k <- 5
th_s <- gen_params(k, rand = TRUE)
th_t <- gen_params(k, rand = FALSE)
cov_evals <- evals_Sigma(th_s$Sigma)

n_iter <- 25
fi <- 0.5
gamma <- gamma_p <- gamma_p1 <- gamma_p2 <- gamma_r <- dist <- rep(0, n_iter)
th_all <- list()

# what is the for loop iterating over?
for (i in seq_len(n_iter)) {
  
  # distance of current theta_t from true theta_s
  dist[i] <- vec_norm(em_params(th_t) - em_params(th_s))
  
  # find maximizer w.r.t. current theta_t
  M_th <- M_theta(th_t, th_s, fi)
  
  # get a random draw from space
  th_rand <- gen_params(k, rand = TRUE)
  
  # evaluated at a random point instead of M(theta_t)
  gamma[i] <- vec_norm(nabla_Q(M_th, th_t, th_s) - 
                       nabla_Q(M_th, th_s, th_s)) / 
              vec_norm(em_params(th_t) - em_params(th_s))
  
  gamma_r[i] <- vec_norm(nabla_Q(th_rand, th_t, th_s) - 
                         nabla_Q(th_rand, th_s, th_s)) / 
                vec_norm(em_params(th_t) - em_params(th_s))
  
  
  # gamma_p[i] <- vec_norm(nabla_Q(th_t, th_t, th_s) - nabla_Q(th_t, th_s, th_s)) / 
  #               vec_norm(em_params(th_t) - em_params(th_s))
  
  gamma_p1[i] <- vec_norm(nabla_Q(th_t, th_t, th_s)) / 
    vec_norm(em_params(th_t) - em_params(th_s))
  
  gamma_p2[i] <- vec_norm(nabla_Q(th_t, th_s, th_s)) / 
    vec_norm(em_params(th_t) - em_params(th_s))
  
  # simplified gamma
  # gamma_p[i] <- sum(diff_nabla_simp(th_rand, th_t, th_s)) / 
  #   vec_norm(em_params(th_t) - em_params(th_s))
  
  # cs_bII(th_t, M_th)
  
  th_t <- M_th
  th_all[[i]] <- th_t
}

# plotting of the results
df <- data.frame(gamma, gamma_p, gamma_r, gamma_p1, gamma_p2, dist, 
                 n_iter = seq_len(n_iter))
df <- reshape2::melt(df, id.vars = c("dist", "n_iter"))
p1 <- ggplot(df, aes(x = n_iter, y = dist)) +
  geom_point() + geom_line() + theme_bw()
p2 <- ggplot(df, aes(x = n_iter, y = value, color = variable)) +
  geom_point() + geom_line() + theme_bw() + 
  theme(legend.position = c(0.8, 0.8)) +
  geom_hline(yintercept = min(cov_evals), linetype = "dashed", color = "red")
cowplot::plot_grid(p1, p2, ncol = 2L)

n_rand <- 1000
gamma_rand <- rep(0, n_rand)

# what is the for loop iterating over?
for (i in seq_len(n_rand)) {
  
  th_rand <- gen_params(k, rand = TRUE)
  th_s_rand <- gen_params(k, rand = TRUE)
  
  gamma_rand[i] <- vec_norm(nabla_Q(th_rand, th_t, th_s_rand) - 
                            nabla_Q(th_rand, th_s, th_s_rand) ) / 
                   vec_norm(em_params(th_t) - em_params(th_s))
}

# project onto 2 dimensions? 

df <- data.frame(gamma, dist, n_iter = seq_len(n_iter))

# how far are theta_t and theta_s? how does \gamma vary with the distance,
# and with steps?
p1 <- ggplot(df, aes(x = n_iter, y = dist)) +
  geom_point() + geom_line() + theme_bw()

p2 <- ggplot(df, aes(x = n_iter, y = gamma)) +
  geom_point() + geom_line() + theme_bw() +
  geom_hline(yintercept = min(cov_evals), linetype = "dashed", color = "red") #  + 
  # geom_hline(yintercept = min(cov_evals))

# how does \gamma compare to the eigenvalues of E[Z Z^T]?
cowplot::plot_grid(p1, p2, ncol = 2L)

# check the CS Term (B-II) - 
cs_bII <- function(th, th_p) {
  
  Sigma <- th$Sigma
  k <- ncol(Sigma)
  
  pz <- binsensate:::cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  
  t_bII <- 0
  for (x in c(0, 1)) {
    
    pxc_z <- px_z(th)
    if (x == 0) pxc_z <- 1 - pxc_z
    pzx <- pz * pxc_z
    
    term_x <- (x - px_z(th_p))^2
    t_bII <- t_bII + sum(pzx * term_x)
  }
  
  t_bII
}

cs_bII(gen_params(k, rand = TRUE), gen_params(k, rand = TRUE))

th <- gen_params(k, rand = FALSE)
th_s <- gen_params(k, rand = TRUE)
(max(c (px_z(th) / px_z(th_s), px_z(th_s) / px_z(th))) - 1) / vec_norm(th$lambda - th_s$lambda)

cs_aI <- function(th, th_s) {
  
  Sigma <- th_s$Sigma
  k <- ncol(Sigma)
  
  all_z <- lapply(
    seq_len(2^k) - 1,
    function(i) binsensate:::cpp_idx_to_bit(i, k)
  )
  all_z <- do.call(rbind, all_z)
  
  pz <- binsensate:::cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  
  pxz_diff <- px_z(th) - px_z(th_s)
  
  # cat("Grad Diff.", (pxz_diff * pz) %*% cbind(1, all_z))
  
  vec_norm((pxz_diff * pz) %*% cbind(1, all_z))
}

th_s <- gen_params(k, rand = FALSE, eps = 0)
for (eps in c(1, 1/2, 1/4, 1/8, 1/16, 1/32)) {
  
  th <- gen_params(k, rand = FALSE, eps = eps)
  
  cat(cs_aI(th, th_s) / vec_norm(th$lambda - th_s$lambda), "\n")
  
  M_th <- M_theta(th, th_s, fi)
  cat("True Gamma",
    vec_norm(nabla_Q(th, th, th_s) - nabla_Q(th, th_s, th_s) ) / 
      vec_norm(em_params(th) - em_params(th_s)), "\n Gamma True v2",
    vec_norm(nabla_Q(M_th, th, th_s) - nabla_Q(M_th, th_s, th_s) ) / 
      vec_norm(em_params(th) - em_params(th_s)), "\n"
  )
}

k <- 2
th <- gen_params(k, FALSE, 0)
th$Sigma[,] <- 0

nabla_Q(th, th, th)

fi <- 0.2
th_s <- gen_params(k, rand = FALSE, eps = 1)
th <- gen_params(k, rand = FALSE, eps = 1+0.25)
vec_norm(p_v(th, fi, 0, 0) / rowSums(p_v(th_s, fi, 0, 0)) - p_v(th_s, fi, 0, 0) / rowSums(p_v(th_s, fi, 0, 0)))
vec_norm(p_v(th, fi, 0, 0) / rowSums(p_v(th, fi, 0, 0)) - p_v(th_s, fi, 0, 0) / rowSums(p_v(th_s, fi, 0, 0)))

# testing of a simplified strategy below
Ti_simp <- function(th, th_s) {
  
  tot <- 0
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      tot <- tot + sum(
        rowSums(
          abs(
            p_v(th, fi, x, y) / rowSums(p_v(th, fi, x, y)) - 
              p_v(th_s, fi, x, y) / rowSums(p_v(th_s, fi, x, y))
          )
        ) * rowSums(p_v(th_s, fi, x, y))
      )
    }
  }
  tot / vec_norm(th$lambda - th_s$lambda)
  # round(
  #   rowSums(
  #     abs(
  #       p_v(th, fi, x, y) / rowSums(p_v(th, fi, x, y)) -
  #         p_v(th_s, fi, x, y) / rowSums(p_v(th_s, fi, x, y))
  #     )
  #   ) / vec_norm(th$lambda - th_s$lambda), digits = 3
  # )
}

diff_nabla_simp <- function(th_p, th, th_s) {
  
  k <- ncol(th_s$Sigma)
  
  all_z <- lapply(
    seq_len(2^k) - 1,
    function(i) binsensate:::cpp_idx_to_bit(i, k)
  )
  all_z <- do.call(rbind, all_z)
  
  # get P_{\theta^*}(r, x, y)
  grad_tot <- rep(0, 2 * k + 3)
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      # get P_{\theta^*}(z, r | x, y)
      pv_s <- p_v(th_s, fi, x, y)
      
      # get P_{\theta^*}(r, x, y)
      p_rxy_s <- rowSums(pv_s)
      
      # get P_{\theta}(z, r | x, y)
      pv <- p_v(th, fi, x, y)
      # get P_{\theta}(z | r, x, y)
      pz_rxy <- abs(pv / rowSums(pv) - pv_s / rowSums(pv_s))
      
      # get x - P_{\theta'}(x_1 | z)
      term_x <- x - px_z(th_p)
      term_x[TRUE] <- 1
      
      pz_rxy_tx <- t(t(pz_rxy) * term_x)
      grad_x <- colSums((pz_rxy_tx * p_rxy_s) %*% cbind(1, all_z))
      
      # cat("Check that", pz_rxy_tx * p_rxy, "=", colSums(p_v(th, fi, x, y)) * term_x, "\n")
      
      # get y - P_{\theta'}(y_1 | z, x)
      term_y <- y - py_xz(th_p, x)
      term_y[TRUE] <- 1
      
      pz_rxy_ty <- t(t(pz_rxy) * term_y)
      grad_y <- colSums((pz_rxy_ty * p_rxy_s) %*% cbind(1, all_z, x))
      
      grad_tot <- grad_tot + c(grad_x, grad_y)
    }
  }
  
  return(grad_tot)
}

k <- 3
fi <- 0.5
for (eps in c(1, 1/2, 1/4, 1/8, 1/16, 1/32)) {
  
  th_s <- gen_params(k, rand = FALSE, eps = 0)
  th <- gen_params(k, rand = FALSE, eps = eps)
  cat(Ti_simp(th, th_s), "\n")
}
