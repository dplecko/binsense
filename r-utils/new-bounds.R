
# attempt of getting the Cauchy-Schwarz bound
csb_bound <- function(th_p, th, th_s, full = FALSE) {
  
  k <- ncol(th_s$Sigma)
  
  all_z <- lapply(
    seq_len(2^k) - 1,
    function(i) binsensate:::cpp_idx_to_bit(i, k)
  )
  all_z <- do.call(rbind, all_z)
  
  # get P_{\theta^*}(r, x, y)
  grad_tot <- rep(0, k + 1)
  grad_termsx <- grad_termsy <- list(list(NULL, NULL), list(NULL, NULL))
  tot_termI <- tot_termII <- 0
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      # get P_{\theta^*}(r, x, y)
      pv_s <- p_v(th_s, fi, x, y)
      p_rxy_s <- rowSums(pv_s)
      
      # get P_{\theta^*}(z | r, x, y)
      pz_rxy_s <- pv_s / rowSums(pv_s)
      
      # get P_{\theta}(z | r, x, y)
      pv <- p_v(th, fi, x, y)
      pz_rxy <- pv / rowSums(pv)
      
      a_div_b <- (1 - pz_rxy / pz_rxy_s)^2
      a_div_b[is.nan(a_div_b)] <- 0
      
      csb_termI <- sum(pv_s * a_div_b)
      tot_termI <- tot_termI + csb_termI
      
      term_x <- colSums(pv_s) * (x - px_z(th_p))^2
      
      csb_termII <- term_x %*% cbind(1, all_z)
      tot_termII <- tot_termII + sum(csb_termII)
    }
  }
  
  sqrt(tot_termI * tot_termII)
}

# computing the second order derivative of Q(\theta^*, \theta^*)
Qhess <- function(th_s) {
  
  th_p <- th_s
  k <- ncol(th_s$Sigma)
  
  all_z <- lapply(
    seq_len(2^k) - 1,
    function(i) binsensate:::cpp_idx_to_bit(i, k)
  )
  all_z <- do.call(rbind, all_z)
  
  
  # get P(z), P(x1 | z)
  pz_s <- binsensate:::cpp_scaling_2d(th_s$Sigma)
  pz_s <- pz_s / sum(pz_s)
  px1_z <- px_z(th_p)
  
  W_z <- diag( pz_s * px1_z * (1 - px1_z) )
  Qhess_x <- t(cbind(1, all_z)) %*% W_z %*% cbind(1, all_z)
  
  # for y this is slightly more complicated since x cannot be integrated out (!)
  Qhess_y <- array(0, dim = c(k+2, k+2))
  for (x in c(0, 1)) {
    
    p_x_z <- px_z(th_p)
    if (x == 0) p_x_z <- 1 - p_x_z
    pxz <- pz_s * p_x_z
    py1_xz <- py_xz(th_p, x)
    W_xz <- diag( pxz * py1_xz * (1 - py1_xz) )
    Qhess_y <- Qhess_y + t(cbind(1, all_z, x)) %*% W_xz %*% cbind(1, all_z, x)
  }
  
  Matrix::bdiag(Qhess_x, Qhess_y)
}
