
explr_traj <- function(k, at_max = FALSE) {
  
  cvx_cmb <- function(th1, th2, scl) {
    
    Map(function(x, y) scl * x + (1-scl) * y, th1, th2)
  }
  
  res <- NULL
  
  th_s <- gen_params(k, rand = TRUE)
  th <- gen_params(k, rand = TRUE)
  th_p <- gen_params(k, rand = TRUE)
  cov_evals <- evals_Sigma(th_s$Sigma)
  
  
  
  scl_seq <- seq(0, 0.9, 0.1)
  for (scl in scl_seq) {
    
    th_cvx <- cvx_cmb(th_s, th, scl)
    dist <- vec_norm(em_params(th_cvx) - em_params(th_s))
    
    if (at_max) {
      
      M_th <- M_theta(th_cvx, th_s, fi)
      th_p <- M_th
    }
    
    
    nab_c <- nabla_Q(th_p, th_cvx, th_s, full = TRUE)
    nab_s <- nabla_Q(th_p, th_s, th_s, full = TRUE)
    gamma <- vec_norm(nab_c$grad_tot - nab_s$grad_tot) / dist
    for (x in c(0, 1)) {
      
      for (y in c(0, 1)) {
        
        diff_mat <- cbind(nab_c$grad_termsx[[1+x]][[1+y]], 
                          nab_c$grad_termsy[[1+x]][[1+y]]) -
          cbind(nab_s$grad_termsx[[1+x]][[1+y]], 
                nab_s$grad_termsy[[1+x]][[1+y]])
        gamma_byt <- sqrt(rowSums(diff_mat^2)) / dist
        
        res <- rbind(
          res,
          cbind(gam_byt = gamma_byt, x = x, y = y, 
                z = seq_along(gamma_byt), gam = gamma,
                scl = scl)
        )
      }
    }
  }
  
  ggplot(as.data.frame(res), aes(x = z, y = gam_byt, color = interaction(x, y))) +
    geom_point() + geom_line() +
    geom_hline(aes(yintercept = gam), color = "red", linetype = "dashed") +
    geom_hline(yintercept = min(cov_evals), color = "orange") +
    facet_wrap(~ scl) + theme_bw()
}

set.seed(2025)
explr_traj(3, at_max = TRUE)
