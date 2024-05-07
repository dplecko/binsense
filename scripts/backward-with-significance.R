
load("inverse-robust.RData")

bd_solver <- function(X, Y, R, fi) {
  
  k <- ncol(R)
  enc_mat <- t(replicate(nrow(R), 2^(seq_len(k) - 1L)))
  
  pr <- pz <- pxy <- list(list(0, 0), list(0, 0))
  
  x0 <- y0 <- 0
  x1 <- y1 <- 1
  
  # expand fi if needed
  for (x in c(T, F))
    for (y in c(T, F))
      if (length(fi[[x + 1]][[y + 1]]) == 1 & k > 1) 
        fi[[x + 1]][[y + 1]] <- rep(fi[[x + 1]][[y + 1]], k)
  
  for (x in c(T, F)) {
    
    for (y in c(T, F)) {
      
      idx <- X == x & Y == y
      
      cprxy <- tabulate(
        rowSums(
          R[idx, , drop = FALSE] * enc_mat[idx, , drop = FALSE]
        ) + 1
      )
      cprxy <- c(cprxy, rep(0L, 2^k - length(cprxy))) # zero-padding
      cprxy <- cprxy / sum(cprxy)
      
      pxy[[x + 1]][[y + 1]] <- nrow(R[idx, , drop = FALSE]) / nrow(R)
      pr[[x + 1]][[y + 1]] <- cprxy
      
      pz[[x + 1]][[y + 1]] <- 
        vapply(
          seq_len(2^k),
          function(z) binsensate:::infer_pz(z, pr[[x + 1]][[y + 1]], 
                                            fi[[x + 1]][[y + 1]], k),
          numeric(1L)
        )
      
      if (binsensate:::check_simplex(pz[[x + 1]][[y + 1]])) {
        
        initial <- pz[[x + 1]][[y + 1]]
        A_inv_xy <- lapply(
          seq_len(2^k),
          function(z) binsensate:::infer_pz(z, pr[[x + 1]][[y + 1]], 
                                            fi[[x + 1]][[y + 1]], k, 
                               TRUE)
        )
        A_inv_xy <- do.call(rbind, A_inv_xy)
        
        # check if ill-posedness is statistically significant
        if (ill_pose_sig(pz[[x + 1]][[y + 1]], pr[[x + 1]][[y + 1]], 
                         A_inv_xy, nsamp = sum(idx))) {
          
          # trade for a well-posed inverse problem
          fi_new <- binsensate:::trade_inv_prob(pz[[x + 1]][[y + 1]], 
                                                pr[[x + 1]][[y + 1]], 
                                                A_inv_xy, fi[[x + 1]][[y + 1]])
          
          # re-run with new fi
          fi_new <- binsensate:::fi_trade(fi, fi_new, x, y)
          return(
            bd_solver(X, Y, R, fi_new)
          )
        }
      }
      
      if (min(pz[[x + 1]][[y + 1]]) < 0 | 
          max(pz[[x + 1]][[y + 1]]) > 1) {
        
        neg_idx <- pz[[x + 1]][[y + 1]] < 0
        pz[[x + 1]][[y + 1]][neg_idx] <- 0
        pz[[x + 1]][[y + 1]] <- pz[[x + 1]][[y + 1]] / sum(pz[[x + 1]][[y + 1]])
      }
    }
  }
  
  pz_inf <- rep(0, 2^k)
  for (x in c(T, F)) {
    for (y in c(T, F)) {
      pz_inf <- pz_inf + pxy[[x + 1]][[y + 1]] * pz[[x + 1]][[y + 1]]
    }
  }
  
  py_x1z <- pz[[x1 + 1]][[y1 + 1]] * pxy[[x1 + 1]][[y1 + 1]] / (
    pz[[x1 + 1]][[y1 + 1]] * pxy[[x1 + 1]][[y1 + 1]] + 
      pz[[x1 + 1]][[y0 + 1]] * pxy[[x1 + 1]][[y0 + 1]]
  )
  
  py_x0z <- pz[[x0 + 1]][[y1 + 1]] * pxy[[x0 + 1]][[y1 + 1]] / (
    pz[[x0 + 1]][[y1 + 1]] * pxy[[x0 + 1]][[y1 + 1]] + 
      pz[[x0 + 1]][[y0 + 1]] * pxy[[x0 + 1]][[y0 + 1]]
  )
  
  py_x1z[is.nan(py_x1z)] <- 0
  py_x0z[is.nan(py_x0z)] <- 0
  ate <- sum((py_x1z - py_x0z) * pz_inf)
  
  return(
    list(pz = pz_inf, ATE = ate, fi = fi, solver = "backward")
  )
}

backward_solver(Xsim[[4]], Ysim[[4]], Rsim[[4]], 
                fi = list(list(0.15, 0.15), list(0, 0)))

th <- list(
  lambda = lambda[-1], l_icept = lambda[1],
  mu = mu[-1], m_icept = mu[1],
  Sigma = Sigma_mim0, beta = dat_sim$beta
)

p_z_xy <- function(th, fi, x, y) {
  
  Sigma <- th$Sigma
  k <- ncol(Sigma)
  
  pz <- binsensate:::cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  
  pxz <- px_z(th)
  if (x == 0) pxz <- 1 - pxz
  
  pyxz <- py_xz(th, x)
  if (y == 0) pyxz <- 1 - pyxz
  
  pz * pxz * pyxz
}
