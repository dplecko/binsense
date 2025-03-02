

likelihood_v <- function(X, Y, R, fi, A, lmb, Sigma) {
  
  # decide whether to do by chunk
  k <- ncol(R)
  if (k > 5) { # possibly consider a chunking approach
    
    # re-order
    enc <- 1 + rowSums(R * t(replicate(nrow(R), 2^(seq_len(k) - 1))))
    reord <- order(enc)
    R <- R[reord, ]
    X <- X[reord]
    Y <- Y[reord]
    enc <- enc[reord]
    
    chunks <- rle(enc)
    enc_vals <- chunks$values
    starts <- c(1, 1 + cumsum(chunks$lengths))
    starts <- starts[-length(starts)]
    stops <- cumsum(chunks$lengths)
    zero_idx <- which(enc == 1) # seq_len(stops[1])
    nzero_idx <- setdiff(seq_len(length(X)), zero_idx)
    
    by_chunk <- if ((length(enc_vals) / length(X)) < 0.2) TRUE else FALSE
  } else by_chunk <- FALSE
  
  # get px_z, py_z
  pz <- cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  px_z <- expit(cpp_p_z(lmb$lambda) + lmb$lambda_icept)
  logity_z <- cpp_p_z(lmb$mu) + lmb$mu_icept
  beta <- lmb$beta
  
  if (by_chunk) {
    
    like_i <- Map(
      function(r_idx, start, stop) {
        
        X_chunk <- X[start:stop]
        Y_chunk <- Y[start:stop]
        ret <- 0
        for (x in c(0, 1)) for (y in c(0, 1)) {
          
          num_xy <- sum(X_chunk == x & Y_chunk == y)
          if (num_xy == 0) next
          
          # P(r | z, x, y) conditionals
          pr_zxy <- A[[1 + x]][[1 + y]][r_idx, ]
          
          if (x == 1) px <- px_z else px <- 1 - px_z
          
          py_z <- expit(logity_z + beta * x)
          if (y == 1) py <- py_z else py <- 1 - py_z
          
          ret <- ret + num_xy * log(sum(pz * px * py * pr_zxy))
        }

        ret
      },
      enc_vals, starts, stops
    )
    like_i <- do.call(c, like_i)
  } else {
    
    like_i <- vapply(
      seq_len(nrow(R)),
      function(j) {
        
        x <- X[j]
        y <- Y[j]
        r <- R[j, ]
        
        r_idx <- 1 + sum(r * 2^(seq_len(k)-1))
        
        # P(r | z, x, y) conditionals
        pr_zxy <- A[[1 + x]][[1 + y]][r_idx, ]
        
        if (x == 1) px <- px_z else px <- 1 - px_z
        
        py_z <- expit(logity_z + beta * x)
        if (y == 1) py <- py_z else py <- 1 - py_z
        
        pv <- pz * px * py * pr_zxy
        log(sum(pv))
      }, numeric(1L)
    )
  }
  
  sum(like_i)
}

incompat_detect_like <- function(X, Y, R, fi, A, lmb, Sigma) {
  
  like0 <- likelihood_v(X, Y, R, fi, A, lmb, Sigma)
  
  nrep <- 100
  for (i in seq_len(nrep)) {
    
    dat_geni <- param_mod_gen(
      n = nrow(R), k = nrow(Sigma),
      seed = i, A = A, 
      Sigma = Sigma, lam = lmb$lambda, mu = lmb$mu,
      icept_x = lmb$lambda_icept, icept_y = lmb$mu_icept,
      beta = lmb$beta
    )
    
    like_i <- likelihood_v(dat_geni$X, dat_geni$Y, dat_geni$R, fi, A, lmb, Sigma)
    if (like_i < like0) return(FALSE)
  }
  
  return(TRUE)
}

incompat_detect <- function(X, Y, R, fi, A, lmb, Sigma) {
  
  tl <- length(lmb)
  lmb <- lmb[[tl]]
  Sigma <- Sigma[[tl]]
  
  incompat_detect_like(X, Y, R, fi, A, lmb, Sigma)
}

zinf_max_fi <- function(X, Y, R) {
  
  n <- length(X)
  k <- ncol(R)
  encR <- colSums(t(R) * 2^(seq_len(k)-1))
  
  max_fi <- list(list(1, 1), list(1, 1))
  for (x in c(0, 1)) for (y in c(0, 1)) {
    
    idx <- (X == x) & (Y == y)
    max_fi[[1+x]][[1+y]] <- mean(encR[idx] == 0)
  }
  
  max_fi
}

zinf_verify <- function(X, Y, R, fi) {
  
  max_fi <- zinf_max_fi(X, Y, R)
  pairs <- c()
  for (x in c(0, 1)) for (y in c(0, 1)) {
    
    if (any(fi[[1+x]][[1+y]] > max_fi[[1+x]][[1+y]])) {
      
      pairs <- rbind(pairs, c(x, y))
    }
  }
  
  if (!is.null(pairs)) {
    
    xy_col <- c()
    for (i in seq_len(nrow(pairs))) {
      
      xy_col <- c(xy_col, paste0("(", pairs[i, 1], ",", pairs[i, 2], ")"))
    }
    xy_col <- paste0(xy_col, collapse = ", ")
    
    message(
      paste0(
        "fi value would imply that ",
        "P(R = (0,...,0) | X = x, Y = y) < 0 for (x,y): ", xy_col, ". \n",
        "Consider adjusting fi value(s) of the ZINF method.\n",
        "For the em solver, you may wish to run with argument ",
        "`detect_viol = TRUE`"
      )
    )
  }
}
