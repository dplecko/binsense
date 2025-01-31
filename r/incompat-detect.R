

likelihood_v <- function(X, Y, R, fi, A, lmb, Sigma) {
  
  # decide whether to do by chunk
  k <- ncol(R)
  if (k > 5) { # possibly consider a chunking approach
    
    # re-order
    enc <- 1 + rowSums(R * t(replicate(nrow(dat$R), 2^(seq_len(k) - 1))))
    reord <- order(enc)
    R <- R[reord, ]
    X <- X[reord]
    Y <- Y[reord]
    enc <- enc[reord]
    
    chunks <- rle(enc)
    enc_vals <- chunks$values
    brks <- c(1, cumsum(chunks$lengths))
    starts <- brks[-length(brks)]
    stops <- brks[-1]
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

cross_entropy <- function(p, y) -mean(log(p) * y + log(1-p) * (1-y))

cv_xgboost <- function(df, y, ret = "preds") {
  
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(df), label = y)
  params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss"
  )
  
  cv <- xgboost::xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 1000,
    nfold = 5,
    early_stopping_rounds = 10,
    prediction = TRUE,
    verbose = FALSE
  )
  
  if (ret == "nrounds") return(cv$best_iteration)
  
  return(cv$pred)
}

incompat_detect_like <- function(X, Y, R, fi, A, lmb, Sigma) {
  
  like0 <- likelihood_v(X, Y, R, fi, A, lmb, Sigma)
  
  nrep <- 100
  for (i in seq_len(nrep)) {
    
    # adversary based
    dat_geni <- synth_data_mid(
      n = nrow(R), k = nrow(Sigma),
      A = A, class = "expfam-2d", seed = i,
      Sigma = Sigma, lam = lmb$lambda, mu = lmb$mu,
      icept_x = lmb$lambda_icept, icept_y = lmb$mu_icept,
      beta = lmb$beta
    )
    
    like_i <- likelihood_v(dat_geni$X, dat_geni$Y, dat_geni$R, fi, A, lmb, Sigma)
    if (like_i < like0) return(FALSE)
  }
  
  return(TRUE)
}

incompat_detect_class <- function(X, Y, R, fi, A, lmb, Sigma) {
  
  # adversary based
  dat_gen <- synth_data_mid(
    n = nrow(R), k = nrow(Sigma),
    A = A, class = "expfam-2d",
    Sigma = Sigma, lam = lmb$lambda, mu = lmb$mu,
    icept_x = lmb$lambda_icept, icept_y = lmb$mu_icept,
    beta = lmb$beta
  )
  
  dt_tru <- cbind(X = X, Y = Y, R = R)
  dt_gen <- cbind(X = dat_gen$X, Y = dat_gen$Y, R = dat_gen$R)
  
  dt <- rbind(dt_tru, dt_gen)
  labs <- rep(c(1, 0), each = nrow(R))
  
  preds <- cv_xgboost(dt, y = labs)
  loss0 <- cross_entropy(preds, labs)
  
  nreps <- 100
  loss_perm <- c()
  for (i in seq_len(nreps)) {
    
    labs_perm <- sample(labs, length(labs))
    preds_perm <- cv_xgboost(dt, y = labs_perm)
    loss_perm <- c(loss_perm, cross_entropy(preds_perm, labs_perm))
  }
  
  # compute the p-value based on permutation testing
  mean( loss_perm < loss0 )
}

incompat_detect <- function(X, Y, R, fi, A, lmb, Sigma) {
  
  tl <- length(lmb)
  lmb <- lmb[[tl]]
  Sigma <- Sigma[[tl]]
  
  incompat_detect_like(X, Y, R, fi, A, lmb, Sigma)
}

zinf_verify <- function(X, Y, R, fi) {
  
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
  
  max_fi <- zinf_max_fi(X, Y, R)
  for (x in c(0, 1)) for (y in c(0, 1)) {
    
    if (any(fi[[1+x]][[1+y]] > max_fi[[1+x]][[1+y]])) {
      
      message(
        paste0(
          "Specified fi value would imply that ",
          "P(R = (0,...,0) | X = ", x, ", Y = ", y, ") < 0. \n",
          "Consider adjusting the parameter value(s) of the ZINF method.\n",
          "For two-stage-em solver, you may wish to run with argument ",
          "`detect_viol = TRUE`"
        )
      )
    }
  }
}


