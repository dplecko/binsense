
binsensate <- function(X, Y, R, fi, 
                       solver = c("backward", "expfam-1d", "expfam-2d"),
                       mc_samples = 10^4, ...) {
  
  solver <- match.arg(solver, c("backward", "expfam-1d", "expfam-2d"))
  
  switch(
    solver,
    backward = backward_solver(X, Y, R, fi, ...),
    `expfam-1d` = expfam_1d_solver(X, Y, R, fi, mc_samples),
    `expfam-2d` = expfam_2d_solver(X, Y, R, fi, mc_samples)
  )
}

expfam_1d_solver <- function(X, Y, R, fi, mc_samples) {
  
  lambda <- em <- list(list(NULL, NULL), list(NULL, NULL))
  Z_mc <- X_mc <- Y_mc <- weights <- NULL
  
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      idx <- X == x & Y == y
      
      R[idx, , drop = FALSE]
      
      em[[1+x]][[1+y]] <- descent_1d(rep(0, ncol(R)), R[idx, , drop = FALSE], 
                       fi = fi[[1+x]][[1+y]],
                       gamma = 10, n_iter = 100)
      
      lambda[[1+x]][[1+y]] <- tail(em[[1+x]][[1+y]]$lambda, n = 1L)[[1]] 
                            # aux$lxy[[1+x]][[1+y]]
      
      # get Monte Carlo samples of Z
      Z_mc <- rbind(Z_mc, z_1d(lambda = lambda[[1+x]][[1+y]], 
                               n_samp = mc_samples))
      X_mc <- c(X_mc, rep(x, mc_samples))
      Y_mc <- c(Y_mc, rep(y, mc_samples))
      weights <- c(weights, rep(sum(idx) / nrow(R), mc_samples))
    }
  }

  logreg_ATE(X_mc, Y_mc, Z_mc, weights)
}

expfam_2d_solver <- function(X, Y, R, fi, mc_samples) {
  
  Sigma <- list(list(NULL, NULL), list(NULL, NULL))
  Z_mc <- X_mc <- Y_mc <- weights <- NULL
  
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      idx <- X == x & Y == y
      
      R[idx, , drop = FALSE]
      
      em <- descent_2d(array(0, dim = c(ncol(R), ncol(R))), 
                       R[idx, , drop = FALSE], 
                       fi = fi[[1+x]][[1+y]],
                       gamma = 5, n_iter = 5)
      
      Sigma[[1+x]][[1+y]] <- tail(em$Sigma, n = 1L)[[1]]
      
      # get Monte Carlo samples of Z
      Z_mc <- rbind(Z_mc, z_2d(Sigma = Sigma[[1+x]][[1+y]], 
                               n_samp = mc_samples))
      X_mc <- c(X_mc, rep(x, mc_samples))
      Y_mc <- c(Y_mc, rep(y, mc_samples))
      weights <- c(weights, rep(sum(idx) / nrow(R), mc_samples))
    }
  }
  
  logreg_ATE(X_mc, Y_mc, Z_mc, weights, "expfam-2d")
}

backward_solver <- function(X, Y, R, fi, verbose = FALSE, ...) {
  
  k <- ncol(R)
  enc_mat <- t(replicate(nrow(R), 2^(seq_len(k) - 1L)))
  
  # get P(y | r, x)
  logreg <- glm(Y ~ ., data = data.frame(Y, X, R), family = "binomial")
  pyrx0 <- pyrx1 <- rep(0, 2^k)
  for (i in seq_len(2^k)) {
    
    i_bin <- as.integer(intToBits(i-1))[1:k]
    pyrx0[i] <- expit(sum(logreg$coefficients * c(1, 0, i_bin)))
    pyrx1[i] <- expit(sum(logreg$coefficients * c(1, 1, i_bin)))
    
  }
  
  pr <- pz <- pxy <- A <- A_inv <- B <- list(list(0, 0), list(0, 0))
  
  x0 <- y0 <- 0
  x1 <- y1 <- 1
  
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
      
      A[[x + 1]][[y + 1]] <- cpp_A_xy(fi_xy = fi[[x + 1]][[y + 1]], k = k)
      
      if (any(A[[x + 1]][[y + 1]] < 0) | any(A[[x + 1]][[y + 1]] < 0)) 
        return(list(ATE = NA, cmb_diff = NA))
      
      
      pz[[x + 1]][[y + 1]] <- 
        vapply(
          seq_len(2^k),
          function(z) infer_pz(z, pr[[x + 1]][[y + 1]], fi[[x + 1]][[y + 1]], k),
          numeric(1L)
        )
      
      # assert_that(sum(abs(pz[[x + 1]][[y + 1]] - pz_cand)) < 10^(-10))
      
      A_inv[[x+1]][[y+1]] <- lapply(
        seq_len(2^k),
        function(z) infer_pz(z, pr[[x + 1]][[y + 1]], fi[[x + 1]][[y + 1]], k, T)
      )
      A_inv[[x+1]][[y+1]] <- do.call(rbind, A_inv[[x+1]][[y+1]])
      
      # B represents P(z | r, x, y)
      B[[x + 1]][[y + 1]] <- cpp_invert_A_xy(A[[x + 1]][[y + 1]], 
                                             pz[[x + 1]][[y + 1]], 
                                             pr[[x + 1]][[y + 1]])
      
      
      if (check_simplex(pz[[x + 1]][[y + 1]])) {
        initial <- pz[[x + 1]][[y + 1]]
        pz[[x + 1]][[y + 1]] <- locate_simplex(pr[[x + 1]][[y + 1]], 
                                               A_inv[[x + 1]][[y + 1]])
        correction <- sum(abs(initial - pz[[x+1]][[y+1]]) * 100)
        if (verbose) {
          cat("Locate simplex correction ", round(correction, 2), "%\n", sep = "")
        }
      }
      
      if (min(pz[[x + 1]][[y + 1]]) < -(10)^(-10) | max(pz[[x + 1]][[y + 1]]) > 1) {
        cat("Walking out of simplex for (x,y)", x, y, "and range", 
            range(pz[[x + 1]][[y + 1]]), "\n")
      }
      
    }
    
  }
  
  pz_inf <- rep(0, 2^k)
  for (x in c(T, F)) {
    for (y in c(T, F)) {
      pz_inf <- pz_inf + pxy[[x + 1]][[y + 1]] * pz[[x + 1]][[y + 1]]
    }
  }
  
  dlt_x0x1 <- rep(0, 2^k)
  for (i in seq_len(2^k)) {
    
    
    ratio_x0 <- B[[x0 + 1]][[y1 + 1]][i, ] / 
      (B[[x0 + 1]][[y1 + 1]][i, ] * pyrx0 + B[[x0+1]][[y0+1]][i, ] * (1-pyrx0))
    ratio_x0[is.nan(ratio_x0)] <- 0
    
    ratio_x1 <- B[[x1+1]][[y1+1]][i, ] / 
      (B[[x1+1]][[y1+1]][i, ] * pyrx1 + B[[x1+1]][[y0+1]][i, ] * (1-pyrx1))
    ratio_x1[is.nan(ratio_x1)] <- 0
    
    # can the fi_xy case not be solved?! additional solving needed?
    przx0 <- A[[x0 + 1]][[y1 + 1]][, i] 
    przx1 <- A[[x1 + 1]][[y1 + 1]][, i]
    
    dlt_x0x1[i] <- sum(pyrx1 * ratio_x1 * przx1 - 
                         pyrx0 * ratio_x0 * przx0) # \sum_r {...}
  }
  
  wgh_vec <- sapply(seq_along(pz_inf) - 1, 
                    function(x) sum(as.integer(intToBits(x))))
  
  list(
    pz = pz_inf,
    ATE = sum(dlt_x0x1 * pz_inf),
    solver = "backward",
    delta = dlt_x0x1
  )
}

logreg_ATE <- function(X, Y, Z, weights, solver = "expfam-1d") {
  
  dat <- data.frame(X, Y, Z)
  
  # bootstrap the thingy
  dat <- dat[sample.int(nrow(dat), replace = TRUE, prob = weights), ]
  mod <- glm(
    Y ~ ., data = dat, family = "binomial"
  )

  dat$X <- 1
  y1 <- predict(mod, newdata = dat, type = "response")
  
  dat$X <- 0
  y0 <- predict(mod, newdata = dat, type = "response")
  
  Z_w <- Z[sample.int(nrow(Z), nrow(Z), TRUE, prob = weights), , drop = FALSE]
  pz_marg <- tabulate(
      rowSums(
        Z_w * t(replicate(nrow(Z_w), 2^(seq_len(ncol(Z)) - 1L)))
      ) + 1
    )
  pz_marg <- c(pz_marg, rep(0L, 2^ncol(Z) - length(pz_marg))) # zero-padding
  pz_marg <- pz_marg / sum(pz_marg)
  
  list(
    pz = pz_marg,
    ATE = mean(y1) - mean(y0),
    solver = solver
  )
}
