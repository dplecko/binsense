
binsensate <- function(X, Y, R, fi, 
                       solver = c("backward", "expfam-1d", "expfam-2d",
                                  "expfam-1d-direct", "backward-direct",
                                  "expfam-2d-mom", "expfam-2d-mom-grad",
                                  "two-stage-em"),
                       mc_samples = 10^4, mc_per_samp = 10, 
                       ret_em = FALSE, ...) {
  
  solver <- match.arg(solver, c("backward", "expfam-1d", "expfam-2d", 
                                "expfam-1d-direct", "backward-direct",
                                "expfam-2d-mom", "expfam-2d-mom-grad",
                                "two-stage-em"))
  
  switch(
    solver,
    backward = backward_solver(X, Y, R, fi, ...),
    `expfam-1d` = expfam_1d_solver(X, Y, R, fi, mc_samples),
    `expfam-2d` = expfam_2d_solver(X, Y, R, fi, mc_samples, ret_em = ret_em),
    `expfam-1d-direct` = expfam_1d_direct(X, Y, R, fi, mc_samples),
    `backward-direct` = backward_solver(X, Y, R, fi, direct = TRUE, ...),
    `expfam-2d-mom` = expfam_2d_mom(X, Y, R, fi, mc_samples),
    `expfam-2d-mom-grad` = expfam_2d_mom_grad(X, Y, R, fi, mc_samples),
    `two-stage-em` = two_stage_em(X, Y, R, fi, mc_per_samp)
  )
}

expfam_1d_solver <- function(X, Y, R, fi, mc_samples) {
  
  lambda <- em <- list(list(NULL, NULL), list(NULL, NULL))
  Z_mc <- X_mc <- Y_mc <- weights <- NULL
  
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      idx <- X == x & Y == y
      em[[1+x]][[1+y]] <- descent_1d(rep(0, ncol(R)), R[idx, , drop = FALSE], 
                                     fi = fi[[1+x]][[1+y]],
                                     gamma = 10, n_iter = 10)
      
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

expfam_1d_direct <- function(X, Y, R, fi, mc_samples) {
  
  lambda <- list(list(NULL, NULL), list(NULL, NULL))
  Z_mc <- X_mc <- Y_mc <- weights <- NULL
  
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      idx <- X == x & Y == y
      lambda[[1+x]][[1+y]] <- logit(
        colMeans(R[idx, , drop = FALSE]) / (1 -  fi[[1+x]][[1+y]])
      )
      
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

expfam_2d_solver <- function(X, Y, R, fi, mc_samples, ret_em = FALSE) {
  
  Sigma <- em_obj <- list(list(NULL, NULL), list(NULL, NULL))
  Z_mc <- X_mc <- Y_mc <- weights <- NULL
  
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      idx <- X == x & Y == y
      
      Rxy <- R[idx, , drop = FALSE]
      fixy <- fi[[1+x]][[1+y]]
      
      erer <- colMeans(Rxy) %*% t(colMeans(Rxy))
      err <- var(Rxy) - erer
      
      fi_mat <- array((1-fixy)^2, dim = dim(err))
      diag(fi_mat) <- 1-fixy 
      varz <- err / fi_mat - erer / (1-fixy)^2
      
      em <- descent_2d(varz, Rxy, fi = fixy,
                       gamma = 5, n_iter = 50,
                       true_Sigma = attr(R, "Sigma"),
                       # true_Sigma = k2_solver(Rxy, fi[[1+x]][[1+y]]),
                       type = "Newton-Raphson")
      
      if (ret_em) em_obj[[1+x]][[1+y]] <- em 
      # browser()
      # vis_descent_2d(em, 30:40)
      
      Sigma[[1+x]][[1+y]] <- tail(em$Sigma, n = 1L)[[1]]
      
      # get Monte Carlo samples of Z
      Z_mc <- rbind(Z_mc, z_2d(Sigma = Sigma[[1+x]][[1+y]], 
                               n_samp = mc_samples))
      X_mc <- c(X_mc, rep(x, mc_samples))
      Y_mc <- c(Y_mc, rep(y, mc_samples))
      weights <- c(weights, rep(sum(idx) / nrow(R), mc_samples))
    }
  }
  
  if (ret_em) return(em_obj)
  logreg_ATE(X_mc, Y_mc, Z_mc, weights, "expfam-2d")
}

expfam_2d_mom <- function(X, Y, R, fi, mc_samples, parallel = TRUE) {
  
  Sigma <- list(list(NULL, NULL), list(NULL, NULL))
  Z_mc <- X_mc <- Y_mc <- weights <- NULL
  k <- ncol(R)
  
  if (parallel) {
    
    mcs <- mclapply(
      0:3,
      function(case) {
        x <- as.integer(case / 2)
        y <- as.integer(case %% 2)
        
        idx <- X == x & Y == y
        Rxy <- R[idx, , drop = FALSE]
        fixy <- fi[[1+x]][[1+y]]
        
        fi_mat <- array((1-fixy)^2, dim = c(k, k))
        diag(fi_mat) <- 1-fixy
        
        mu_hat <- t(Rxy) %*% Rxy / nrow(Rxy) / fi_mat
        
        if (max(mu_hat) > 1) {
          
          message("Data likely not generated from this \\phi value...\n")
        }
        
        Sigma_zero <- array(0, dim = c(k, k))
        Sigma[[1+x]][[1+y]] <- mom_gradient(Sigma_zero, mu_hat, -1, nrow(Rxy))
        
        list(
          Z = z_2d(Sigma = Sigma[[1+x]][[1+y]],
                   n_samp = mc_samples),
          X = rep(x, mc_samples),
          Y = rep(y, mc_samples),
          weights = rep(sum(idx) / nrow(R), mc_samples)
        )
      }, mc.cores = 4
    )

    Z_mc <- do.call(rbind, lapply(mcs, `[[`, "Z"))
    X_mc <- do.call(c, lapply(mcs, `[[`, "X"))
    Y_mc <- do.call(c, lapply(mcs, `[[`, "Y"))
    weights <- do.call(c, lapply(mcs, `[[`, "weights"))
  } else {
    
    for (x in c(0, 1)) {

      for (y in c(0, 1)) {

        idx <- X == x & Y == y
        Rxy <- R[idx, , drop = FALSE]
        fixy <- fi[[1+x]][[1+y]]

        fi_mat <- array((1-fixy)^2, dim = c(k, k))
        diag(fi_mat) <- 1-fixy

        mu_hat <- t(Rxy) %*% Rxy / nrow(Rxy) / fi_mat

        if (max(mu_hat) > 1) {

          message("Data likely not generated from this \\phi value...\n")
        }

        Sigma_zero <- array(0, dim = c(k, k))
        Sigma[[1+x]][[1+y]] <- mom_gradient(Sigma_zero, mu_hat, -1, nrow(Rxy))

        # get Monte Carlo samples of Z
        Z_mc <- rbind(Z_mc, z_2d(Sigma = Sigma[[1+x]][[1+y]],
                                 n_samp = mc_samples))
        X_mc <- c(X_mc, rep(x, mc_samples))
        Y_mc <- c(Y_mc, rep(y, mc_samples))
        weights <- c(weights, rep(sum(idx) / nrow(R), mc_samples))
      }
    }
  }

  logreg_ATE(X_mc, Y_mc, Z_mc, weights, "expfam-2d")
}

expfam_2d_mom_grad <- function(X, Y, R, fi, mc_samples, parallel = FALSE) {
  
  Sigma <- pxy <- list(list(NULL, NULL), list(NULL, NULL))
  Z_mc <- X_mc <- Y_mc <- weights <- NULL
  k <- ncol(R)
  
  if (parallel) {
    
    mcs <- mclapply(
      0:3,
      function(case) {
        x <- as.integer(case / 2)
        y <- as.integer(case %% 2)
        
        idx <- X == x & Y == y
        Rxy <- R[idx, , drop = FALSE]
        fixy <- fi[[1+x]][[1+y]]
        
        Sigma_zero <- array(0, dim = c(k, k))
        Sigma_opt <- r_descent_2d(Sigma_zero, Rxy, fixy, multistep = TRUE)
        Sigma_opt <- tail(Sigma_opt, n = 1L)[[1]]
        Sigma[[1+x]][[1+y]] <- Sigma_opt
        pxy[[1+x]][[1+y]] <- mean(idx)
        
        list(
          Z = z_2d(Sigma = Sigma[[1+x]][[1+y]],
                   n_samp = mc_samples),
          X = rep(x, mc_samples),
          Y = rep(y, mc_samples),
          weights = rep(sum(idx) / nrow(R), mc_samples)
        )
      }, mc.cores = 4
    )
    Z_mc <- do.call(rbind, lapply(mcs, `[[`, "Z"))
    X_mc <- do.call(c, lapply(mcs, `[[`, "X"))
    Y_mc <- do.call(c, lapply(mcs, `[[`, "Y"))
    weights <- do.call(c, lapply(mcs, `[[`, "weights"))
  } else {
    
    for (x in c(0, 1)) {
      
      for (y in c(0, 1)) {
        
        idx <- X == x & Y == y
        Rxy <- R[idx, , drop = FALSE]
        fixy <- fi[[1+x]][[1+y]]
        
        Sigma_zero <- array(0, dim = c(k, k))
        Sigma_opt <- r_descent_2d(Sigma_zero, Rxy, fixy, multistep = TRUE)
        Sigma_opt <- tail(Sigma_opt, n = 1L)[[1]]
        Sigma[[1+x]][[1+y]] <- Sigma_opt
        pxy[[1+x]][[1+y]] <- mean(idx)
        
        # get Monte Carlo samples of Z
        Z_mc <- rbind(Z_mc, z_2d(Sigma = Sigma[[1+x]][[1+y]],
                                 n_samp = mc_samples))
        X_mc <- c(X_mc, rep(x, mc_samples))
        Y_mc <- c(Y_mc, rep(y, mc_samples))
        weights <- c(weights, rep(sum(idx) / nrow(R), mc_samples))
      }
    }
  }
  
  ateor_Sigma_xy(Sigma, pxy)
  # logreg_ATE_flex(X_mc, Y_mc, Z_mc, weights, "expfam-2d", Sigma)
}

ateor_Sigma_xy <- function(Sigma, pxy) {
  
  pz <- rep(0, 2^ncol(Sigma[[1]][[1]]))
  
  py_x <- pz_xy <- pxy # inherits the 2 x 2 list structure
  
  x1 <- y1 <- 1
  x0 <- y0 <- 0
  
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      pz_xy[[1+x]][[1+y]] <- cpp_scaling_2d(Sigma[[1+x]][[1+y]])
      pz_xy[[1+x]][[1+y]] <- pz_xy[[1+x]][[1+y]] / sum(pz_xy[[1+x]][[1+y]])
      pz <- pz + pz_xy[[1+x]][[1+y]] * pxy[[1+x]][[1+y]]
      
      py_x[[1+x]][[1+y]] <- pxy[[1+x]][[1+y]] / (pxy[[1+x]][[1]] + pxy[[1+x]][[2]])
    }
  }
  
  py1_x1z <- py_x[[2]][[2]] * pz_xy[[2]][[2]] / 
    (py_x[[1+x1]][[1]] * pz_xy[[1+x1]][[1]] + py_x[[1+x1]][[2]] * pz_xy[[2]][[2]]) 
  
  py1_x0z <- py_x[[1]][[2]] * pz_xy[[1]][[2]] / 
    (py_x[[1+x0]][[1]] * pz_xy[[1+x0]][[1]] + py_x[[1+x0]][[2]] * pz_xy[[1+x0]][[2]])
  
  ods_x1 <- py1_x1z / (1 - py1_x1z)
  ods_x0 <- py1_x0z / (1 - py1_x0z)
  
  list(
    ATE = sum((py1_x1z - py1_x0z) * pz),
    beta = log ( sum(ods_x1 / ods_x0 * pz) ),
    pz = pz
  )
}

backward_solver <- function(X, Y, R, fi, verbose = FALSE, 
                            direct = FALSE, ...) {
  
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
  
  if (direct) {
    
    py_x1z <- pz[[x1+1]][[y1+1]] * pxy[[x1+1]][[y1+1]] / (
      pz[[x1+1]][[y1+1]] * pxy[[x1+1]][[y1+1]] + pz[[x1+1]][[y0+1]] * pxy[[x1+1]][[y0+1]]
    )
    
    py_x0z <- pz[[x0+1]][[y1+1]] * pxy[[x0+1]][[y1+1]] / (
      pz[[x0+1]][[y1+1]] * pxy[[x0+1]][[y1+1]] + pz[[x0+1]][[y0+1]] * pxy[[x0+1]][[y0+1]]
    )
    
    py_x1z[is.nan(py_x1z)] <- 0
    py_x0z[is.nan(py_x0z)] <- 0
    ate <- sum((py_x1z - py_x0z) * pz_inf)
    
    return(
      list(pz = pz_inf, ATE = ate, solver = "backward-direct")
    )
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

logreg_ATE <- function(X, Y, Z, weights, solver = "expfam-1d", Sigma) {
  
  dat <- data.frame(X, Y, Z)
  
  # bootstrap the thingy
  dat <- dat[sample.int(nrow(dat), replace = TRUE, prob = weights), ]
  tryCatch(
    mod <- glm(
      Y ~ ., data = dat, family = "binomial"
    ),
    warning = function(w) {
      
      # Print the warning message
      message("Warning message: ", conditionMessage(w))
      
      # Enter browser mode for interactive debugging
      browser()
    }
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
    solver = solver,
    Sigma = Sigma,
    beta = mod$coefficients["X"]
  )
}


logreg_ATE_flex <- function(X, Y, Z, weights, solver = "expfam-1d", Sigma) {
  
  mc_samp <- 10^4
  wgh <- c(weights[1], weights[1 + mc_samp], weights[1 + 2 * mc_samp], 
           weights[1 + 3 * mc_samp])
  
  get_boot_samp <- function(Sigma, wgh, nsamp) {
    
    Z_mc <- X_mc <- Y_mc <- NULL
    for (x in c(0, 1)) {
      
      for (y in c(0, 1)) {
        
        nsamp_xy <- round(nsamp * wgh[1 + 2 * x + y])
        Z_mc <- rbind(Z_mc, z_2d(Sigma = Sigma[[1+x]][[1+y]],
                                 n_samp = nsamp_xy))
        X_mc <- c(X_mc, rep(x, nsamp_xy))
        Y_mc <- c(Y_mc, rep(y, nsamp_xy))
      }
    }
    
    dat <- data.frame(X = X_mc, Y = Y_mc, Z = Z_mc)
  }
  
  # nsamp_seq <- c(5, 10, 25, 50, 100) * 1000
  # df <- data.frame(nsamp = nsamp_seq)
  # df$coef <- df$std <- 0
  # for (i in seq_along(nsamp_seq)) {
  #   
  #   dat <- get_boot_samp(Sigma, wgh, nsamp_seq[i])
  #   mod <- glm(
  #     Y ~ ., data = dat, family = "binomial"
  #   )
  #   
  #   df$coef[i] <- coef(mod)["X"]
  #   df$std[i] <- summary(mod)$coefficients["X", "Std. Error"]
  #   df$nsamp[i] <- nsamp_seq[i]
  # }
  # 
  # browser()
  # 
  # ggplot(df, aes(x = log(nsamp), y = coef)) +
  #   geom_point() + geom_line() +
  #   geom_ribbon(aes(ymin = coef - 1.96 * std, ymax = coef + 1.96 * std),
  #               alpha = 0.2, linewidth = 0) +
  #   theme_minimal()
  # browser()
  
  # bootstrap the thingy
  dat <- get_boot_samp(Sigma, wgh, 10^6) # dat[sample.int(nrow(dat), replace = TRUE, prob = weights), ]
  tryCatch(
    mod <- glm(
      Y ~ ., data = dat, family = "binomial"
    ),
    warning = function(w) {
      
      # Print the warning message
      message("Warning message: ", conditionMessage(w))
      
      # Enter browser mode for interactive debugging
      browser()
    }
  )
  OR_std <- summary(mod)$coefficients["X", "Std. Error"]
  OR_hat <- mod$coefficients["X"]
  list(OR = exp(OR_hat), Sigma = Sigma,
       OR_lwr = exp(OR_hat - 1.96 * OR_std), 
       OR_upr = exp(OR_hat + 1.96 * OR_std))
}

two_stage_em <- function(X, Y, R, fi, mc_per_samp) {
  
  n_epoch <- 10
  # Step (0) - infer Sigma
  lmb <- list()
  
  Sigma <- infer_Sigma(X, Y, R, fi)
  
  pz <- cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  
  # Step (0+) - infer initial lambda, mu, beta
  fit_lmb <- function(X, Y, R) {
    
    xfit <- glm(X ~ ., data = data.frame(X = X, R = R), family = "binomial")
    lambda0 <- xfit$coefficients[-1]
    lambda_icept <- xfit$coefficients[1]
    
    yfit <- glm(Y ~ ., data = data.frame(Y = Y, X = X, R = R), family = "binomial")
    r_coeff <- grepl("^R", names(yfit$coefficients))
    
    mu_icept <- yfit$coefficients[1]
    mu0 <- yfit$coefficients[r_coeff]
    beta0 <- yfit$coefficients["X"]
    
    list(lambda = lambda0, mu = mu0, beta = beta0,
         lambda_icept = lambda_icept, mu_icept = mu_icept)
  }
  
  lmb[[1]] <- fit_lmb(X, Y, R)
  
  for (i in seq_len(n_epoch)) {
    
    # get px_z, py_z
    
    px_z <- expit(cpp_p_z(lmb[[i]]$lambda) + lmb[[i]]$lambda_icept)
    logity_z <- cpp_p_z(lmb[[i]]$mu) + lmb[[i]]$mu_icept
    beta <- lmb[[i]]$beta
    
    # sample mc_samples for each actual sample
    mc_twist <- lapply(
      seq_len(nrow(R)),
      function(k) {
        
        x <- X[k]
        y <- Y[k]
        r <- R[k, ]
        
        pr_zxy <- cpp_bit_fi_hash(r, fi[[1 + x]][[1 + y]])
        
        if (x == 1) px <- px_z else px <- 1 - px_z
        
        py_z <- expit(logity_z + beta * x)
        if (y == 1) py <- py_z else py <- 1 - py_z
        
        probs <- pz * px * py * pr_zxy
        probs <- probs / sum(probs)
        
        r_mc_idx <- sample(length(probs), size = mc_per_samp, prob = probs,
                           replace = TRUE)
        r_mc <- cpp_idx_to_bit(r_mc_idx - 1, length(r)) # returns a matrix
        
        cbind(X = x, Y = y, R = r_mc)
      }
    )
    
    mc_twist <- as.data.frame(do.call(rbind, mc_twist))
    
    # re-fit for lambda, mu, beta
    lmb[[i+1]] <- fit_lmb(mc_twist[, 1], mc_twist[, 2], mc_twist[, -c(1, 2)])
    cat("\r", i)
  }
  
  px_z <- expit(cpp_p_z(lmb[[i+1]]$lambda) + lmb[[i+1]]$lambda_icept)
  logity_z <- cpp_p_z(lmb[[i+1]]$mu) + lmb[[i+1]]$mu_icept
  beta <- lmb[[i+1]]$beta
  
  list(
    pz = pz,
    ATE = sum( (expit(logity_z + beta) - expit(logity_z)) * pz),
    solver = "two-stage-em",
    Sigma = Sigma,
    beta = beta
  )
}
