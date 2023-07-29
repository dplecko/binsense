
synth_data_mid <- function(n = 10^4, k = 5, seed = 22, 
                           fi = list(list(0.3, 0.3), list(0.3, 0.3)),
                           class = c("latent_u", "expfam-1d", "expfam-2d",
                                     "expfam-2d*")) {
  
  
  class <- match.arg(class, c("latent_u", "expfam-1d", "expfam-2d", 
                              "expfam-2d*"))
  assertthat::assert_that(class %in% c("latent_u", "expfam-1d", "expfam-2d", 
                                       "expfam-2d*"),
                          msg = "Unknown model class specified.")
  
  set.seed(seed)
  
  Z <- switch(class,
     latent_u = z_latent_u(k, n, seed),
     `expfam-1d` = z_1d(k, n, seed),
     `expfam-2d` = z_2d(k, n, seed),
     `expfam-2d*` = z_2d(k, n, seed, Sigma_adv = TRUE)
  )
  
  lam <- -rep(2/3, k) # comorbidities make obesity less likely!
  mu <- rep(2/3, k) # comorbidities make mortality more likely!
  beta <- -0.15 # obesity protects!
  
  X <- rbinom(n, 1, expit(Z %*% lam - sum(lam) / 2))
  
  u <- runif(n)
  Y <- as.integer(u < expit(Z %*% mu + X * beta - sum(mu) / 2))
  
  Y0 <- u < expit(Z %*% mu + 0 * beta - sum(mu) / 2)
  Y1 <- u < expit(Z %*% mu + 1 * beta - sum(mu) / 2)
  ate <- mean(Y1) - mean(Y0)
  
  R <- Z
  pz_all <- lxy <- list(list(NULL, NULL), list(NULL, NULL))
  for (x in c(T, F)) {
    
    for (y in c(T, F)) {
      
      idx <- X == x & Y == y
      for (i in seq_len(k)) {
        
        R[idx, i] <- Z[idx, i] * rbinom(sum(idx), 1, prob = 1 - fi[[1+x]][[1+y]])
        lxy[[1+x]][[1+y]] <- logit(colMeans(Z[idx, ]))
      }
    }
  }
  
  attr(R, "lambda") <- lam
  attr(R, "mu") <- mu
  
  enc_mat <- t(replicate(nrow(Z), 2^(seq_len(k) - 1L)))
  pz <- tabulate(rowSums(Z * enc_mat) + 1, nbins = 2^k)
  pz <- pz / sum(pz)
  
  for (x in c(T, F)) {
    
    for (y in c(T, F)) {
      
      idx <- X == x & Y == y
      pz_all[[x+1]][[y+1]] <- tabulate(rowSums(Z[idx, ] * enc_mat[idx, ]) + 1, 
                                       nbins = 2^k)
      pz_all[[x+1]][[y+1]] <- pz_all[[x+1]][[y+1]] / sum(pz_all[[x+1]][[y+1]])
    }
  }
  
  list(
    X = X, Y = Y, R = R, Z = Z, pz = pz, beta = beta, ATE = ate, pz_all = pz_all,
    lxy = lxy
  )
}

real_data <- function(src = c("miiv", "mimic_demo")) {
  
  src <- match.arg(src, c("miiv", "mimic_demo"))
  
  ind <- load_concepts(c("bmi", "death"), src)
  ind[, c(index_var(ind)) := NULL]
  ind <- ind[!is.na(bmi)]
  
  dat <- merge(ind, get_icd(index = "charlson"), all.x = TRUE)
  dat[is.na(death), death := FALSE]
  dat[, bmi_bin := bmi > 25]
  
  for (col in setdiff(names(dat), "death")) {
    dat[is.na(get(col)), c(col) := FALSE]
  }
  
  keep_vars <- c("stay_id", "bmi_bin", "death")
  cmb <- c("CHF", "Pulmonary", "LiverSevere", "Renal", "Stroke", "MI")
  
  dat <- dat[, c(keep_vars, cmb), with=FALSE]
  dat <- setnames(dat, cmb, paste("cmb", seq_along(cmb), sep = "_"))
  dat
}




