
synth_data_mid <- function(n = 10^4, k = 5, seed = 22, 
                           fi = list(list(0.3, 0.3), list(0.3, 0.3)),
                           class = c("latent_u", "expfam-1d", "expfam-2d",
                                     "expfam-2d*"),
                           Sigma = NULL,
                           lam = NULL, mu = NULL, 
                           icept_x = NULL, icept_y = NULL,
                           beta = NULL) {
  
  
  class <- match.arg(class, c("latent_u", "expfam-1d", "expfam-2d", 
                              "expfam-2d*"))
  assertthat::assert_that(class %in% c("latent_u", "expfam-1d", "expfam-2d", 
                                       "expfam-2d*"),
                          msg = "Unknown model class specified.")
  
  set.seed(seed)
  
  Z <- switch(class,
     latent_u = z_latent_u(k, n, seed),
     `expfam-1d` = z_1d(k, n, seed),
     `expfam-2d` = z_2d(k, n, seed, Sigma),
     `expfam-2d*` = z_2d(k, n, seed, Sigma, Sigma_adv = TRUE)
  )
  
  use_icept <- TRUE # whether to use a specified intercept
  if (is.null(beta) | is.null(lam) | is.null(mu)) {
    
    lam <- -rep(2/3, k) # comorbidities make obesity less likely!
    mu <- rep(2/3, k) # comorbidities make mortality more likely!
    beta <- -0.15 # obesity protects!
    use_icept <- FALSE
  }
  
  u <- runif(n)
  
  if (use_icept) {
    
    X <- rbinom(n, 1, expit(Z %*% lam + icept_x))
    Y <- as.integer(u < expit(Z %*% mu + X * beta + icept_y))
    Y0 <- u < expit(Z %*% mu + 0 * beta + icept_y)
    Y1 <- u < expit(Z %*% mu + 1 * beta + icept_y)
  } else {
    
    X <- rbinom(n, 1, expit(Z %*% lam - sum(lam) / 2))
    
    Y <- as.integer(u < expit(Z %*% mu + X * beta - sum(mu) / 2))
    Y0 <- u < expit(Z %*% mu + 0 * beta - sum(mu) / 2)
    Y1 <- u < expit(Z %*% mu + 1 * beta - sum(mu) / 2)
  }
  
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
  fl_pth <- file.path(root, "data", paste0(src, "_real_data.rda"))
  
  if (!file.exists(fl_pth)) {
    
    data <- load_concepts(c("death", "bmi_all", "elixhauser", "age"), src, 
                          id_type = "patient")
    save(data, file = fl_pth)
  } else {
    
    load(fl_pth)
  }
  
  data <- replace_na(data, val = FALSE, vars = "death")
  data[, c(index_var(data)) := NULL]
  data[, Obesity := NULL]
  data <- data[complete.cases(data)]
  
  cmb_sel <- c("CHF", "Arrhythmia", "PVD", "Paralysis", "NeuroOther", "Pulmonary", 
               "Liver", "Lymphoma", "Mets", "Coagulopathy", "FluidsLytes")
  
  R <- as.matrix(data[, cmb_sel, with=FALSE])
  X <- as.integer(data$bmi_all > 25)
  Y <- as.integer(data$death)
  
  list(R = R, X = X, Y = Y)
}




