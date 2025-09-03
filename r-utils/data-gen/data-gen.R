
zw_2d <- function(k, n_samp, seed = 2022, SLO) {
  
  set.seed(seed)
  
  d_index <- k + seq_len(ncol(SLO) - k)
  Sigma <- SLO[-d_index, -d_index]
  Sigma0 <- solve(SLO[d_index, d_index, drop = FALSE])
  Lambda <- SLO[d_index, -d_index, drop = FALSE]
  
  pz <- binsensate:::theta_Ez(SLO, d_index)
  z <- sample.int(2^k, size = n_samp, prob = pz, replace = TRUE)
  z <- binsensate:::cpp_idx_to_bit(z-1, k)
  w <- z %*% t(Sigma0 %*% Lambda) + 
    MASS::mvrnorm(n_samp, mu = rep(0, length(d_index)), Sigma = Sigma0)
  cbind(z, w)
}

synth_data_mid <- function(n = 10^4, k = 5, seed = 2025, 
                           fi = list(list(0.3, 0.3), list(0.3, 0.3)),
                           method = c("IF", "ZINF"),
                           A = NULL,
                           class = c("latent_u", "expfam-1d", "expfam-2d",
                                     "expfam-2d*", "expfam-2d-cts"),
                           Sigma = NULL,
                           SLO = NULL,
                           lam = NULL, mu = NULL, 
                           icept_x = NULL, icept_y = NULL,
                           beta = NULL, ate = NULL) {
  
  method <- match.arg(method, c("IF", "ZINF"))
  class <- match.arg(class, c("latent_u", "expfam-1d", "expfam-2d", 
                              "expfam-2d*", "expfam-2d-cts"))
  assertthat::assert_that(class %in% c("latent_u", "expfam-1d", "expfam-2d", 
                                       "expfam-2d*", "expfam-2d-cts"),
                          msg = "Unknown model class specified.")
  
  set.seed(seed)
  
  if (!is.null(A)) {
    
    fi <- NULL
    method <- "matrix"
    if (is.matrix(A)) A <- list(list(A, A), list(A, A))
  } else {
    
    # if agnostic fi, create a list
    if (!is.list(fi)) fi <- list(list(fi, fi), list(fi, fi))
    fi <- binsensate:::expand_fi(fi, k, method)
    A <- binsensate:::fi_to_A(fi, k, method)
  }
  
  Z <- switch(class,
              latent_u = z_latent_u(k, n, seed),
              `expfam-1d` = z_1d(k, n, seed),
              `expfam-2d` = z_2d(k, n, seed, Sigma),
              `expfam-2d*` = z_2d(k, n, seed, Sigma, Sigma_adv = TRUE),
              `expfam-2d-cts` = zw_2d(k, n, seed, SLO)
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
  
  R <- Z[, seq_len(k)]
  pz_xy <- pr_xy <- lxy <- list(list(NULL, NULL), list(NULL, NULL))
  
  if (!is.null(A) & method != "IF") { # code for generalized missingness
    
    for (x in c(T, F)) {
      
      for (y in c(T, F)) {
        
        idx <- X == x & Y == y
        for (i in seq_len(k)) {
          
          # get the cpp bit to idx, idx to bit
          z_idx <- 1 + binsensate:::cpp_bit_to_idx_mat(Z[idx, seq_len(k)])
          r_idx <- vapply(
            z_idx, 
            function(z_id) {
              
              sample.int(n = 2^k, size = 1, prob = A[[1+x]][[1+y]][, z_id])
            }, numeric(1L)
          )
          R[idx, ] <- binsensate:::cpp_idx_to_bit(r_idx-1, k)
        }
      }
    }
  } else {
    
    for (x in c(T, F)) {
      
      for (y in c(T, F)) {
        
        idx <- X == x & Y == y
        for (i in seq_len(k)) {
          
          R[idx, i] <- Z[idx, i] * rbinom(sum(idx), 1, prob = 1 - fi[[1+x]][[1+y]][i])
        }
        lxy[[1+x]][[1+y]] <- logit(colMeans(Z[idx, seq_len(k), drop = FALSE]))
      }
    }
  }
  
  attr(R, "lambda") <- lam
  attr(R, "mu") <- mu
  
  enc_mat <- t(replicate(nrow(Z[, seq_len(k)]), 2^(seq_len(k) - 1L)))
  if (k == 1) enc_mat <- t(enc_mat)
  pz <- tabulate(rowSums(Z[, seq_len(k)] * enc_mat) + 1, nbins = 2^k)
  pz <- pz / sum(pz)
  
  for (x in c(T, F)) {
    
    for (y in c(T, F)) {
      
      idx <- X == x & Y == y
      pz_xy[[x+1]][[y+1]] <- tabulate(rowSums(Z[idx, seq_len(k), drop = FALSE] * 
                                                 enc_mat[idx, , drop = FALSE]) + 1, 
                                       nbins = 2^k)
      pz_xy[[x+1]][[y+1]] <- pz_xy[[x+1]][[y+1]] / sum(pz_xy[[x+1]][[y+1]])
      
      pr_xy[[x+1]][[y+1]] <- tabulate(rowSums(R[idx, , drop = FALSE] * 
                                                enc_mat[idx, , drop = FALSE]) + 1, 
                                      nbins = 2^k)
      pr_xy[[x+1]][[y+1]] <- pr_xy[[x+1]][[y+1]] / sum(pr_xy[[x+1]][[y+1]])
    }
  }
  
  dat <- list(
    X = X, Y = Y, R = R, Z = Z[, seq_len(k)], pz = pz, beta = beta, ATE = ate, 
    pz_xy = pz_xy, pr_xy = pr_xy,
    lxy = lxy
  )
  
  if (!is.null(SLO)) dat[["W"]] <- Z[, -seq_len(k)]
  dat
}

real_data <- function(src = c("miiv", "mimic_demo"), n = Inf, k = Inf,
                      trim_ends = TRUE, add_W = FALSE, ret_dt = FALSE,
                      ngroups = NULL) {
  
  src <- match.arg(src, c("miiv", "mimic_demo"))
  fl_pth <- file.path(root, "data", paste0(src, "_real_data2.rda"))
  
  if (!file.exists(fl_pth)) {
    
    data <- load_concepts(c("death", "bmi_all", "elix_wide", "age"), src, 
                          id_type = "patient")
    save(data, file = fl_pth)
  } else {
    
    load(fl_pth)
  }
  
  # data2 <- load_concepts(c("death", "bmi_all", "elix_wide", "age"), src, 
  #                       id_type = "patient")
  # iid <- intersect(id_col(data), id_col(data2))
  # data <- data[subject_id %in% iid, c(id_vars(data), names(icd9_elix_map)[-7]), with=F]
  # data2 <- data2[subject_id %in% iid, c(id_vars(data), names(icd9_elix_map)), with=F]
  # 
  # data <- replace_na(data, val = 0)
  # data2 <- replace_na(data2, val = 0)
  # 
  # for (col in names(data[, -c(1)])) {
  #   
  #   cat(col, "\n")
  #   print(table(data[[col]], data2[[col]]))
  #   cat("\n")
  # }
  
  data <- data[age >= 18]
  data <- data[is.na(death), death := FALSE]
  data[, c(meta_vars(data), "elix_wide") := NULL]
  data <- data[complete.cases(data)]

  # remove Obesity and Weight Loss as they are closely related to BMI
  data[, c("Obesity", "WeightLoss") := NULL]
  
  # combine DM and DMcx
  data[, DM := pmax(DM, DMcx)]
  data[, DMcx := NULL]
  
  # remove extremes (underweight < 18.5, class II/III obese > 35)
  if (trim_ends) data <- data[bmi_all >= 18.5 & bmi_all < 35]
  
  if (ret_dt) return(data)
   
  # regress
  logreg <- glm(death ~ . - age - bmi_all + (bmi_all >= 25), data = data, 
                family = "binomial")
  lr_coef <- summary(logreg)$coefficients
  lr_coef <- lr_coef[order(lr_coef[, 1], decreasing = TRUE), ]
  sel1 <- lr_coef[, 1] > 0 & lr_coef[, 4] < 0.05
  sel2 <- lr_coef[, 1] > 0
  cmb_sel1 <- rownames(lr_coef)[sel1]
  cmb_sel2 <- rownames(lr_coef)[sel2]
  
  cmb_sel <- cmb_sel1
  
  n <- min(n, nrow(data))
  if (!is.null(ngroups)) {
    
    # 5 groups
    g5 <- list(
      CardioVascular = c("Arrhythmia", "CHF", "PVD", "PHTN"),
      RenalMetabolic = c("Renal", "FluidsLytes", "Liver", "Coagulopathy"),
      Pulmonary      = c("Pulmonary"),
      Neuro          = c("Paralysis", "NeuroOther"),
      Malignancy     = c("Tumor", "Mets", "Lymphoma")
    )
    
    # 6 groups
    g6 <- list(
      CardioVascular = c("Arrhythmia", "CHF", "PVD", "PHTN"),
      Renal          = c("Renal", "FluidsLytes"),
      Pulmonary      = c("Pulmonary"),
      Neuro          = c("Paralysis", "NeuroOther"),
      Malignancy     = c("Tumor", "Mets", "Lymphoma"),
      LiverCoag      = c("Liver", "Coagulopathy")
    )
    
    grp <- if (ngroups == 6) g6 else g5
    grp_enc <- lapply(grp, function(x) as.integer(rowMeans(data[, x, with=FALSE]) > 0))
    R <- do.call(cbind, grp_enc)[1:n, ]
  } else {
    
    k <- min(k, length(cmb_sel))
    R <- as.matrix(data[1:n, cmb_sel[1:k], with=FALSE])[1:n, ]
  }
  
  X <- as.integer(data[1:n]$bmi_all >= 25)
  Y <- as.integer(data[1:n]$death)
  W <- as.matrix(as.integer(data[1:n]$age))
  
  ret <- list(R = R, X = X, Y = Y)
  if (add_W) ret[["W"]] <- (W - mean(W)) / sd(W)
  
  ret
}

