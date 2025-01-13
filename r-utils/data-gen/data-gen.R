
synth_data_mid <- function(n = 10^4, k = 5, seed = 22, 
                           fi = list(list(0.3, 0.3), list(0.3, 0.3)),
                           method = c("IM", "ZINF"),
                           A = NULL,
                           class = c("latent_u", "expfam-1d", "expfam-2d",
                                     "expfam-2d*"),
                           Sigma = NULL,
                           lam = NULL, mu = NULL, 
                           icept_x = NULL, icept_y = NULL,
                           beta = NULL) {
  
  method <- match.arg(method, c("IM", "ZINF"))
  class <- match.arg(class, c("latent_u", "expfam-1d", "expfam-2d", 
                              "expfam-2d*"))
  assertthat::assert_that(class %in% c("latent_u", "expfam-1d", "expfam-2d", 
                                       "expfam-2d*"),
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
  pz_xy <- pr_xy <- lxy <- list(list(NULL, NULL), list(NULL, NULL))
  
  if (!is.null(A) & method != "IM") { # code for generalized missingness
    
    for (x in c(T, F)) {
      
      for (y in c(T, F)) {
        
        idx <- X == x & Y == y
        for (i in seq_len(k)) {
          
          # get the cpp bit to idx, idx to bit
          z_idx <- 1 + binsensate:::cpp_bit_to_idx_mat(Z[idx, ])
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
          lxy[[1+x]][[1+y]] <- logit(colMeans(Z[idx, , drop = FALSE]))
        }
      }
    }
  }
  
  attr(R, "lambda") <- lam
  attr(R, "mu") <- mu
  
  enc_mat <- t(replicate(nrow(Z), 2^(seq_len(k) - 1L)))
  if (k == 1) enc_mat <- t(enc_mat)
  pz <- tabulate(rowSums(Z * enc_mat) + 1, nbins = 2^k)
  pz <- pz / sum(pz)
  
  for (x in c(T, F)) {
    
    for (y in c(T, F)) {
      
      idx <- X == x & Y == y
      pz_xy[[x+1]][[y+1]] <- tabulate(rowSums(Z[idx, , drop = FALSE] * 
                                                 enc_mat[idx, , drop = FALSE]) + 1, 
                                       nbins = 2^k)
      pz_xy[[x+1]][[y+1]] <- pz_xy[[x+1]][[y+1]] / sum(pz_xy[[x+1]][[y+1]])
      
      pr_xy[[x+1]][[y+1]] <- tabulate(rowSums(R[idx, , drop = FALSE] * 
                                                enc_mat[idx, , drop = FALSE]) + 1, 
                                      nbins = 2^k)
      pr_xy[[x+1]][[y+1]] <- pr_xy[[x+1]][[y+1]] / sum(pr_xy[[x+1]][[y+1]])
    }
  }
  
  list(
    X = X, Y = Y, R = R, Z = Z, pz = pz, beta = beta, ATE = ate, 
    pz_xy = pz_xy, pr_xy = pr_xy,
    lxy = lxy
  )
}

real_data <- function(src = c("miiv", "mimic_demo"), n = Inf, k = Inf) {
  
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
  
  cmb_sel <- c("FluidsLytes", "Mets", "Liver", "Coagulopathy", "Arrhythmia", 
               "NeuroOther", "Paralysis", "Lymphoma", "CHF", "PVD", "Pulmonary")
  
  k <- min(k, length(cmb_sel))
  n <- min(n, nrow(data))
  
  R <- as.matrix(data[1:n, cmb_sel[1:k], with=FALSE])
  X <- as.integer(data[1:n]$bmi_all > 25)
  Y <- as.integer(data[1:n]$death)
  
  list(R = R, X = X, Y = Y)
}




