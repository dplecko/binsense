
synth_data <- function(n = 10^5, pz = c(0.5, 0.1, 0.1, 0.3), k = 2, seed = 22) {
  
  set.seed(seed)
  Z <- sample(seq_len(2^k)-1, n, prob = pz, replace = TRUE)
  
  lambda1 <- lambda2 <- 1
  mu1 <- mu2 <- 1
  beta <- 0.1
  fi <- 0.3
  
  
  Z <- lapply(Z, function(z) as.integer(intToBits(z)[1:k]))
  Z <- do.call(rbind, Z)
  
  X <- rbinom(n, 1, expit(Z %*% c(lambda1, lambda2) - 2))
  
  Y <- rbinom(n, 1, expit(Z %*% c(mu1, mu2) + X * beta - 2))
  
  ate <- mean(expit(Z %*% c(mu1, mu2) + beta - 2) - expit(Z %*% c(mu1, mu2) - 2))
  
  R <- Z
  for (i in seq_len(k)) {
    
    R[, i] <- Z[, i] * rbinom(n, 1, prob = 1 - fi)
    
  }
  
  dat <- data.frame(X, R, Y)
  dat$X <- ifelse(dat$X == 1, 30, 20)
  dat$Y <- ifelse(dat$Y == 1, TRUE, FALSE)
  names(dat) <- c("bmi", "CHF", "Pulmonary", "death")
  dat <- as.data.table(dat)
  
  list(
    ate = ate, dat = dat, Z = Z
  )
}

synth_data_mid <- function(n = 10^4, k = 5, seed = 22,
                           class = c("latent_u", "2nd_order_exp")) {
  
  
  class <- match.arg(class, c("latent_u", "2nd_order_exp"))
  assertthat::assert_that(class %in% c("latent_u", "2nd_order_exp"),
                          msg = "Unknown model class specified.")
  
  set.seed(seed)
  fi <- 0.3
  
  if (class == "latent_u") {
    
    u <- runif(n)
    ulam <- rep(1/4, 5)
    
    Z <- c()
    for (i in seq_len(k)) {
      
      Z <- cbind(Z, rbinom(n, 1, prob = expit(ulam * u - 1)))
    }
    
  } else if (class == "2nd_order_exp") {
    
    pz <- rep(0, 2^k)
    sig <- runif(k) - 2 / 3
    Sig <- matrix((runif(k^2) - 1/2) / 3 , nrow = k)
    
    pz <- vapply(
      seq_len(2^k),
      function(i) {
        z <- as.integer(intToBits(i-1)[seq_len(k)])
        exp(sum(z * sig) + t(z) %*% Sig %*% z)
      }, numeric(1L)
    )
    
    Z <- sample.int(2^k, size = n, prob = pz, replace = TRUE)
    Z <- lapply(Z, function(z) as.integer(intToBits(z-i)[seq_len(k)]))
    Z <- do.call(rbind, Z)
    
  }
  
  lam <- -rep(2/3, k) # comorbidities make obesity less likely!
  mu <- rep(2/3, k) # comorbidities make mortality more likely!
  beta <- -0.05 # obesity protects!
  
  X <- rbinom(n, 1, expit(Z %*% lam - sum(lam) / 2))
  
  Y <- rbinom(n, 1, expit(Z %*% mu + X * beta - sum(mu) / 2))
  
  # ate <- mean(expit(Z %*% mu + beta - 0.66 * k * mu[1]) - 
  #               expit(Z %*% mu - 0.5 * k * lam[1]))
  
  R <- Z
  for (i in seq_len(k)) {
    
    R[, i] <- Z[, i] * rbinom(n, 1, prob = 1 - fi)
    
  }
  
  dat <- data.frame(X, R, Y)
  dat$X <- ifelse(dat$X == 1, 30, 20)
  dat$Y <- ifelse(dat$Y == 1, TRUE, FALSE)
  names(dat) <- c("bmi", paste0("cmb_", seq_len(k)), "death")
  dat <- as.data.table(dat)
  
  enc_mat <- t(replicate(nrow(dat), 2^(seq_len(k) - 1L)))
  pz <- tabulate(rowSums(Z * enc_mat) + 1, nbins = 2^k)
  pz <- pz / sum(pz)
  
  list(
    dat = dat, Z = Z, pz = pz, beta = beta#, ate = ate 
  )
}

