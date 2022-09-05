
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

synth_data_mid <- function(n = 10^5, k = 5, seed = 22) {
  
  set.seed(seed)
  
  u <- runif(n)
  ulam <- rep(1/4, 5)
  
  Z <- c()
  for (i in seq_len(k)) {
    
    Z <- cbind(Z, rbinom(n, 1, prob = expit(ulam * u - 1)))
  }
  
  lam <- rep(1/2, k)
  mu <- rep(1/2, k)
  beta <- 0.1
  fi <- 0.3
  
  X <- rbinom(n, 1, expit(Z %*% lam - k/2))
  
  Y <- rbinom(n, 1, expit(Z %*% mu + X * beta - k/2))
  
  ate <- mean(expit(Z %*% mu + beta - k/2) - expit(Z %*% mu - k/2))
  
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
  pz <- table(rowSums(Z * enc_mat))
  pz <- pz / sum(pz)
  
  list(
    ate = ate, dat = dat, Z = Z, pz = pz
  )
}
