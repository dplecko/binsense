
#' * 2d exponential family generators *
gen_Sigma <- function(k, adv = FALSE) {
  
  Sigma <- array(runif(k^2, -1, 1), dim = c(k, k))
  
  Sigma[upper.tri(Sigma)] <- 0
  
  if (adv) {
    
    Sigma[!(upper.tri(Sigma)) & !(lower.tri(Sigma))] <- 0
    Sigma <- Sigma * sign(Sigma)
    Sigma <- t(t(Sigma) * (-1)^(2:(k+1)))
  }
  
  (Sigma + t(Sigma) * upper.tri(t(Sigma)))
}

scaling_2d <- function(Sigma) {
  
  k <- ncol(Sigma)
  
  pz <- vapply(
    seq_len(2^k), 
    function(x) {
      x_bit <- to_bits(x-1, k)
      exp(sum(t(x_bit) %*% Sigma %*% x_bit))
    },
    numeric(1L)
  )
  
  pz
}

z_2d <- function(k, n_samp, seed = 2022, Sigma = NULL, mode = "cpp",
                 Sigma_adv = FALSE) {
  
  set.seed(seed)
  
  if (is.null(Sigma)) {
    
    Sigma <- gen_Sigma(k, adv = Sigma_adv)
  } else {
    
    k <- ncol(Sigma)
  }
  
  if (mode == "cpp") {
    pz <- scaling_2d(Sigma)
  } else {
    pz <- cpp_scaling_2d(Sigma)
  }
  
  z <- sample.int(2^k, size = n_samp, prob = pz, replace = TRUE)
  
  z <- do.call(rbind, lapply(z, function(zi) to_bits(zi-1, k = k)))
  attr(z, "Sigma") <- Sigma
  attr(z, "pz") <- pz
  
  z
}

r_2d <- function(k, n_samp, fi, Sigma = NULL, seed = 2022) {
  
  z <- z_2d(k, n_samp, seed, Sigma = Sigma)
  
  mask <- array(rbinom(k * n_samp, 1, 1-fi), dim = dim(z))
  
  r <- z * mask
  
  attr(r, "Sigma") <- attr(z, "Sigma")
  
  r
}

#' * 1d exponential family generators *
z_1d <- function(k = 3, n_samp, seed = 2022, lambda = NULL) {
  
  set.seed(seed)
  
  if (is.null(lambda)) {
    
    lambda <- runif(k, -1, 1)
  } else {
    
    k <- length(lambda)
  }
  
  pz <- scaling_1d(lambda)
  
  pz <- pz / sum(pz)
  
  # can speed up with binary search? -> not really a bottleneck
  z <- sample.int(2^k, size = n_samp, prob = pz, replace = TRUE) 
  
  z <- do.call(rbind, lapply(z, function(zi) to_bits(zi-1, k = k)))
  attr(z, "lambda") <- lambda
  
  z
}

scaling_1d <- function(lambda) {
  
  k <- length(lambda)
  pz <- vapply(
    seq_len(2^k), 
    function(x) exp(sum(to_bits(x-1, k) * lambda)),
    numeric(1L)
  )
  
  pz
}

r_1d <- function(k, n_samp, fi, seed = 2022) {
  
  z <- z_1d(k, n_samp, seed)
  
  mask <- array(rbinom(k * n_samp, 1, 1-fi), dim = dim(z))
  
  r <- z * mask
  
  attr(r, "lambda") <- attr(z, "lambda")
  
  r
}

#' * latent u generator *
z_latent_u <- function(k, nsamp, seed = 2022) {
  
  set.seed(seed)
  
  u <- runif(nsamp)
  ulam <- rep(1/4, 5)
  
  replicate(k, rbinom(nsamp, 1, expit(ulam * u - 1)))
}


