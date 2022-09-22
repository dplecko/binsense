
naive_bayes <- function(dat, fi, cmb = setdiff(names(dat), c("bmi", "death"))) {
  
  add_cmb <- function(col, bmi_bin, death, fi) {
    
    pi <- fii <- rep(0, length(col))
    for (x in c(T, F)) {
      for (y in c(T, F)) {
        idx <- bmi_bin == x & death == y
        pi[idx] <- mean(col[idx]) / (1 - fi[[x+1]][[y+1]])
        fii[idx] <- fi[[x+1]][[y+1]]
      }
    }
    
    prob <- pi * fii / (1 - pi + pi*fii)
    
    add <- rbinom(n = length(col), size = 1, prob = prob)
    col | add
  }
  
  for (col in cmb) {
    dat[, c(col) := add_cmb(get(col), get("bmi") < 25, get("death"), fi)]
  }
  
  enc_mat <- t(replicate(nrow(dat), 2^(seq_along(cmb) - 1L)))
  
  pz <- table(rowSums(dat[, cmb, with=FALSE] * enc_mat))
  pz <- pz / sum(pz)
  cat("Recovered P(z) = ", pz)
  
} 

spcol_bayes <- function(dat, fi, cmb = setdiff(names(dat), c("bmi", "death"))) {
  
  add_cmb <- function(col, spcol, bmi_bin, death, fi) {
    
    pi <- fii <- rep(0, length(col))
    for (x in c(T, F)) {
      
      for (y in c(T, F)) {
        
        for (sp in c(0, 1)) {
          
          idx <- bmi_bin == x & death == y & spcol == sp
          pi[idx] <- mean(col[idx]) / (1 - fi[[x+1]][[y+1]])
          fii[idx] <- fi[[x+1]][[y+1]]
        }
      }
    }
    
    prob <- pi * fii / (1 - pi + pi*fii)
    
    
    
    add <- rbinom(n = length(col), size = 1, prob = prob)
    col | add
  }
  
  for (col in cmb) {
    dat[, c(col) := add_cmb(get(col), get(setdiff(cmb, col)), 
                            get("bmi") < 25, get("death"), fi)]
  }
  
  enc_mat <- t(replicate(nrow(dat), 2^(seq_along(cmb) - 1L)))
  
  pz <- table(rowSums(dat[, cmb, with=FALSE] * enc_mat))
  pz <- pz / sum(pz)
  cat("Recovered P(z) = ", pz)
  
} 

model_bayes <- function(dat, fi, cmb = setdiff(names(dat), c("bmi", "death")),
                        niter = 1L, pz_ground = NULL, Z_ground = NULL,
                        seed = 22, display = FALSE) {
  
  mod_and_pred <- function(x, y, fid) {
    glm(y ~ ., cbind(y = y, x = x), family = binomial())$fitted.values / (1-fid)
  }
  
  add_cmb <- function(col, aux_data, u, bmi_bin, death, fi) {
    
    pi <- fii <- rep(0, length(col))
    for (x in c(T, F)) {
      
      for (y in c(T, F)) {
        
        # subset on bmi, death
        idx <- bmi_bin == x & death == y
        
        # build model and get predictions
        pi[idx] <- mod_and_pred(aux_data[idx], col[idx], fid = fi[[x+1]][[y+1]])
        fii[idx] <- fi[[x+1]][[y+1]]
      }
    }
    
    prob <- pi * fii / (1 - pi + pi * fii)
    add <- u < prob
    list(
      update = col | add,
      thresh = prob
    )
  }
  
  enc_mat <- t(replicate(nrow(dat), 2^(seq_along(cmb) - 1L)))
  
  Ui <- lapply(
    seed,
    function(sid) {
      with_seed(sid, array(runif(length(enc_mat)), dim = dim(enc_mat)))
    }
  )
  
  enc_mat <- t(replicate(nrow(dat) * length(seed), 2^(seq_along(cmb) - 1L)))

  Rt <- lapply(seed, function(x) copy(dat[, c(cmb), with=FALSE]))
  bmi_bin <- lapply(seed, function(x) dat$bmi > 25)
  death <- lapply(seed, function(x) dat$death)
  
  Ui <- do.call(rbind, Ui)
  Rt <- do.call(rbind, Rt)
  bmi_bin <- do.call(c, bmi_bin)
  death <- do.call(c, death)
  
  thresh <- beta_hat <- nZ_hat <- list()
  cnt <- 1L
  for (nit in seq_len(niter)) {
    
    pz <- table(rowSums(Rt * enc_mat))
    pz <- pz / sum(pz)
    
    if (!is.null(pz_ground) & display) {
      cat("Iter", nit-1, "->", round(100 * sum(abs(pz - pz_ground)^2), 3), "\n")
    }
    
    for (i in seq_len(ncol(Rt))) {
      
      step <- add_cmb(dat[[cmb[i]]], Rt[, -i, with=FALSE], Ui[, i], 
                      bmi_bin, death, fi)
      Rt[, c(cmb[i]) := step$update]
      
      # need to do several iterations of Gibbs before refitting
      fit <- glm(death ~ ., family = binomial(), 
                 data = cbind(death, Rt, bmi_bin))
      
      beta_hat[[cnt]] <- fit$coefficients["bmi_binTRUE"]
      nZ_hat[[cnt]] <- colSums(Z_ground) - colSums(Rt)
      cnt <- cnt + 1
      
      if (i == 1) {
        
        thresh[[nit]] <- step$thresh
        if (nit > 1 & display) {
          
          cat("Threshold L2 change", sum((thresh[[nit]] - thresh[[nit-1]])^2), 
              "\n") 
        }
      }
    }
  }
  
  # for (i in seq_along(cmb)) dat[, c(cmb[i]) := Rt[[i]]]
  
  list(
    pz = pz,
    thresholds = thresh,
    beta = beta_hat,
    nZ_diff = nZ_hat,
    seed = seed
  )
  
}
