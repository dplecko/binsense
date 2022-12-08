
naive_bayes <- function(dat, fi, cmb = setdiff(names(dat), c("bmi_bin", "death"))) {
  
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

spcol_bayes <- function(dat, fi, cmb = setdiff(names(dat), c("bmi_bin", "death"))) {
  
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

model_bayes <- function(dat, fi, cmb = setdiff(names(dat), c("bmi_bin", "death")),
                        niter = 1L, pz_ground = NULL, Z_ground = NULL,
                        seed = 22, gibbs = 10, gibbs_slope = 0, tail_samples = 1,
                        display = FALSE) {
  
  mod <- function(x, y) {
    glm(as.formula(paste(names(y), "~ .")), 
        cbind(y, x), family = binomial())
  }
  
  pred <- function(mod, x) {
    predict(mod, x, type = "response")
  }
  
  add_cmb <- function(mod, col, aux_data, u, bmi_bin, death, fi) {
    
    pi <- fii <- rep(0, length(col))
    for (x in c(T, F)) {
      
      for (y in c(T, F)) {
        
        # subset on bmi, death
        idx <- bmi_bin == x & death == y
        
        # build model and get predictions
        pi[idx] <- pred(mod[[1+x]][[1+y]], aux_data[idx])
        fii[idx] <- fi[[x+1]][[y+1]]
      }
    }
    
    prob <- pi * fii / (1 - pi + pi * fii)
    add <- runif(length(u)) < prob # allow variability (for Gibbs)
    list(
      update = as.integer(col | add),
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
  bmi_bin <- lapply(seed, function(x) dat$bmi_bin)
  death <- lapply(seed, function(x) dat$death)
  
  Ui <- do.call(rbind, Ui)
  Rt <- do.call(rbind, Rt)
  bmi_bin <- do.call(c, bmi_bin)
  death <- do.call(c, death)
  
  if (!is.null(Z_ground)) {

    Z_ground <- as.data.table(Z_ground)
    Z_ground <- setnames(Z_ground, names(Z_ground), names(Rt))
  }
  
  thresh <- beta_hat <- ate_hat <- nZ_hat <- fit <- adv_acc <- pz <- list()
  cnt <- 1L
  for (nit in seq_len(niter)) {
    
    # train adversarial forest to predict which is which
    if (!is.null(Z_ground)) {
      
      adv <- rbind(Z_ground, Rt)
      adv <- cbind(adv, label = c(rep(0, nrow(Z_ground)), rep(1, nrow(Rt))))
      adv_rf <- ranger(label ~ ., data = adv)
      adv_acc[[nit]] <- mean(abs(adv$label - adv_rf$predictions) < 0.5)
    }
    
    pz[[nit]] <- tabulate(rowSums(Rt * enc_mat) + 1, nbins = 2^length(cmb))
    pz[[nit]] <- pz[[nit]] / sum(pz[[nit]])
    
    if (!is.null(pz_ground) & display) {
      cat("Iter", nit-1, "->", round(100 * sum(abs(pz[[nit]] - pz_ground)^2), 3), 
          "\n")
    }
    
    # update the model (M-step) / or initialize in first round
    #'* Currently taking the MLE over a single Gibbs sample *
    #'* -> need to take over a tail (after burn-in period) *
    #'* Need to compare stability -> how to do it? *
    if (nit == 1L) {
      psi_dat <- Rt
      num_rep <- 1L
    }
    for (i in seq_len(ncol(Rt))) {
      
      fiti <- list(list(0, 0), list(0, 0))
      
      for (x in c(T, F)) {
        
        for (y in c(T, F)) {
          
          idx <- rep(bmi_bin == x & death == y, num_rep)
          
          #'* expand the previous_steps list *
          fiti[[1+x]][[1+y]] <- mod(psi_dat[idx, -i, with=FALSE], 
                                    psi_dat[idx, i, with=FALSE])
        }
      }
      fit[[i]] <- fiti
    }
    
    # Gibbs sampling until convergence (Monte Carlo for E-step)
    prev_step <- list()
    for (gb in seq_len(gibbs + gibbs_slope * (nit - 1))) {
      for (i in seq_len(ncol(Rt))) {
        
        # need to do several iterations of Gibbs before refitting
        step <- add_cmb(fit[[i]], dat[[cmb[i]]], Rt[, -i, with=FALSE], Ui[, i], 
                        bmi_bin, death, fi)
        Rt[, c(cmb[i]) := step$update]
      }
      #'* save the different steps in a previous_steps list! *
      prev_step[[gb]] <- copy(Rt)
    }
    
    # estimate Beta_hat
    psi_dat <- do.call(rbind, tail(prev_step, n = tail_samples))
    num_rep <- nrow(psi_dat) / nrow(Rt)
    assert_that(num_rep == round(num_rep), msg = "Non-integer repetition value.")
    
    fit_dat <- cbind(death = rep(death, num_rep), psi_dat, 
                     bmi_bin = rep(bmi_bin, num_rep))
    fit <- glm(death ~ ., family = binomial(), data = fit_dat)
    
    fit_dat$bmi_bin <- FALSE
    y_b0 <- predict(fit, fit_dat)
    
    fit_dat$bmi_bin <- TRUE
    y_b1 <- predict(fit, fit_dat)
    
    beta_hat[[cnt]] <- fit$coefficients["bmi_binTRUE"]
    ate_hat[[cnt]] <- mean(y_b1) - mean(y_b0)
    nZ_hat[[cnt]] <- colSums(Z_ground) - colSums(Rt)
    cnt <- cnt + 1
    
  }
  
  list(
    pz = pz,
    thresholds = thresh,
    beta = beta_hat,
    ATE = ate_hat,
    nZ_diff = nZ_hat,
    adverse = adv_acc,
    seed = seed
  )
  
}
