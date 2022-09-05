
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

dat <- synth_data()[["dat"]]
fixy <- list(list(0, 0), list(0, 0))
fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3

naive_bayes(dat, fixy)
# conclusion: naive bayes fails when Z1,...,Zk are not independent.

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

dat <- synth_data()[["dat"]]
fixy <- list(list(0, 0), list(0, 0))
fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3
spcol_bayes(dat, fixy)

model_bayes <- function(dat, fi, cmb = setdiff(names(dat), c("bmi", "death")),
                        niter = 1L, pz_ground = NULL) {
  
  mod_and_pred <- function(x, y) {
    glm(y ~ ., cbind(y = y, x = x), family = binomial())$fitted.values
  }
  
  add_cmb <- function(col, aux_data, u, bmi_bin, death, fi) {
    
    pi <- fii <- rep(0, length(col))
    for (x in c(T, F)) {
      
      for (y in c(T, F)) {
        
        # subset on bmi, death
        idx <- bmi_bin == x & death == y
        
        # build model and get predictions
        pi[idx] <- mod_and_pred(aux_data[idx], col[idx])
        fii[idx] <- fi[[x+1]][[y+1]]
      }
    }
    
    prob <- pi * fii / (1 - pi + pi*fii)
    add <- u < prob
    list(
      update = col | add,
      thresh = prob
    )
  }
  
  enc_mat <- t(replicate(nrow(dat), 2^(seq_along(cmb) - 1L)))
  Ui <- array(runif(length(enc_mat)), dim = dim(enc_mat))
  Rt <- copy(dat[, c(cmb), with=FALSE])
  bmi_bin <- dat$bmi < 25
  death <- dat$death
  thresh <- list()
  for (nit in seq_len(niter)) {
    
    pz <- table(rowSums(Rt * enc_mat))
    pz <- pz / sum(pz)
    
    if (!is.null(pz_ground)) {
      cat("Iter", nit-1, "->", round(100 * sum(abs(pz - pz_ground)^2), 3), "\n")
    }
    
    for (i in seq_len(ncol(Rt))) {
      
      step <- add_cmb(dat[[cmb[i]]], Rt[, -i, with=FALSE], Ui[, i], 
                      bmi_bin, death, fi)
      Rt[, c(cmb[i]) := step$update]
      
      
      if (i == 1) {
        
        thresh[[nit]] <- step$thresh
        if (nit > 1) {
          
          cat("Threshold L2 change", sum((thresh[[nit]] - thresh[[nit-1]])^2), 
              "\n") 
        }
      }
    }
  }
  
  for (i in seq_along(cmb)) dat[, c(cmb[i]) := Rt[[i]]]
  
  cat("Recovered P(z) = ", pz)
  
} 

# k = 2 (manual) example
# dat <- synth_data()[["dat"]]
# fixy <- list(list(0, 0), list(0, 0))
# fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3
# model_bayes(dat, fixy)

# k = 5 (u-model) example
aux <- synth_data_mid()
dat <- aux[["dat"]]

fixy <- list(list(0, 0), list(0, 0))
fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3
model_bayes(dat, fixy, pz_ground = aux[["pz"]], niter = 15L)




