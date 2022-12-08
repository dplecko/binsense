
sensitivity_ATE <- function(fi, dat, order = Inf, verbose = FALSE) {
  
  cmb <- setdiff(names(dat), c("stay_id", "bmi_bin", "death"))
  k <- length(cmb)
  enc_mat <- t(replicate(nrow(dat), 2^(seq_along(cmb) - 1L)))
  
  # get P(y | r, x)
  logreg <- glm(
    as.formula(paste("factor(death) ~", paste(c("bmi_bin", cmb), 
                                              collapse = "+"))),
    data = dat, family = "binomial"
  )
  
  pyrx0 <- pyrx1 <- rep(0, 2^k)
  for (i in seq_len(2^k)) {
    
    i_bin <- as.integer(intToBits(i-1))[1:k]
    pyrx0[i] <- expit(sum(logreg$coefficients * c(1, 0, i_bin)))
    pyrx1[i] <- expit(sum(logreg$coefficients * c(1, 1, i_bin)))
    
  }
  
  pr <- pz <- pxy <- A <- A_inv <- B <- list(list(0, 0), list(0, 0))
  
  x0 <- y0 <- 0
  x1 <- y1 <- 1
  
  for (x in c(T, F)) {
    
    for (y in c(T, F)) {
      
      idx <- dat$bmi_bin == x & dat$death == y
      
      cprxy <- tabulate(
        rowSums(
          dat[idx, c(cmb), with=F] * enc_mat[idx, ]
        ) + 1
      )
      cprxy <- c(cprxy, rep(0L, 2^k - length(cprxy))) # zero-padding
      cprxy <- cprxy / sum(cprxy)
      
      pxy[[x + 1]][[y + 1]] <- nrow(dat[idx, c(cmb), with=F]) / nrow(dat)
      
      pr[[x + 1]][[y + 1]] <- cprxy
      
      A[[x + 1]][[y + 1]] <- A_xy(fi_xy = fi[[x + 1]][[y + 1]], k = k)
      
      if (any(A[[x + 1]][[y + 1]] < 0) | any(A[[x + 1]][[y + 1]] < 0)) 
        return(list(ATE = NA, cmb_diff = NA))
      
      
      pz[[x + 1]][[y + 1]] <- 
        # solve(A[[x + 1]][[y + 1]]) %*%  pr[[x + 1]][[y + 1]]
        vapply(
          seq_len(2^k),
          function(z) infer_pz(z, pr[[x + 1]][[y + 1]], fi[[x + 1]][[y + 1]], k),
          numeric(1L)
      )
      
      # assert_that(sum(abs(pz[[x + 1]][[y + 1]] - pz_cand)) < 10^(-10))
      
      A_inv[[x+1]][[y+1]] <- lapply(
        seq_len(2^k),
        function(z) infer_pz(z, pr[[x + 1]][[y + 1]], fi[[x + 1]][[y + 1]], k, T)
      )
      A_inv[[x+1]][[y+1]] <- do.call(rbind, A_inv[[x+1]][[y+1]])
      
      # max(abs(solve(A[[x + 1]][[y + 1]]) - A_inv)) (works!)
      
      # B represents P(z | r, x, y)
      B[[x + 1]][[y + 1]] <- invert_A_xy(A[[x + 1]][[y + 1]], 
                                         pz[[x + 1]][[y + 1]], 
                                         pr[[x + 1]][[y + 1]])
      
      
      if (check_simplex(pz[[x + 1]][[y + 1]])) {
        initial <- pz[[x + 1]][[y + 1]]
        pz[[x + 1]][[y + 1]] <- locate_simplex(pr[[x + 1]][[y + 1]], 
                                               A_inv[[x + 1]][[y + 1]])
        correction <- sum(abs(initial - pz[[x+1]][[y+1]]) * 100)
        if (verbose) {
          cat("Locate simplex correction ", round(correction, 2), "%\n", sep = "")
        }
      }
      
      if (min(pz[[x + 1]][[y + 1]]) < -(10)^(-10) | max(pz[[x + 1]][[y + 1]]) > 1) {
        cat("Walking out of simplex for (x,y)", x, y, "and range", 
            range(pz[[x + 1]][[y + 1]]), "\n")
      }
      
    }
    
  }
  
  pz_inf <- rep(0, 2^k)
  for (x in c(T, F)) {
    for (y in c(T, F)) {
      pz_inf <- pz_inf + pxy[[x + 1]][[y + 1]] * pz[[x + 1]][[y + 1]]
    }
  }
  
  dlt_x0x1 <- rep(0, 2^k)
  for (i in seq_len(2^k)) {
    
    
    ratio_x0 <- B[[x0 + 1]][[y1 + 1]][i, ] / 
      (B[[x0 + 1]][[y1 + 1]][i, ] * pyrx0 + B[[x0+1]][[y0+1]][i, ] * (1-pyrx0))
    ratio_x0[is.nan(ratio_x0)] <- 0
    
    ratio_x1 <- B[[x1+1]][[y1+1]][i, ] / 
      (B[[x1+1]][[y1+1]][i, ] * pyrx1 + B[[x1+1]][[y0+1]][i, ] * (1-pyrx1))
    ratio_x1[is.nan(ratio_x1)] <- 0
    
    # can the fi_xy case not be solved?! additional solving needed?
    przx0 <- A[[x0 + 1]][[y1 + 1]][, i] 
    przx1 <- A[[x1 + 1]][[y1 + 1]][, i]
    
    dlt_x0x1[i] <- sum(pyrx1 * ratio_x1 * przx1 - 
                         pyrx0 * ratio_x0 * przx0) # \sum_r {...}
  }

  wgh_vec <- sapply(seq_along(pz_inf) - 1, 
                    function(x) sum(as.integer(intToBits(x))))

  list(
    delta = dlt_x0x1,
    pz = pz_inf,
    ATE = sum(dlt_x0x1 * pz_inf) # \sum_z {...} P(z)
    #cmb_diff = sum(pzx1 * wgh_vec) - sum(pzx0 * wgh_vec)
  )
}
