
r_idx <- function(rk, i, j, k) {
  
  rk_bit <- to_bits(rk - 1, k)
  vapply(
    seq_len(2^k),
    function(x) {
      
      x_bit <- to_bits(x-1, k)
      
      if (any(x_bit < rk_bit)) return(FALSE)
      
      if (x_bit[i] == 0 | x_bit[j] == 0) return(FALSE)
      
      return(TRUE)
    }, logical(1L)
  )
}

# computes P(r | z)
r_fi_hash <- function(rk, fi, k) {
  
  rk_bit <- to_bits(rk - 1, k)
  vapply(
    seq_len(2^k),
    function(z) {
      
      z_bit <- to_bits(z-1, k)
      
      if (any(z_bit < rk_bit)) return(0)
      
      fi ^ sum(z_bit - rk_bit) * (1 - fi) ^ sum(rk_bit == 1)
      
    }, numeric(1L)
  )
}

compute_grad_ij_old <- function(i, j, r, pz, fi, mode = "cpp") {
  
  k <- as.integer(log(length(pz), 2))
  
  rhash <- 1 + rowSums(r * t(replicate(nrow(r), 2^(seq_len(k) - 1))))
  
  k_grad_ij <- parallel::mclapply(
    rhash, function(rk) {
      
      if (mode == "cpp") {
        
        fi_wgh_hash <- cpp_fi_hash(rk, fi, k)
        idx <- cpp_idx(rk, i, j, k)
      } else {
        
        fi_wgh_hash <- r_fi_hash(rk, fi, k)
        idx <- r_idx(rk, i, j, k) 
      }
      
      sum(fi_wgh_hash[idx] * pz[idx]) / sum(fi_wgh_hash * pz)
    }, mc.cores = n_cores()
  )
  
  k_grad_ij <- unlist(k_grad_ij)
  
  if (i != j) k_grad_ij <- k_grad_ij * 2
  
  mean(k_grad_ij)
}

compute_grad_ij <- function(i, j, r, pz, fi, mc_cores = 1L) {
  
  k <- as.integer(log(length(pz), 2))
  
  rhash <- 1 + rowSums(r * t(replicate(nrow(r), 2^(seq_len(k) - 1))))
  
  # mc_cores <- 1 # n_cores()
  folds <- rep(seq_len(mc_cores), length(rhash) / 6 + 1)[seq_along(rhash)]
  
  if (mc_cores == 1) {
    
    k_grad_ij <- cpp_grad_ij(rhash, i, j, k, fi, pz)
  } else {
    
    k_grad_ij <- parallel::mclapply(
      seq_len(mc_cores), function(fold) {
        
        rhash_fold <- rhash[folds == fold]
        cpp_grad_ij(rhash_fold, i, j, k, fi, pz)
      }, mc.cores = mc_cores
    )
    
    k_grad_ij <- do.call(c, k_grad_ij)
  }
  
  if (i != j) k_grad_ij <- k_grad_ij * 2
  
  mean(k_grad_ij)
}

compute_grad_logc_ij <- function(i, j, pz) {
  
  dim <- as.integer(log(length(pz), 2))
  idx <- cpp_idx(1, i, j, dim)
  (1 + (i != j)) * 1 / sum(pz) * sum(pz[idx])
}

likelihood_2d <- function(Sigma, r, fi, mode = "cpp") {
  
  if (fi == 0) return(likelihood_2d_fi0(Sigma, r))
  pz <- cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  
  k <- as.integer(log(length(pz), 2))
  
  rhash <- 1 + rowSums(r * t(replicate(nrow(r), 2^(seq_len(k) - 1))))
  
  ll_k <- vapply(
    rhash, function(rk) {
      
      if (mode == "cpp") {
        
        fi_wgh_hash <- cpp_fi_hash(rk, fi, k)
      } else {
        
        fi_wgh_hash <- r_fi_hash(rk, fi, k)
      }
      
      log(sum(fi_wgh_hash * pz))
      
    }, numeric(1L)
  )
  
  mean(ll_k)
  
}

gradient_2d <- function(Sigma, r, fi) {
  
  if (fi == 0) return(gradient_2d_fi0(Sigma, r))
  pz <- cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  grad <- array(0, dim = dim(Sigma))
  k <- ncol(Sigma)
  
  for (i in seq_len(k)) {
    for (j in seq.int(i, k, 1)) {
      
      grad[i, j] <- compute_grad_ij(i, j, r, pz, fi) - 
        compute_grad_logc_ij(i, j, pz)
      
      grad[j, i] <- grad[i, j]
    }
  }
  
  grad
}

descent_2d <- function(Sigma0, r, fi, gamma = 1, n_iter = 50, true_Sigma = NULL,
                       verbose = FALSE, 
                       type = c("vanishing", "normalized", "backtracking", 
                                "RMSprop", "Rprop", "Newton-Raphson"),
                       eps = 10^(-log(nrow(r), 10) - 2)) {
  
  type <- match.arg(type, c("vanishing", "normalized", "backtracking", 
                            "RMSprop", "Rprop", "Newton-Raphson"))
  Sigma_dist <- loglc <- cor <- rep(0, n_iter)
  alpha <- gsc <- lam <- mu <- rep(0, n_iter)
  Sigma <- grads <- list()
  Sigma[[1]] <- Sigma0
  t_init <- 1
  beta <- 3 / 4
  k <- ncol(Sigma0)
  hess_sing <- FALSE
  
  loglc[1] <- likelihood_2d(Sigma[[1]], r, fi)
  if (!is.null(true_Sigma)) Sigma_dist[1] <- sum((Sigma[[1]] - true_Sigma)^2)
  gamma_i <- gamma
  
  if (ncol(Sigma0) == 2L & !is.null(true_Sigma)) {
    Sigma_mome <- k2_solver(r, fi)
  } else Sigma_mome <- NULL
  # cat("Iter:")
  for (i in seq.int(2, n_iter)) {
    
    # compute gradient
    curr_grad <- gradient_2d(Sigma[[i-1]], r, fi)
    grads[[i]] <- curr_grad
    
    if (!is.null(true_Sigma)) gsc[i] <- nabq_gamma(true_Sigma, Sigma[[i-1]], fi)
    
    # determine step size
    if (type == "vanishing") {
      
      gamma_i <- gamma / (1 + i)
    } else if (type == "normalized") {
      
      gamma_i <- gamma / (vec_norm(curr_grad) * (1 + i))
    } else if (type == "backtracking") {
     
      tuned <- FALSE
      
      Sigma_cand <- Sigma[[i-1]] + gamma_i * curr_grad
      Delta <- likelihood_2d(Sigma_cand, r, fi) - loglc[i - 1]
      
      if (Delta < 0) {
        
        while (Delta < 0) {
          
          gamma_i <- gamma_i / 2
          Sigma_cand <- Sigma[[i-1]] + gamma_i * curr_grad
          Delta <- likelihood_2d(Sigma_cand, r, fi) - loglc[i - 1]
        }
      } else if (Delta < gamma_i * 1 / 8 * vec_norm(curr_grad)^2) {
        
        while (Delta < gamma_i * 1 / 8 * vec_norm(curr_grad)^2) {
          
          gamma_i <- gamma_i * 2
          Sigma_cand <- Sigma[[i-1]] + gamma_i * curr_grad
          Delta <- likelihood_2d(Sigma_cand, r, fi) - loglc[i - 1]
          if (Delta < 0) break
        }
        if (Delta < 0) gamma_i <- gamma_i / 2
      } else gamma_i <- gamma_i # / 2
      
    } else if (type == "RMSprop") {
      
      if (i == 2) {
        
        vt <- vec_norm(curr_grad)^2
      } else {
        
        vt <- beta * vt + (1 - beta) * vec_norm(curr_grad)^2
      }
      gamma_i <- gamma / sqrt(vt)
    } else if (type == "Rprop") {
      
      if (i > 2) {
        
        if (cor(as.vector(curr_grad), as.vector(grads[[i-1]])) < 0) {
          
          gamma_i <- gamma_i / 2
        } else if (cor(as.vector(curr_grad), as.vector(grads[[i-1]])) > 0.5) {
          
          gamma_i <- gamma_i * 2
        }
      }
      
      Sigma_cand <- Sigma[[i-1]] + gamma_i * curr_grad
      Delta <- likelihood_2d(Sigma_cand, r, fi) - loglc[i - 1]
      
      if (Delta < 0) {
        
        while (Delta < 0) {
          
          gamma_i <- gamma_i * 1 / 2
          Sigma_cand <- Sigma[[i-1]] + gamma_i * curr_grad
          Delta <- likelihood_2d(Sigma_cand, r, fi) - loglc[i - 1]
          if (!is.numeric(Delta)) browser()
        }
      }
    } else if (type == "Newton-Raphson") {
      
      # browser()
      Sigma_comp <- Sigma[[i-1]][!upper.tri(Sigma[[i-1]])]
      pz <- cpp_scaling_2d(Sigma[[i-1]])
      # compute the Hessian matrix
      Hess <- cpp_hessian_2d(pz / sum(pz), which(!upper.tri(Sigma[[i-1]]))-1, 
                             ncol(Sigma[[i-1]]))
      
      update <- tryCatch(
        solve(Hess, curr_grad[!upper.tri(curr_grad)]),
        error = function(e) e
      )
      
      if (inherits(update, "error")) {
        
        if (!hess_sing) {
          cat("Hessian singular, moving to simple gradients\n")
          hess_sing <- TRUE
        } else cat(".")
        Sigma_cand <- Sigma[[i-1]] + gamma_i * curr_grad
        Delta <- likelihood_2d(Sigma_cand, r, fi) - loglc[i - 1]
        
        if (Delta < 0 | is.nan(Delta)) {
          
          gamma_i <- gamma / 2
          while (Delta < 0 | is.nan(Delta)) {
            
            gamma_i <- gamma_i * 1 / 2
            Sigma_cand <- Sigma[[i-1]] + gamma_i * curr_grad
            Delta <- likelihood_2d(Sigma_cand, r, fi) - loglc[i - 1]
            if (!is.numeric(Delta)) browser()
            if (gamma_i < 10^(-16)) {
              
              cat("No update in direction of gradient\n")
              break
            }
          }
        }

        Sigma[[i]] <- Sigma_cand
      } else {
        
        hess_sing <- FALSE
        if (!is.null(true_Sigma)) {
          
          evals <- eigen(Hess, only.values = TRUE)$values
          lam[i] <- min(evals)
          mu[i] <- max(evals)
        }
        
        halving <- 1
        Delta <- likelihood_2d(Sigma_pf(Sigma_comp, update, k), r, fi) - 
                 loglc[i - 1]
        while (Delta < 0 | is.nan(Delta)) {

          halving <- halving * 2
          Delta <- likelihood_2d(Sigma_pf(Sigma_comp, update/halving, k), r, fi) - 
                   loglc[i - 1]
        }
        Sigma[[i]] <- Sigma_pf(Sigma_comp, update/halving, k)
      }
    }
    
    if (!is.null(true_Sigma)) {
      
      a <- Sigma[[i-1]]
      b <- a + gamma_i * curr_grad
      c <- true_Sigma
      
      lambda_proj <- sum((c - a) * (b - a)) / vec_norm(b - a)
      cor[[i]] <- sum((c - a) * (b - a)) / (vec_norm(b - a) * vec_norm(c - a))
      alpha[[i]] <- lambda_proj / vec_norm(b - a)
    }
    
    if (type != "Newton-Raphson") Sigma[[i]] <- Sigma[[i-1]] + gamma_i * curr_grad
    
    loglc[i] <- likelihood_2d(Sigma[[i]], r, fi)
    
    Delta <- loglc[i] - loglc[i-1]
    
    if (!is.null(true_Sigma)) Sigma_dist[i] <- sum((Sigma[[i]] - true_Sigma)^2)
    
    if (i %% 10 == 0) {
      cat("\nLog(L_c) =", loglc[i], "at iter", i, "\n")
    }
    
    if (log(Delta, 10) < -10) break
  }
  
  if (!is.null(true_Sigma)) {
    loglc_true <- likelihood_2d(true_Sigma, r, fi)
    if (!is.null(Sigma_mome)) {
      mome_llc <- likelihood_2d(Sigma_mome, r, fi)
    } else mome_llc <- NA
  } else {
    loglc_true <- NA
    mome_llc <- NA
  }
  
  list(Sigma = Sigma[seq_len(i)], Sigma_dist = Sigma_dist[seq_len(i)], 
       loglc = loglc[seq_len(i)], loglc_true = loglc_true, 
       cor = cor[seq_len(i)], alpha = alpha[seq_len(i)],
       gsc = gsc[seq_len(i)], lam = lam[seq_len(i)], mu = mu[seq_len(i)],
       mome_llc = mome_llc)
}

vis_descent_2d <- function(object, iter_subset = TRUE) {
  
  df <- data.frame(Sigma_dist = object[["Sigma_dist"]], 
                   loglc = object[["loglc"]],
                   cor = object[["cor"]],
                   alpha = object[["alpha"]],
                   gsc = object[["gsc"]],
                   iter = seq_along(object[["Sigma"]]))
  df <- df[iter_subset, ]
  df <- reshape2::melt(df, id.vars = "iter")
  dummy2 <- data.frame(variable = c("Sigma_dist", "loglc", "cor", "alpha", 
                                    "loglc"), 
                       Z = c(0, object$loglc_true, 1, 1, object$mome_llc),
                       clr = c(1, 1, 1, 1, 2))

  ggplot(df, aes(x = iter, y = value)) +
    geom_point() + geom_line() + theme_bw() +
    facet_wrap(vars(variable), scales = "free_y") +
    geom_hline(data = dummy2, aes(yintercept = Z, color = factor(clr),
                                  linetype = factor(clr)), 
               linewidth=0.5)
}

k2_solver <- function(R, fi) {
  
  p11 <- mean(R[,1] == 1 & R[,2] == 1)
  p01 <- mean(R[,1] == 0 & R[,2] == 1)
  p10 <- mean(R[,1] == 1 & R[,2] == 0)
  p00 <- 1 - p11 - p01 - p10
  
  b2d <- -log(p10 / p11 * (1-fi) - fi)
  b2a <- -log(p01 / p11 * (1-fi) - fi)
  
  a2bd <- -log ( (1-fi)^2 * p00/p11 - fi * exp(- b2d) - fi * exp(- b2a) -
                   fi^2)
  
  a <- a2bd - b2d
  d <- a2bd - b2a
  b <- (a2bd - a - d) / 2
  
  matrix(c(a, b, b, d), ncol = 2L)
}

nabq_gamma <- function(Sigma_s, Sigma_c, fi) {
  
  Sigma_norm <- function(Sigma) vec_norm(Sigma[!lower.tri(Sigma)])
  
  k <- ncol(Sigma_s)
  rhash_all <- seq_len(2^k)
  
  pz_s <- cpp_scaling_2d(Sigma_s)
  pz_s <- pz_s / sum(pz_s)
  
  pz_c <- cpp_scaling_2d(Sigma_c)
  pz_c <- pz_c / sum(pz_c)
  
  A <- do.call(rbind, lapply(rhash_all, function(ri) cpp_fi_hash(ri, fi, k)))
  pr <- as.vector(A %*% pz_s)
  
  grad_diff <- Sigma_s
  grad_diff[,] <- 0
  for (j in seq.int(1, k)) {
    
    for (l in seq.int(j, k)) {
      
      ex_zrts <- cpp_grad_ij(rhash_all, j, l, k, fi, pz_s)
      ex_zrtc <- cpp_grad_ij(rhash_all, j, l, k, fi, pz_c)
      grad_diff[j, l] <- sum( (ex_zrts - ex_zrtc) * pr )
    }
  }
  
  Sigma_norm(grad_diff) / Sigma_norm(Sigma_s - Sigma_c)
}

vis_conv_radius <- function(object) {

  dat <- data.table(gam = object$gsc[-1], 
                    # mu = object$mu[-1], 
                    lam = object$lam[-1], 
                    iter = seq_along(object$gsc)[-1])
  dat <- melt(dat, id.vars = "iter")
  
  ggplot(dat, aes(x = iter, y = value, color = variable)) +
    geom_point() + geom_line() + theme_bw() +
    theme(legend.position = "bottom")
}
