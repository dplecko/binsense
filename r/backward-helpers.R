
#' @importFrom utils combn
infer_pz <- function(z, pr, phi, k, weights = FALSE) {
  
  zbit <- as.integer(intToBits(z-1)[seq.int(1, k, 1)])
  idx <- which(zbit == 0)
  
  if (length(idx) == 0) {
    r_geq_z <- list(integer(0L))
  } else if (length(idx) == 1) {
    r_geq_z <- list(integer(0L), idx)
  } else {
    
    r_geq_z <- unlist(lapply(0:length(idx),
                             combn, 
                             x = idx,
                             simplify = FALSE), 
                      recursive = FALSE)
  }
  
  if (length(phi) == 1) {
    
    r_wgh <- unlist(
      lapply(r_geq_z, function(r) (-1)^length(r) * phi^length(r) * 
               (1-phi)^(-length(r) - sum(zbit > 0)))
    )
  } else if (length(phi) == k) {
    
    r_wgh <- unlist(
      lapply(r_geq_z, function(r) (-1)^length(r) * 
               prod(phi[r]) / 
               prod(1 - phi[c(r, which(zbit > 0))]))
    )
  }
  
  r_idx <- unlist(lapply(r_geq_z, function(r) sum(2^((1:k) - 1)[r]) + z))
  
  if (weights) {
    wgh <- rep(0, 2^k)
    wgh[r_idx] <- r_wgh
    return(wgh)
  }
  
  sum(r_wgh * pr[r_idx])
}

trade_inv_prob_ZINF <- function(pz, pr, fi) {
  
  k <- as.integer(log(length(pz), 2))
  while(TRUE) {
    
    pz_neg_init <- pz[1]
    
    # compute the gradient for the zero-inflated (ZINF) setting
    fi_grad <- cpp_fi_grad_ZINF(pr, fi)
    
    # fix any undesirable gradient values
    fi_grad[fi_grad > 0] <- 0 # don't allow increases in fi
    fi_grad[fi == 0] <- 0 # vanishing gradient for fi equal to zero
    
    # apply Armijo-Goldstein to ensure step size not too large
    alpha <- min(min(fi / abs(fi_grad), na.rm = TRUE), 1)
    gain <- -Inf
    c <- 1/2
    while( gain < c * alpha * sum(fi_grad^2) ) {
      
      alpha <- alpha / 2
      fi_cand <- fi + alpha * fi_grad 
      
      # compute P(z) for the fi_cand
      
      pz_new <- cpp_Ainv_xy_zinf(fi_cand, k) %*% pr
      gain <- pz_new[1] - pz_neg_init
      
      if (alpha < 10^(-10)) break
      if (gain > 0 & (gain < 1/5 * abs(sum(pz_neg_init)) | gain < 10^(-3))) break
    }
    
    # update the fi that was selected
    fi <- fi + alpha * fi_grad
    fi[fi < 0.001] <- 0
    
    # update the pz
    pz <- cpp_Ainv_xy_zinf(fi, k) %*% pr
    
    # check if the stopping criterion is achieved
    if(abs(sum(pz[pz < 0])) < 10^(-5)) break
  }
  
  fi
}

trade_inv_prob_IF <- function(pz, pr, Ainv, fi) {
  
  k <- length(fi)
  while(TRUE) {
    
    neg_idx <- which(pz < 0)
    if (length(neg_idx) == 0) break
    
    pz_neg_init <- pz[neg_idx]
    
    fi_grad <- cpp_fi_grad_IF(Ainv, neg_idx, pr, fi)
    
    # expecting all fi_grad to have negative direction (!)
    fi_grad[fi_grad > 0] <- 0 # don't allow increases in fi
    fi_grad[fi == 0] <- 0 # vanishing gradient for fi equal to zero
    
    # apply Armijo-Goldstein to ensure step size not too large
    alpha <- min(min(fi / abs(fi_grad), na.rm = TRUE), 1)
    gain <- -Inf
    c <- 1/2
    while( gain < c * alpha * sum(fi_grad^2) ) {
      
      alpha <- alpha / 2
      fi_cand <- fi + alpha * fi_grad 
      pz_neg <- vapply(
        neg_idx,
        function(z) infer_pz(z, pr, fi_cand, k), # update infer_pz computation
        numeric(1L)
      )
      
      gain <- sum(pz_neg) - sum(pz_neg_init)
      
      if (alpha < 10^(-10)) break
      
      if (gain > 0 & (gain < 1/5 * abs(sum(pz_neg_init)) | gain < 10^(-3))) break
    }
    fi <- fi + alpha * fi_grad
    fi[fi < 0.001] <- 0
    
    # -- re-compute entire pz --
    pz <- vapply(
      seq_len(2^k),
      function(z) infer_pz(z, pr, fi, k),
      numeric(1L)
    )
    Ainv <- do.call(rbind, lapply(
      seq_len(2^k),
      function(z) infer_pz(z, pr, fi, k, TRUE)
    ))
    
    if(abs(sum(pz[pz < 0])) < 10^(-5)) break
  }
  
  fi
}

cmp_vec <- function(a, b) all(a == b)

fi_trade <- function(fi, fi_new, x, y) {
  
  type <- NULL
  if (!cmp_vec(fi[[1]][[1]], fi[[2]][[1]]) |
      !cmp_vec(fi[[1]][[2]], fi[[2]][[2]])) type <- paste0(type, "x")
  if (!cmp_vec(fi[[1]][[1]], fi[[1]][[2]]) |
      !cmp_vec(fi[[2]][[1]], fi[[2]][[2]])) type <- paste0(type, "y")
  
  if (is.null(type)) {
    
    fi <- list(list(fi_new, fi_new), list(fi_new, fi_new))
  } else if (type == "x") {
    
    fi[[x + 1]][[1]] <- fi[[x + 1]][[2]] <- fi_new
  } else if (type == "y") {
    
    fi[[1]][[y + 1]] <- fi[[2]][[y + 1]] <- fi_new
  } else {
    
    fi[[x + 1]][[y + 1]] <- fi_new
  }
  
  fi
}

check_simplex <- function(pzx) {
  
  if (min(pzx) < 0 || max(pzx) > 1) {
    return(TRUE)
  }
  FALSE
}

#' @importFrom stats pbinom
ill_pose_sig <- function(pz, pr, A_inv_xy, nsamp, sig_lvl = 0.01) {
  
  # compute lower/upper bound based on a finite sample
  cmp_lub <- function(phat, n, alpha, lwr) {
    
    tol <- 0.001
    low <- 0
    high <- 1
    while (high - low > tol) {
      mid <- (low + high) / 2
      
      if (pbinom(n * phat, size = n, prob = mid, lower.tail = !lwr) > alpha) {
        
        if (!lwr) low <- mid else high <- mid
      } else {
        
        if (!lwr) high <- mid else low <- mid
      }
    }

    return((low + high) / 2)
  }
  
  neg_idx <- which(pz < 0)
  signif_viol <- vapply(
    neg_idx,
    function(i) {
      
      coefs <- A_inv_xy[i, ]
      n_terms <- sum(coefs != 0)
      
      pr_best <- vapply(
        seq_along(coefs),
        function(j) {
          if (coefs[j] == 0) return(0)
          cmp_lub(phat = pr[j], n = nsamp, alpha = sig_lvl, lwr = coefs[j] < 0)
        }, numeric(1L)
      )
      
      sum(coefs * pr_best) < 0
    }, logical(1L)
  )
  
  if (any(signif_viol)) return(TRUE)
  return(FALSE)
}

