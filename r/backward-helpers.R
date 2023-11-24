
trade_inv_prob <- function(pz, pr, Ainv, fi) {
  
  k <- length(fi)
  while(TRUE) {
    
    neg_idx <- which(pz < 0)
    if (length(neg_idx) == 0) browser() # return(fi)
    
    pz_neg_init <- pz[neg_idx]
    
    fi_grad <- cpp_fi_grad(Ainv, neg_idx, pr, fi)
    
    # expecting all fi_grad to have negative direction (!)
    fi_grad[fi_grad > 0] <- 0 # don't allow increases in fi
    fi_grad[fi == 0] <- 0 # vanishing gradient for fi equal to zero
    
    # apply Armijo-Goldstein to ensure step size not too large
    
    alpha <- min(min(fi / abs(fi_grad), na.rm = TRUE), 1)
    gain <- -Inf
    c <- 1/2
    cat("Init alpha =", alpha, ";")
    while( gain < c * alpha * sum(fi_grad^2) ) {
      
      alpha <- alpha / 2
      #' re-compute pz -- * but only for a select few, speeds up! *
      fi_cand <- fi + alpha * fi_grad 
      pz_neg <- vapply(
        neg_idx,
        function(z) infer_pz(z, pr, fi_cand, k), # update infer_pz computation
        numeric(1L)
      )
      
      gain <- sum(pz_neg) - sum(pz_neg_init)
      
      if (alpha < 10^(-10)) browser()
      
      if (gain > 0 & (gain < 1/5 * abs(sum(pz_neg_init)) | gain < 10^(-3))) break
      cat("Armijo-Goldstein sets alpha to", alpha, "\n")
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
    cat("Sum of negatives", sum(pz[pz < 0]), "\n")
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
    
    fi[[1]][[y + 1]] <- fi[[2]][[y + 1]] <- fi_new
  } else if (type == "y") {
    
    fi[[x + 1]][[1]] <- fi[[x + 1]][[2]] <- fi_new
  } else {
    
    fi[[x + 1]][[y + 1]] <- fi_new
  }
  
  fi
}

fi_prof <- function(fi) {
  
  k <- length(fi[[1]][[1]])
  fi_prof <- data.frame(
    rbind(
      cbind(fi[[1]][[1]], 0, 0, seq_len(k)),
      cbind(fi[[1]][[2]], 0, 1, seq_len(k)),
      cbind(fi[[2]][[1]], 1, 0, seq_len(k)),
      cbind(fi[[2]][[2]], 1, 1, seq_len(k))
    )
  )
  names(fi_prof) <- c("val", "x", "y", "feature")
  fi_prof$x <- paste("X =", fi_prof$x)
  fi_prof$y <- paste("Y =", fi_prof$y)
  ggplot(fi_prof, aes(x = feature, y = val)) +
    geom_col(position = "dodge") + theme_bw() +
    facet_grid(x ~ y) +
    ylab(latex2exp::TeX("$\\phi$ value")) + xlab("Feature Number")
}
