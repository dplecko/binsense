
#'  * can also test for more complex patterns * 
man_inv <- function(fi, k) {
  
  Ainv <- array(0, dim = c(2^k, 2^k))
  for (i in seq_len(2^k)) {
    
    i_bit <- to_bits(i-1, k)
    for (j in seq_len(2^k)) {
      
      j_bit <- to_bits(j-1, k)
      
      if (any(j_bit < i_bit)) next
      
      r0 <- sum(j_bit)
      rmz0 <- sum(j_bit - i_bit)
      
      Ainv[i, j] <- (-fi)^(rmz0) * (1-fi)^(-r0)
    }
  }
  
  Ainv
}

fi <- 0.2
k <- 7
max(abs(solve(cpp_A_xy(fi, k)) - man_inv(fi, k)))

# explore the agnostic case...

dat_lst <- real_data()
R <- dat_lst$R[, 1:5]
X <- dat_lst$X
Y <- dat_lst$Y

# vary fi \in [0, 0.2]

# compute for each x,y: (i) number of entries < 0; (ii) sum < 0; ()
backward_diag <- function(X, Y, R, fi, verbose = FALSE, ...) {
  
  k <- ncol(R)
  enc_mat <- t(replicate(nrow(R), 2^(seq_len(k) - 1L)))
  
  pr <- pz <- pxy <- list(list(0, 0), list(0, 0))
  
  x0 <- y0 <- 0
  x1 <- y1 <- 1
  
  df <- data.frame()
  
  for (x in c(T, F)) {
    
    for (y in c(T, F)) {
      
      idx <- X == x & Y == y
      
      cprxy <- tabulate(
        rowSums(
          R[idx, , drop = FALSE] * enc_mat[idx, , drop = FALSE]
        ) + 1
      )
      cprxy <- c(cprxy, rep(0L, 2^k - length(cprxy))) # zero-padding
      cprxy <- cprxy / sum(cprxy)
      
      pxy[[x + 1]][[y + 1]] <- nrow(R[idx, , drop = FALSE]) / nrow(R)
      pr[[x + 1]][[y + 1]] <- cprxy
      
      pz[[x + 1]][[y + 1]] <- 
        vapply(
          seq_len(2^k),
          function(z) infer_pz(z, pr[[x + 1]][[y + 1]], fi[[x + 1]][[y + 1]], k),
          numeric(1L)
        )
      
      sum(pr[[x + 1]][[y + 1]] == 0)
      
      df <- rbind(
        df, 
        data.frame(x = x, y = y, fi = fi[[x + 1]][[y + 1]], 
                   cneg = sum(pz[[x + 1]][[y + 1]] < 0), 
                   sneg = -sum(pz[[x + 1]][[y + 1]][pz[[x + 1]][[y + 1]] < 0]),
                   bzero = sum(pr[[x + 1]][[y + 1]] == 0))
      )
    }
  }
  
  df
}

fixy <- list(list(0, 0), list(0, 0))
res <- NULL
fi_seq <- seq(0, 0.2, 0.01)
fi_seq <- c(fi_seq[1], 0.001, fi_seq[-1])
for (fi in fi_seq) {
  
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- fi
  res <- rbind(
    res,
    backward_diag(X, Y, R, fixy)
  )
  cat("\r", fi)
}

ggplot(res, aes(x = fi, y = cneg)) +
  geom_point() + geom_line() +
  theme_bw() +
  facet_grid(x ~ y, labeller = label_both) # +
  # geom_hline(aes(yintercept = bzero), color = "red")
