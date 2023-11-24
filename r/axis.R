
axis_search <- function(X, Y, R, fi_seq, type, solver, scale = "ATE", 
                        nboot = 1L) {
  
  df <- expand.grid(fi = fi_seq, boot_id = seq_len(nboot))
  if (is.list(X)) df$id <- as.numeric(factor(df$fi, levels = fi_seq))
  df$val <- 0
  fixy <- list(list(NULL, NULL), list(NULL, NULL))
  Sigma_list <- list()
  
  res <- mclapply(
    seq_len(nrow(df)),
    function(i) {
      
      if (is.list(X)) Xi <- X[[df$id[i]]] else Xi <- X
      if (is.list(Y)) Yi <- Y[[df$id[i]]] else Yi <- Y
      if (is.list(R)) Ri <- R[[df$id[i]]] else Ri <- R
      
      if (type == "x") {
        
        fixy[[1]][[1]] <- fixy[[1]][[2]] <- df$fi[i]
        fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0
      } else if (type == "y") {
        
        fixy[[1]][[1]] <- fixy[[2]][[1]] <- 0
        fixy[[1]][[2]] <- fixy[[2]][[2]] <- df$fi[i]
      } else {
        
        fixy[[1]][[1]] <- fixy[[1]][[2]] <- df$fi[i]
        fixy[[2]][[1]] <- fixy[[2]][[2]] <- df$fi[i]
      }
      
      if (df$boot_id[i] == 1) idx <- seq_along(Xi) else {
        
        idx <- sample(length(Xi), replace = TRUE)
      }
      
      binsensate(Xi[idx], Yi[idx], Ri[idx, ], fi = fixy, solver = solver)
    }, mc.cores = detectCores()
  )

  df$val <- vapply(res, `[[`, numeric(1L), scale)
  
  df <- as.data.table(df)
  if (nboot > 1L) {
    df <- df[, list(mid = val[boot_id == 1], 
                    lwr = val[boot_id == 1] - 1.96 * sd(val[boot_id != 1]),
                    upr = val[boot_id == 1] + 1.96 * sd(val[boot_id != 1])), 
             by = c("fi")]
  }
  
  list(df = df, res = res)
}
