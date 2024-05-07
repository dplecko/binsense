
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


rng_to_df <- function(type, ...) {
  
  if (type == "agnostic") {
    
    df <- data.frame(type = rep(type, length(list(...)[[1]])))
    df$fi_x0y0 <- df$fi_x1y0 <- df$fi_x0y1 <- df$fi_x1y1 <- list(...)[[1]]
  } else if (type == "x") {
    
    fi_grid <- expand.grid(a = list(...)[[1]], b = list(...)[[2]])
    
    df <- data.frame(type = rep(type, nrow(fi_grid)))
    df$fi_x0y0 <- df$fi_x0y1 <- fi_grid$a
    df$fi_x1y0 <- df$fi_x1y1 <- fi_grid$b
  } else {
    
    fi_grid <- expand.grid(a = list(...)[[1]], b = list(...)[[2]])
    
    df <- data.frame(type = rep(type, nrow(fi_grid)))
    df$fi_x0y0 <- df$fi_x1y0 <- fi_grid$a
    df$fi_x0y1 <- df$fi_x1y1 <- fi_grid$b
  }
  
  df
}

grid_search <- function(fi_grid, solver) {
  
  grid <- cbind(fi_grid, solver = solver)
  grid$ate <- 0
  
  i_range <- seq_len(nrow(grid))
  fixy <- list(list(0, 0), list(0, 0))
  
  grid_runs <- mclapply(
    i_range,
    function(i) {
      
      fixy[[1]][[1]] <- grid$fi_x0y0[i]
      fixy[[1]][[2]] <- grid$fi_x0y1[i]
      fixy[[2]][[1]] <- grid$fi_x1y0[i]
      fixy[[2]][[2]] <- grid$fi_x1y1[i]
      sns <- binsensate(X, Y, R, fi = fixy, solver = grid$solver[i])
      sns
    }, mc.cores = detectCores()
  )
  
  grid$ate <- vapply(grid_runs, `[[`, numeric(1L), "ATE")
  
  attr(grid, "profs") <- lapply(grid_runs, `[[`, "fi")
  grid
}

grid_to_plt <- function(lgrid, slvr, type) {
  
  if (type == "agnostic") {
    
    agn <- ggplot(subset(lgrid, type == "agnostic"), 
                  aes(x = fi_x0y0, y = ate)) + 
      geom_point() + geom_line() +
      theme_bw() +
      xlab(latex2exp::TeX("$\\phi$")) + ylab("ATE") +
      scale_y_continuous(labels = scales::percent)
    
    return(agn)
  }
  
  if (type == "x") {
    
    trt <- ggplot(subset(lgrid, type == "x"), 
                  aes(x = fi_x0y0, y = fi_x1y0)) + 
      geom_tile(aes(fill = ate), color = 'white') +
      # scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
      theme_bw() + 
      labs(fill="ATE") +
      scale_fill_viridis_c(label = scales::label_percent()) +
      geom_text(aes(label=sprintf("%.2f%%", 100 * ate)), vjust=1) +
      xlab(latex2exp::TeX("$\\phi_{x_0}$")) + 
      ylab(latex2exp::TeX("$\\phi_{x_1}$"))
    
    return(trt)
  }
  
  if (type == "y") {
    
    out <- ggplot(subset(lgrid, type == "y"), 
                  aes(x = fi_x0y0, y = fi_x0y1)) + 
      geom_tile(aes(fill = ate), color = 'white') +
      # scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
      theme_bw() + 
      labs(fill="ATE") + 
      scale_fill_viridis_c(label = scales::label_percent()) +
      geom_text(aes(label=sprintf("%.2f%%", 100 * ate)), vjust=1) +
      xlab(latex2exp::TeX("$\\phi_{y_0}$")) + 
      ylab(latex2exp::TeX("$\\phi_{y_1}$"))
    
    return(out)
  }
}

grid_multi_plt <- function(grid, slvr) {

  cowplot::plot_grid(
    grid_to_plt(grid, slvr, "agnostic"), 
    grid_to_plt(grid, slvr, "x"), 
    grid_to_plt(grid, slvr, "y"),
    ncol = 3L
  )
}
