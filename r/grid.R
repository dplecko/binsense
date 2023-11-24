
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
                  aes(x = fi_x0y0, y = 100 * ate)) + 
      geom_point() + geom_line() +
      theme_minimal() +
      xlab(latex2exp::TeX("$\\phi$")) + ylab("ATE") +
      scale_y_continuous(labels = scales::percent)
    
    return(agn)
  }
  
  if (type == "x") {
    
    trt <- ggplot(subset(lgrid, type == "x"), 
                  aes(x = fi_x0y0, y = fi_x1y0)) + 
      geom_tile(aes(fill = 100 * ate), color = 'white') +
      # scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
      theme_minimal() + 
      labs(fill="ATE") +
      scale_fill_viridis_c() +
      geom_text(aes(label=sprintf("%.2f%%", 100 * ate)), vjust=1) +
      xlab(latex2exp::TeX("$\\phi_{x_0}$")) + 
      ylab(latex2exp::TeX("$\\phi_{x_1}$"))
    
    return(trt)
  }
  
  if (type == "y") {
    
    out <- ggplot(subset(lgrid, type == "y"), 
                  aes(x = fi_x0y0, y = fi_x0y1)) + 
      geom_tile(aes(fill = 100 * ate), color = 'white') +
      # scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
      theme_minimal() + 
      labs(fill="ATE") + 
      scale_fill_viridis_c() +
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
