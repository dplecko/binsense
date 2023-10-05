
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))


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

# load data
dat_lst <- real_data()
R <- dat_lst$R
X <- dat_lst$X
Y <- dat_lst$Y

fi_seq <- seq(0, 0.2, 0.05)
fi_grid <- rbind(
  rng_to_df("agnostic", fi_seq),
  rng_to_df("x", fi_seq, fi_seq),
  rng_to_df("y", fi_seq, fi_seq)
)

grid <- cbind(fi_grid, solver = "backward-direct")
# grid <- rbind(cbind(fi_grid, solver = "backward-direct"), 
#               cbind(fi_grid, solver = "two-stage-em"))
grid$ate <- 0

i_range <- seq_len(nrow(grid))
evl <- list()
fixy <- list(list(0, 0), list(0, 0))
invisible(
  foreach(i = i_range) %dopar% { 
    
    fixy[[1]][[1]] <- grid$fi_x0y0[i]
    fixy[[1]][[2]] <- grid$fi_x0y1[i]
    fixy[[2]][[1]] <- grid$fi_x1y0[i]
    fixy[[2]][[2]] <- grid$fi_x1y1[i]
    sns <- binsensate(X, Y, R, fi = fixy, solver = grid$solver[i])
    grid$ate[i] <- sns$ATE
  }
)


grid_to_plt <- function(grid, solver) {
  
  grid <- subset(grid, solver == solver)
  agn <- ggplot(subset(grid, type == "agnostic"), 
                aes(x = fi_x0y0, y = 100 * ate)) + 
    geom_point() + geom_line() +
    theme_minimal() +
    xlab(latex2exp::TeX("$\\phi$")) + ylab("ATE")
  
  trt <- ggplot(subset(grid, type == "x"), 
                aes(x = fi_x0y0, y = fi_x1y0)) + 
    geom_tile(aes(fill = 100 * ate), color = 'white') +
    scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
    theme_minimal() + 
    labs(fill="ATE") + 
    geom_text(aes(label=sprintf("%.2f", 100 * ate)), vjust=1) +
    xlab(latex2exp::TeX("$\\phi_{x_0}$")) + ylab(latex2exp::TeX("$\\phi_{x_1}$"))
  
  out <- ggplot(subset(grid, type == "y"), 
                aes(x = fi_x0y0, y = fi_x0y1)) + 
    geom_tile(aes(fill = 100 * ate), color = 'white') +
    scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
    theme_minimal() + 
    labs(fill="ATE") + 
    #scale_fill_viridis_c() +
    geom_text(aes(label=sprintf("%.2f", 100 * ate)), vjust=1) +
    xlab(latex2exp::TeX("$\\phi_{y_0}$")) + ylab(latex2exp::TeX("$\\phi_{y_1}$"))
  
  cowplot::plot_grid(agn + ggtitle(solver), trt + ggtitle(solver), 
                     out + ggtitle(solver), ncol = 3L)
}

cowplot::plot_grid(grid_to_plt(grid, "backward-direct"), 
                   grid_to_plt(grid, "two-stage-em"))
