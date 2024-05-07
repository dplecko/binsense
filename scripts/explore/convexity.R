library(plotly)
library(data.table)
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))
expit <- function(x) exp(x) / (1 + exp(x))

gen_data <- function(n, alpha, fix0, fix1) {
  
  icept <- -1
  
  # cat("fi(x0) =", fix0, "; fi(x1) =", fix1)
  
  z <- rbinom(n, 1, 0.5)
  x <- rbinom(n, 1, prob = expit(alpha * z + icept))
  mask <- rbinom(n, 1, prob = ifelse(x, 1 - fix1, 1 - fix0))
  r <-  z * mask
  
  data.frame(x, r)
}

est_alpha <- function(data) {
  
  # cat("; P(X = 1 | R = 1) =", mean(data$x[data$r == 1]), ";",
  #     "P(X = 1 | R = 0) =", mean(data$x[data$r == 0]))
  
  logreg <- glm(x ~ ., data = data, family = "binomial")
  # cat("; alpha =", logreg$coefficients["r"], "\n")
  logreg$coefficients["r"]
}

fi_range <- seq(0, 0.3, 0.01)
fi_grid <- expand.grid(fi0 = fi_range, fi1 = fi_range)
fi_grid$alpha <- 0

for (i in seq_len(nrow(fi_grid))) {
  
  data <- gen_data(10000, 0.5, fix0 = fi_grid$fi0[i], fix1 = fi_grid$fi1[i])
  fi_grid$alpha[i] <- est_alpha(data)
}

# plot_ly(data = fi_grid, x = ~fi0, y = ~fi1, z = ~alpha, type = "surface")

# Convert data to matrix form for surface plotting
alpha_matrix <- with(fi_grid, tapply(alpha, list(fi0, fi1), sum, default = NA))

# Create the surface plot
plot_ly(x = ~fi_range, y = ~fi_range, z = ~alpha_matrix, type = "surface") %>%
  layout(
    scene = list(
      xaxis = list(title = "fi0"),
      yaxis = list(title = "fi1"),
      zaxis = list(title = "alpha")
    )
  )

# does it get absorbed into the intercept?! cannot quite be...?
ggplot(fi_grid, aes(x = fi0, y = fi1, fill = alpha)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw()

fi_dt <- as.data.table(fi_grid)

delta <- 0.14
ggplot(
  fi_dt[fi0 <= (0.3 - delta)][fi1 == (fi0 + delta)], aes(x = fi0, y = alpha)
) +
  geom_line() + theme_bw()


### convexity testing

# custom backward solver
bwd_custom <- function(X, Y, R, fi, verbose = FALSE, g_truth = NULL, k = NULL, 
                       ...) {
  
  pr <- pz <- pxy <- A <- A_inv <- B <- list(list(0, 0), list(0, 0))
  
  x0 <- y0 <- 0
  x1 <- y1 <- 1
  
  for (x in c(T, F)) {
    
    for (y in c(T, F)) {
      
      if (!is.null(g_truth)) {
        
        pxy <- g_truth$pxy
        pr <- g_truth$pr
      } else {
        
        enc_mat <- t(replicate(nrow(R), 2^(seq_len(k) - 1L)))
        
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
      }
      
      A[[x + 1]][[y + 1]] <- cpp_A_xy(fi_xy = fi[[x + 1]][[y + 1]], k = k)
      
      if (any(A[[x + 1]][[y + 1]] < 0) | any(A[[x + 1]][[y + 1]] < 0)) 
        return(list(ATE = NA, cmb_diff = NA))
      
      
      pz[[x + 1]][[y + 1]] <- 
        vapply(
          seq_len(2^k),
          function(z) infer_pz(z, pr[[x + 1]][[y + 1]], fi[[x + 1]][[y + 1]], k),
          numeric(1L)
        )
    }
  }
  
  pz_inf <- rep(0, 2^k)
  for (x in c(T, F)) {
    for (y in c(T, F)) {
      pz_inf <- pz_inf + pxy[[x + 1]][[y + 1]] * pz[[x + 1]][[y + 1]]
    }
  }
  
  py_x1z <- pz[[x1+1]][[y1+1]] * pxy[[x1+1]][[y1+1]] / (
    pz[[x1+1]][[y1+1]] * pxy[[x1+1]][[y1+1]] + 
      pz[[x1+1]][[y0+1]] * pxy[[x1+1]][[y0+1]]
  )
  
  py_x0z <- pz[[x0+1]][[y1+1]] * pxy[[x0+1]][[y1+1]] / (
    pz[[x0+1]][[y1+1]] * pxy[[x0+1]][[y1+1]] + 
      pz[[x0+1]][[y0+1]] * pxy[[x0+1]][[y0+1]]
  )
  
  py_x1z[is.nan(py_x1z)] <- 0
  py_x0z[is.nan(py_x0z)] <- 0
  ate <- sum((py_x1z - py_x0z) * pz_inf)
  
  # return(
  #   list(pz = pz_inf, ATE = ate, solver = "backward-direct")
  # )
  ate
}

# Step 1 - generate p(r | x, y) for all x, y
mono_test <- function() {
  
  k <- 5
  Sigma <- gen_Sigma(k)
  lambda <- runif(k, -1, 1)
  mu <- runif(k, -1, 1)
  beta <- 0
  
  pxy <- pr <- fixy <- list(list(NULL, NULL), list(NULL, NULL))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.15
  
  pz <- cpp_scaling_2d(Sigma)
  pz <- pz / sum(pz)
  
  px_z <- expit(cpp_p_z(lambda))
  logity_z <- cpp_p_z(mu)
  
  for (x in c(0, 1)) {
    
    for (y in c(0, 1)) {
      
      if (x == 1) px <- px_z else px <- 1 - px_z
      
      py_z <- expit(logity_z + beta * x)
      if (y == 1) py <- py_z else py <- 1 - py_z
      
      pz_xy <- pz * px * py / sum(pz * px * py)
      
      pxy[[1 + x]][[1 + y]] <- sum(pz * px * py)
      
      A_xy <- cpp_A_xy(fi_xy = fixy[[x + 1]][[y + 1]], k = k)
      
      pr[[1 + x]][[1 + y]] <- as.vector(A_xy %*% pz_xy)
    }
  }
  
  g_truth <- list(pr = pr, pxy = pxy)
  
  # Step 4 - search the grid [0, 0.2] x [0, 0.2]
  fi_range <- seq(0, 0.2, 0.005)
  fi_grid <- expand.grid(fi_x0 = fi_range, fi_x1 = fi_range)
  fi_grid$ate <- 0
  
  for (i in seq_len(nrow(fi_grid))) {
    
    fixy[[1]][[1]] <- fixy[[1]][[2]] <- fi_grid$fi_x0[i]
    fixy[[2]][[1]] <- fixy[[2]][[2]] <- fi_grid$fi_x1[i]
    
    # fixy[[1]][[1]] <- fixy[[2]][[1]] <- fi_grid$fi_x0[i]
    # fixy[[1]][[2]] <- fixy[[2]][[2]] <- fi_grid$fi_x1[i]
    
    fi_grid$ate[i] <- bwd_custom(X=NULL, Y = NULL, R = NULL, fi = fixy,
                                 g_truth = g_truth, k = k)
    
    # cat("\r", i, "of", nrow(fi_grid))
  }
  
  # make a plot
  p1 <- ggplot(
    fi_grid, aes(x = fi_x0, y = fi_x1, color = ate)
  ) +
    geom_point(size = 3.5) + theme_bw() +
    scale_color_viridis_c() +
    ggtitle(paste("Sign =", sign(sum(lambda * mu))))
  
  p2 <- ggplot(
    fi_grid, aes(x = fi_x0, y = ate, color = fi_x1, group = fi_x1)
  ) +
    geom_line() + theme_bw() +
    scale_color_viridis_c() +
    ggtitle(paste("Sign =", sign(sum(lambda * mu))))
  
  cowplot::plot_grid(p1, p2)
}

mono_plts <- list()
for (i in 1:6) {
  
  mono_plts[[i]] <- mono_test() + theme(plot.background = 
                                        element_rect(color = "black"))
}

cowplot::plot_grid(plotlist = mono_plts, ncol = 3L)
