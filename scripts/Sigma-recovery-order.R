
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

kseq <- seq.int(3L, 7L, 1L)
nsampseq <- c(500, 1000, 2000, 5000, 10000, 20000, 50000)
famseq <- c("expfam-2d")
seedseq <- seq.int(35L, 44L)
solver <- c("expfam-2d-mom-grad")

df <- expand.grid(kseq, nsampseq, famseq, seedseq, solver, stringsAsFactors = FALSE)
names(df) <- c("dimension", "samples", "family", "seed", "solver")
df$Sigma_MSE <- df$ATE_MSE <- 0
i_range <- seq_len(nrow(df))
evl <- list()
for (i in i_range) {
  
  fixy <- list(list(0, 0), list(0, 0))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.1
  
  aux <- synth_data_mid(n = 100000, k = df$dimension[i], class = df$family[i],
                        fi = fixy, seed = df$seed[i])
  
  # recover \Sigma from R (marginally; ignore x, y conditioning)
  des_mom <- r_descent_2d(var(aux$R), aux$R, 0.1, n_iter = 100, multistep = TRUE)
  Sigma_hat <- tail(des_mom, 1L)[[1]]
  Sigma_err <- vec_norm(Sigma_hat - attr(aux$R, "Sigma"))^2
  
  # recover ATE from R, X, Y
  mcb <- microbenchmark::microbenchmark(
    sns <- binsensate(aux$X, aux$Y, aux$R, fixy, solver = df$solver[i]), 
    times = 1L
  )
  sns$runtime <- mcb$time
  df$Sigma_MSE[i] <- Sigma_err
  df$ATE_MSE[i] <- vec_norm(sns$ATE - aux$ATE)^2
  cat("Finished run", i, "\n")
}

# manipulate the samples order (^2, sqrt) and dimension (k, 1 / k)
ggplot(df, aes(x = log(samples), y = Sigma_MSE * samples / (dimension^2))) +
  geom_point() + theme_minimal() +
  facet_grid(cols = vars(dimension), scales = "free") +
  geom_smooth() + ylab("Scaled Error")

# the \Sigma MSE rate does not seem to be O(k^2 / n)
# is it because it is a method of moments?
# think about the transformation, it may have inverses etc.
# probably due to the additional error from Newton-Raphson?
# although the procedure seems quite stable; converges quickly


