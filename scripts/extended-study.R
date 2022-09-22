
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

kseq <- seq.int(3L, 11L, 2L)
mcseq <- c(1, 3, 5, 10)
famseq <- c("latent_u")
seedseq <- seq.int(2022L, 2027L)

df <- expand.grid(kseq, mcseq, famseq, seedseq, stringsAsFactors = FALSE)
names(df) <- c("dimension", "MC", "family", "seed")
df$MSE <- df$beta <- df$beta_hat <- 0

for (i in seq_len(nrow(df))) {
  
  aux <- synth_data_mid(n = 10^4, k = df$dimension[i], class = df$family[i],
                        seed = df$seed[i])
  dat <- aux[["dat"]]
  
  fixy <- list(list(0, 0), list(0, 0))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3
  
  fwd <- model_bayes(copy(dat), fixy, pz_ground = aux[["pz"]], niter = 5L,
                     Z_ground = aux$Z, seed = seq_len(df$MC[i]))
  
  df$beta[i] <- aux$beta
  df$beta_hat[i] <- tail(fwd$beta, n = 1L)[[1]]
  df$MSE[i] <- abs(aux$beta - tail(fwd$beta, n = 1L)[[1]]) / 0.05 
  
  cat("\r", i)
}


library(plotly)
# volcano is a numeric matrix that ships with R
contour <- reshape2::acast(df, dimension ~ MC, value.var = "MSE", 
                           fun.aggregate = mean)
fig <- plot_ly(x = kseq, y = mcseq, z = ~contour)
fig <- fig %>% add_surface()
fig <- fig %>% layout(
  scene = list(xaxis=list(title = "Dimensionality"),
               yaxis=list(title = "# MC samples"),
               zaxis=list(title = "relative L1 error")
  )
)

fig
