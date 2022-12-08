
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

kseq <- seq.int(3L, 11L, 2L)
mcseq <- c(0, 1, 5)
famseq <- c("latent_u")
seedseq <- seq.int(22L, 24L)

df <- expand.grid(kseq, mcseq, famseq, seedseq, stringsAsFactors = FALSE)
names(df) <- c("dimension", "gibbs", "family", "seed")
df$MSE <- df$beta <- df$beta_hat <- 0

# i_range <- seq_len(nrow(df))
i_range <- 41L #' * go back to this, apparently interesting? *
evl <- list()
for (i in i_range) {
  
  aux <- synth_data_mid(n = 10^4, k = df$dimension[i], class = df$family[i],
                        seed = df$seed[i])
  dat <- aux[["dat"]]
  
  fixy <- list(list(0, 0), list(0, 0))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3
  
  fwd <- model_bayes(copy(dat), fixy, pz_ground = aux[["pz"]], niter = 20L,
                     Z_ground = aux$Z, gibbs_slope = df$gibbs[i])
  
  df$beta[i] <- aux$beta
  df$beta_hat[i] <- tail(fwd$beta, n = 1L)[[1]]
  df$MSE[i] <- abs(aux$beta - tail(fwd$beta, n = 1L)[[1]]) / 0.05
  evl[[i]] <- evl_frm(fwd, aux)
  
  cat("\r", i)
}

vis_evl(evl[[41L]])

bwd <- sensitivity_ATE(fixy, aux$dat)

library(plotly)
# volcano is a numeric matrix that ships with R
contour <- reshape2::acast(df, dimension ~ gibbs, value.var = "MSE", 
                           fun.aggregate = mean)
fig <- plot_ly(x = ~kseq, y = ~mcseq, z = ~contour)
fig <- fig %>% add_surface()
fig <- fig %>% layout(
  scene = list(xaxis=list(title = "Dimensionality"),
               yaxis=list(title = "# Gibbs slope"),
               zaxis=list(title = "relative L1 error")
  )
)

fig
