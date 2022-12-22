
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

kseq <- seq.int(2L, 8L, 6L)
famseq <- c("latent_u")
seedseq <- seq.int(22L, 22L)
niter <- 10L

df <- expand.grid(kseq, famseq, seedseq, stringsAsFactors = FALSE)
names(df) <- c("dimension", "family", "seed")
df$fMSE <- df$bMSE <- df$ate <- df$ate_fhat <-  df$ate_bhat <- df$pz_ferr <- df$pz_berr <- 0

i_range <- seq_len(nrow(df))
evl <- list()
for (i in i_range) { 
  
  aux <- synth_data_mid(n = 10^4, k = df$dimension[i], class = df$family[i],
                        seed = df$seed[i])
  dat <- aux[["dat"]]
  
  fixy <- list(list(0, 0), list(0, 0))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3
  
  bwd <- sensitivity_ATE(fixy, aux$dat)
  fwd <- model_bayes(copy(dat), fixy, pz_ground = aux$pz, niter = niter,
                     Z_ground = aux$Z, gibbs_slope = 1L,
                     tail_samples = 3L)
  
  df$ate[i] <- aux$ate
  df$ate_bhat[i] <- bwd$ATE
  df$ate_fhat[i] <- tail(fwd$ATE, n = 1L)[[1]]
  df$bMSE[i] <- abs(aux$ate - bwd$ATE)^2
  df$fMSE[i] <- abs(aux$ate - tail(fwd$ATE, n = 1L)[[1]])^2
  df$pz_berr[i] <- sum(abs(aux$pz - as.vector(bwd$pz)))
  df$pz_ferr[i] <- sum(abs(aux$pz - tail(fwd$pz, n = 1L)[[1]]))
  
  cat("\r", i)
}

# p1 <- ggplot(df, aes(x = dimension, y = bMSE)) +
#   geom_point() + theme_bw() + ylab("ATE error") + ggtitle("ATE error") +
#   geom_smooth()
# 
# p2 <- ggplot(df, aes(x = dimension, y = pz_err)) +
#   geom_point() + theme_bw() + ylab("P(z) error") + ggtitle("P(z) error") +
#   geom_smooth()
# 
# cowplot::plot_grid(p1, p2, ncol = 2L)

p1 <- ggplot(df, aes(x = ate_bhat, y = ate, color = dimension)) +
  geom_point() + theme_bw() + 
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  ggtitle("Backward vs. True ATE") +
  theme(legend.position = "bottom") + scale_color_viridis_b() + xlim(c(-0.2, 0.3))

p11 <- ggplot(df, aes(x = ate_fhat, y = ate, color = dimension)) +
  geom_point() + theme_bw() + 
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  ggtitle("Forward vs. True ATE") +
  theme(legend.position = "bottom") + scale_color_viridis_b() + xlim(c(-0.2, 0.3))

p2 <- ggplot(df, aes(x = bMSE, y = fMSE)) +
  geom_point() + theme_bw() + 
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  ggtitle("Forward vs. Backward MSE")

p3 <- ggplot(df, aes(x = pz_berr, y = pz_ferr)) +
  geom_point() + theme_bw() + 
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  ggtitle("Forward vs. Backward P(z) error")

cowplot::plot_grid(p1, p11, p2, p3, ncol = 2L)
