
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

kseq <- seq.int(3L, 11L, 4L)
mcseq <- c(0, 1)#, 5)
famseq <- c("latent_u")
seedseq <- 22L #seq.int(22L, 24L)
tail_samp <- c(1L, 10L)

df <- expand.grid(kseq, mcseq, famseq, seedseq, tail_samp, 
                  stringsAsFactors = FALSE)
evl <- list()
names(df) <- c("dimension", "gibbs", "family", "seed", "tail")
df$MSE <- df$beta <- df$beta_hat <- df$tail_fluc <- df$tail_error <- 0

niter <- 30L
i_target <- 8L # seq_len(nrow(df))
for (i in i_target) {
  
  aux <- synth_data_mid(n = 10^4, k = df$dimension[i], class = df$family[i],
                        seed = df$seed[i])
  dat <- aux[["dat"]]
  
  fixy <- list(list(0, 0), list(0, 0))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3
  
  fwd <- model_bayes(copy(dat), fixy, pz_ground = aux$pz, niter = niter,
                     Z_ground = aux$Z, gibbs_slope = df$gibbs[i],
                     tail_samples = df$tail[i])
  
  df$beta[i] <- aux$beta
  df$beta_hat[i] <- tail(fwd$beta, n = 1L)[[1]]
  df$MSE[i] <- abs(aux$beta - tail(fwd$beta, n = 1L)[[1]]) / 0.05
  df$tail_fluc[i] <- sum(abs(diff(unlist(tail(fwd$beta, n = niter - 10L)))))
  df$tail_error[i] <- sum(abs(aux$beta - unlist(tail(fwd$beta, n = niter - 10L))))
  evl[[i]] <- evl_frm(fwd, aux)
  
  cat("\r", i)
}

# cowplot::plot_grid(
#   vis_evl(evl[[8]]), vis_evl(evl[[44]]), ncol = 2L,
#   labels = c("01", "10")
# )

ggplot(df, aes(x = tail_fluc, y = tail_error, color = factor(tail))) +
  geom_point() + theme_bw() +
  scale_color_discrete(name = "Tail Samples")

vis_evl(evl[[8]])
# 
bwd <- sensitivity_ATE(fixy, aux$dat)

sum(100 * abs(as.vector(bwd$pz) - aux$pz)) # backward error
sum(100 * abs(as.vector(fwd$pz[[niter]]) - aux$pz)) # forward error


#' * Umap of P(z) * 
# library(umap)
# pzf <- do.call(rbind, fwd$pz)
# pzf <- rbind(pzf, aux$pz) 
# ump <- umap(pzf, n_neighbors = 3)
# 
# ump_lay <- data.frame(ump$layout)
# ump_lay <- cbind(ump_lay, niter = c(seq_len(nrow(ump_lay)-1), -1))
# 
# ggplot(ump_lay, aes(x = X1, y = X2, color = niter > 0)) +
#   geom_point() + geom_text(aes(label = niter), vjust = 1.5) +
#   theme_bw() + ggtitle("UMAP representation") + geom_path()

