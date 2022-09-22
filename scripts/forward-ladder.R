
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE), source))

dat <- synth_data()[["dat"]]
fixy <- list(list(0, 0), list(0, 0))
fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3

naive_bayes(dat, fixy)
# conclusion: naive bayes fails when Z1,...,Zk are not independent.

dat <- synth_data()[["dat"]]
fixy <- list(list(0, 0), list(0, 0))
fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3
spcol_bayes(dat, fixy)

# k = 2 (manual) example
# dat <- synth_data()[["dat"]]
# fixy <- list(list(0, 0), list(0, 0))
# fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3
# model_bayes(dat, fixy)

# k = 5 (u-model) example
aux <- synth_data_mid()
dat <- aux[["dat"]]

fixy <- list(list(0, 0), list(0, 0))
fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0.3
fwd <- lapply(
  list(22, c(25, 26, 27), c(28, 29, 30, 31, 32)),
  function(sid) {
    model_bayes(copy(dat), fixy, pz_ground = aux[["pz"]], niter = 5L,
                Z_ground = aux$Z, seed = sid)
  }
)

# make plots for ATE changing with the iteration
df <- lapply(fwd, function(dt) {
  ret <- data.frame(beta = unlist(dt$beta), seed = paste(dt$seed, collapse = "_"))
  cbind(ret, update = seq_len(nrow(ret)))
})
df <- do.call(rbind, df)

ggplot(df, aes(x = update, y = -beta, color = factor(seed))) +
  geom_line() + geom_point() + theme_bw() + 
  geom_hline(yintercept = -aux[[5]], linetype = "dashed", color = "red") +
  xlab("# Updates") + ylab(latex2exp::TeX("$\\hat{\\beta}_{i}$")) +
  theme(
    axis.title.y = element_text(size = 20),
    legend.position = "bottom"
  ) + 
  scale_color_discrete(name = "Monte Carlo Seed")

# ggplot(df, aes(x = update, y = nZ)) +
#   geom_line() + geom_point() + theme_bw() + 
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
#   xlab("# Updates") + ylab(latex2exp::TeX("$n_Z$ difference")) +
#   theme(
#     axis.title.y = element_text(size = 20)
#   )
