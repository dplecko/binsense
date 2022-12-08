
expit <- function(x) exp(x) / (1 + exp(x))

synth <- synth_data()
dat <- synth[["dat"]]

fixy <- list(list(0, 0), list(0, 0))
fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- fi
sensitivity_ATE(fixy, dat, cmb)

res <- list()
p_range <- seq(0.0, fi, 0.05)
df <- expand.grid(fi = p_range)
for (j in seq_len(nrow(df))) {
  
  fixy <- list(list(0, 0), list(0, 0))
  fixy[[1]][[1]] <- fixy[[1]][[2]] <- fixy[[2]][[1]] <- fixy[[2]][[2]] <- df$fi[j]
  
  res[[j]] <- sensitivity_ATE(fixy, dat, cmb)
  df$adjusted_ATE[j] <- res[[j]][["ATE"]]
  
}

for (x in c(T, F)) {
  for (y in c(T, F)) {
    idx <- X == x & Y == y
    zsum <- table(Z[idx, ] %*% c(1, 2))
    cat(mean(idx), "*")
    cat(zsum / sum(zsum), "\n")
  }
}

sensitivity_ATE(fixy, dat, cmb)

ggplot(df, aes(x = fi, y = adjusted_ATE)) +
  geom_line() + theme_bw() + geom_point() +
  xlab(latex2exp::TeX("Fidelity parameter \\phi")) +
  ylab("Adjusted ATE") +
  ggtitle(latex2exp::TeX("Basic case: $\\phi = \\phi_{x} = \\phi_{xy}$")) +
  scale_y_continuous(labels = scales::percent)
