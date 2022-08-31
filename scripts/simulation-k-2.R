
expit <- function(x) exp(x) / (1 + exp(x))

n <- 100000
pz <- c(0.5, 0.2, 0.2, 0.1)
k <- 2

set.seed(2022)

Z <- sample(seq_len(2^k)-1, n, prob = pz, replace = TRUE)

lambda1 <- lambda2 <- 1
mu1 <- mu2 <- 1
beta <- 0.1
fi <- 0.3


Z <- lapply(Z, function(z) as.integer(intToBits(z)[1:k]))
Z <- do.call(rbind, Z)

X <- rbinom(n, 1, expit(Z %*% c(lambda1, lambda2) - 2))

Y <- rbinom(n, 1, expit(Z %*% c(mu1, mu2) + X * beta - 2))

ate <- mean(expit(Z %*% c(mu1, mu2) + beta - 2) - expit(Z %*% c(mu1, mu2) - 2))



R <- Z
for (i in seq_len(k)) {
  
  R[, i] <- Z[, i] * rbinom(n, 1, prob = 1 - fi)
  
}

dat <- data.frame(X, R, Y)
dat$X <- ifelse(dat$X == 1, 30, 20)
dat$Y <- ifelse(dat$Y == 1, TRUE, FALSE)
names(dat) <- c("bmi", "CHF", "Pulmonary", "death")
dat <- as.data.table(dat)


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
