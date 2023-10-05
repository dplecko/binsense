
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE), 
                 source))
cpp_dir <- file.path(root, "cpp")
invisible(lapply(list.files(cpp_dir, full.names = TRUE), sourceCpp))

# load data
dat_lst <- real_data()
R <- dat_lst$R
X <- dat_lst$X
Y <- dat_lst$Y

solver <- "backward-direct"
type <- "x"
scale <- "ATE" # or "OR"

# axis paths (x0, y1)
df <- data.frame(fi = c(seq(0, 0.2, 0.01)))
df$val <- 0
fixy <- list(list(NULL, NULL), list(NULL, NULL))
Sigma_list <- list()
invisible(
  for (i in seq_len(nrow(df))) {
  # foreach (i = seq_len(nrow(df))) %dopar% {
    
    if (type == "x") {
      
      fixy[[1]][[1]] <- fixy[[1]][[2]] <- df$fi[i]
      fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0
    } else {
      
      fixy[[1]][[1]] <- fixy[[2]][[1]] <- 0
      fixy[[1]][[2]] <- fixy[[2]][[2]] <- df$fi[i]
    }
    
    #' * TODO: replace this line with binsensate_boot(..., nboot = 100) * 
    sns <- binsensate(X, Y, R, fi = fixy, solver = solver)
    df$val[i] <- sns[[scale]]
  }
)


#' * TODO: add geom_ribon() for uncertainty * 
ggplot(df, aes(x = fi, y = val)) +
  geom_point() + geom_line() + theme_minimal() +
  xlab(latex2exp::TeX(
    paste0("$\\phi", ifelse(type == "x", "_{x_0}", "_{y_1}"), "$")
  )) +
  geom_hline(yintercept = ifelse(scale == "ATE", 0, 1), color = "red", 
             linetype = "dotdash", linewidth = 0.75)
  

df$smooth <- loess(val ~ fi, df)$fitted

ggplot(df, aes(x = fi, y = smooth)) +
  geom_point() + geom_line() + theme_minimal() +
  xlab(latex2exp::TeX("$\\phi$")) +
  geom_hline(yintercept = ifelse(scale == "ATE", 0, 1), color = "red", 
             linetype = "dotdash", linewidth = 0.75)

# infer the Sigma assuming \Phi = (0, 0, 0, 0)
Sigma_mim0 <- infer_Sigma(X, Y, R, fi = list(list(0, 0), list(0, 0)))

# infer the Z coefficients at X, Y
lambda <- glm(X ~ ., data = data.frame(X, R), family = "binomial")$coef
mu <- glm(Y ~ ., data = data.frame(Y, X, R), family = "binomial")$coef[-2]
mu_mod <- glm(Y ~ ., data = data.frame(Y, X, R), family = "binomial")

ate_to_or <- function(ate, mu_mod) {
  
  logits <- predict(mu_mod, data.frame(X = 0, R))
  
  a <- 0
  b <- 1
  while (b - a < 0.005) {
    
    beta <- (b+a) / 2
    ate_b <- mean(expit(logits + beta) - expit(logits))
    
    if (ate_b > ate) b <- beta else a <- beta
  }
  
  exp((b+a) / 2)
}

# walk along the trajectory
nrep <- 5
df_new <- as.data.table(expand.grid(fi = df$fi, seed = seq_len(nrep)))
df_new[, value := ifelse(scale == "ATE", 0, 1)]
df_new <- merge(df_new, df[, c("fi", "smooth")], by = "fi")
for (i in seq_len(nrow(df_new))) {
    
    # generate simulated data 
    dat_sim <- synth_data_mid(
      n = nrow(R), k = nrow(Sigma_mim0),
      seed = df_new$seed[i],
      fi = list(list(df_new$fi[i], df_new$fi[i]), list(0, 0)),
      class = "expfam-2d",
      Sigma = Sigma_mim0, lam = lambda[-1], mu = mu[-1],
      icept_x = lambda[1], icept_y = mu[1],
      beta = ifelse(scale == "ATE", ate_to_or(df_new$smooth[i], mu_mod), 
                    log(df_new$smooth[i])) #' * would need an ate_to_or function *
    )
    
    # perform inference
    if (type == "x") {
      
      fixy[[1]][[1]] <- fixy[[1]][[2]] <- df$fi[i]
      fixy[[2]][[1]] <- fixy[[2]][[2]] <- 0
    } else {
      
      fixy[[1]][[1]] <- fixy[[2]][[1]] <- 0
      fixy[[1]][[2]] <- fixy[[2]][[2]] <- df$fi[i]
    }
    
    sns <- binsensate(dat_sim$X, dat_sim$Y, dat_sim$R, fi = fixy, 
                      solver = solver)
    df_new$value[i] <- sns[[scale]]
    cat("\r", i, "out of", nrow(df_new))
}

# get the confidence interval
plot_dat <- df_new[, list(mid = mean(value), 
                          lwr = mean(value) - 1.96 * sd(value),
                          upr = mean(value) + 1.96 * sd(value)), 
                   by = c("fi", "smooth")]

# check the correctness
ggplot(plot_dat, aes(x = fi, y = mid)) +
  geom_line() + geom_point() +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              linewidth = 0, alpha = 0.2) +
  geom_hline(yintercept = 1, color = "red", linetype = "dotdash",
             linewidth = 0.75) +
  theme_minimal() +
  geom_point(data = df, aes(x = fi, y = smooth), color = "red") +
  geom_line(data = df, aes(x = fi, y = smooth), color = "red")
  theme(legend.position = "bottom")
