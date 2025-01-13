
root <- rprojroot::find_root(rprojroot::has_dir(".git"))
r_dir <- file.path(root, "r-utils")
invisible(lapply(list.files(r_dir, full.names = TRUE, recursive = TRUE),
                 source))
options(scipen = 999)

# load the real data
dat_lst <- real_data()
k <- 11
samp <- TRUE
R <- dat_lst$R[samp, seq_len(k)]
X <- dat_lst$X[samp]
Y <- dat_lst$Y[samp]

# get the initial estimates of \Sigma, \lambda, \mu, \beta
Sigma <- binsensate:::infer_Sigma_IM(X, Y, R, list(list(0, 0), list(0, 0)))
fit_lmb <- function(X, Y, R) {
  xfit <- glm(X ~ ., data = data.frame(X = X, R = R), family = "binomial")
  lambda0 <- xfit$coefficients[-1]
  lambda_icept <- xfit$coefficients[1]
  yfit <- glm(Y ~ ., data = data.frame(Y = Y, X = X, R = R), 
              family = "binomial")
  r_coeff <- grepl("^R", names(yfit$coefficients))
  mu_icept <- yfit$coefficients[1]
  mu0 <- yfit$coefficients[r_coeff]
  beta0 <- yfit$coefficients["X"]
  list(lambda = lambda0, mu = mu0, beta = beta0, l_icept = lambda_icept, 
       m_icept = mu_icept)
}
th_s <- fit_lmb(X, Y, R)
th_s$Sigma <- Sigma

# compute the minimal eigenvalue of second order \nabla^2 Q()
nabla2Q <- Qhess(th_s)
lambda_min <- min(eigen(nabla2Q)$values)

# explore different random directions around the \Sigma_hat
fi <- 0.2

plan(multisession, workers = detectCores() - 1L)  
ret <- future_lapply(
  c(0.0001, 0.001, 0.01, 0.1, 0.5, 1),
  function(size) {
    
    res <- NULL
    for (i in 1:50) {
      
      th_t <- perturb_th(th_s, size = size)
      M_th <- M_theta(th_t, th_s, fi)
      
      gamma_fos <- 
        vec_norm(nabla_Q(M_th, th_t, th_s, fi) - nabla_Q(M_th, th_s, th_s, fi)) /
        vec_norm(em_params(th_t) - em_params(th_s))
      
      gamma_gs <- 
        vec_norm(nabla_Q(th_t, th_t, th_s, fi) - nabla_Q(th_t, th_s, th_s, fi)) /
        vec_norm(em_params(th_t) - em_params(th_s))
      
      res <- rbind(res, c(size, i, gamma_fos, gamma_gs))
    }
    res
  }, future.seed = TRUE
)
ret <- do.call(rbind, ret)
ret <- as.data.frame(ret)
names(ret) <- c("size", "rep", "FOS", "GS")

ret <- reshape2::melt(ret, id.vars = c("size", "rep"))

# compare the \lambda_min (\Sigma) vs. \gamma values obtained
custom_labeller <- function(variable, value) {
  TeX(paste0("$||\\delta|| = ", value))
}

library(latex2exp)
ret$delta <- paste("||δ|| =", ret$size)

ggplot(ret, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.4) +
  theme_bw() +
  facet_wrap(~ delta, nrow = 2) +
  geom_vline(xintercept = lambda_min, linetype = "dashed", color = "red") +
  xlim(0, NA) +
  scale_fill_discrete(name = "Condition") +
  theme(legend.position = "bottom") +
  xlab("Parameter Value") + ylab("Probability Density")

ggsave(file.path(root, "paper", "figures", "fos-gs.png"), 
       width = 8, height = 5)
